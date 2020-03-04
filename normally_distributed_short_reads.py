import numpy as np
from scipy.stats import uniform, norm
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import truncnorm
import re
import json


fileLocation = open('data.fastq','r')
fileContents = fileLocation.read()

fileElements = fileContents.split('\n')

fileLocation.close()

allReadDict = {}

# Store the ID and reads(long) in a dictionary

oneReadDict = {}
for i in range(len(fileElements)):
    if re.match(r'^@', fileElements[i]):
        oneReadDict[fileElements[i]] = []
        if re.match(r'[A-Z]', fileElements[i+1]):
            oneReadDict[fileElements[i]].append(fileElements[i+1])
    else:
        continue

reads = oneReadDict.values()
ids = list(oneReadDict.keys())

# Sample short reads from normal distribution of read lengths

def get_truncated_norm(m=0, sd=1, low=0, upp=10):
    return truncnorm((low-m)/sd, (upp-m)/sd, loc=m, scale=sd)


def sample_short_reads(mean, standard_dev, sreadLen):
    temp = 0
    for r in reads:
        longReadLength = len(r[0])
        normDist = []
        readList = []
        counter = 0
        readCounter = 0
        boundary = 0

        while counter<=0:
            insertSize = get_truncated_norm(m=mean, sd=standard_dev, low=200, upp=300).rvs()
            normDist.append(insertSize)  # Normally distributed short read lengths
            if (readCounter + 2*sreadLen + int(insertSize)) <= longReadLength:
                r1 = r[0][readCounter:readCounter+sreadLen]    # The 0 is used to index the element stored as [[]] in r.
                r2 = r[0][readCounter + sreadLen + int(insertSize): readCounter + sreadLen + int(insertSize) + sreadLen]
                readList.append([r1+'_'+ids[temp]+'_'+str(readCounter),r2+'_'+ids[temp]+'_'+str(readCounter + sreadLen + int(insertSize))])
                boundary = readCounter + sreadLen + int(insertSize) + sreadLen
                readCounter += 1
                continue
            else:
                counter=1
        allReadDict[ids[temp]] = readList
        print(ids[temp])
        temp+=1
        print(len(readList))
        print(longReadLength, boundary)
        print(longReadLength-boundary)
    return allReadDict, normDist


readDict, ndist = sample_short_reads(250,1,95)

# PLot the distribution of short read lengths for one long read

ax = sns.distplot(ndist,
                  bins=100,
                  kde=True,
                  color='skyblue',
                  hist_kws={"linewidth": 15,'alpha':1})
ax.set(xlabel='Distribution of short read lengths', ylabel='Frequency')
plt.show()

# Store the contents into a json file

with open('normallydistributedreads.json','w') as jsonfile:
    json.dump(readDict,jsonfile)