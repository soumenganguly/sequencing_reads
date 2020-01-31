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

def get_truncated_norm(mean=0, sd=1, low=0, upp=10):
    return truncnorm((low-mean)/sd, (upp-mean)/sd, loc=mean, scale=sd)

temp = 0
for r in reads:
    longReadLength = len(r[0])
    normDist = []
    sum = 0
    normReads = {}
    counter = 0

    while counter == 0:
        normShortReadLen = get_truncated_norm(mean=95, sd=1, low=90, upp=110).rvs()
        normDist.append(normShortReadLen)
        if sum+int(normShortReadLen)<longReadLength:
            normReads[r[0][sum:sum+int(normShortReadLen)]] = sum
            sum+=int(normShortReadLen)
            continue
        else:
            counter=1
    allReadDict[ids[temp]] = normReads
    print(ids[temp])
    temp+=1
    print(len(normReads))
    print(longReadLength-sum)

# PLot the distribution of short read lengths for one long read

ax = sns.distplot(normDist,
                  bins=100,
                  kde=True,
                  color='skyblue',
                  hist_kws={"linewidth": 15,'alpha':1})
ax.set(xlabel='Distribution of short read lengths', ylabel='Frequency')
plt.show()

# Store the contents into a json file

with open('normallydistributedreads.json','w') as jsonfile:
    json.dump(allReadDict,jsonfile)