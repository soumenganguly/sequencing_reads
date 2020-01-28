import numpy as np
from scipy.stats import uniform, norm
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import truncnorm

# Normal distribution in the range 90-110
def get_truncated_norm(mean=0, sd=1, low=0, upp=10):
    return truncnorm((low-mean)/sd, (upp-mean)/sd, loc=mean, scale=sd)

long_read_length = 12000
sum = 0
short_read_length = []
counter=0

while counter==0:
    data = get_truncated_norm(mean=95, sd=1, low=90, upp=110).rvs()
    sum +=int(data)
    if sum<long_read_length:
        short_read_length.append(int(data))
        continue
    else:
        sum=sum-int(data)
        counter=1

print(len(short_read_length))
print(long_read_length-sum)
print(short_read_length)

ax = sns.distplot(short_read_length,
                  bins=100,
                  kde=True,
                  color='skyblue',
                  hist_kws={"linewidth": 15,'alpha':1})
ax.set(xlabel='Distribution of short read lengths', ylabel='Frequency')
plt.show()