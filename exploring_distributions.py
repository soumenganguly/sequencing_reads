# Create a random distribution, uniform distribution and normal distribution 
import numpy as np
from scipy.stats import uniform, norm
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import truncnorm

sns.set(rc={'figure.figsize':(10,10)})

# Random distribution
list_of_rand_numbers = []
for i in range(0,100):
    list_of_rand_numbers.append(np.random.randint(90,110))

ax = sns.distplot(list_of_rand_numbers,
                  bins=100,
                  kde=True)
ax.set(xlabel='Random numbers', ylabel='Probability density', title='Probability distribution of random numbers in the range 90-110')
plt.show()

# Uniform distribution
list_of_uniform_random_numbers = uniform.rvs(size=10000, loc = 90, scale=20)
ax = sns.distplot(list_of_uniform_random_numbers,
                  bins=100,
                  kde=True,
                  color='skyblue',
                  hist_kws={"linewidth": 15,'alpha':1})
ax.set(xlabel='Uniform Distribution ', ylabel='Frequency')
plt.show()

# Normal distribution
data_normal = norm.rvs(size=1000, loc=0, scale=1)
ax = sns.distplot(data_normal,
                  bins=100,
                  kde=True,
                  color='skyblue',
                  hist_kws={"linewidth": 15,'alpha':1})
ax.set(xlabel='Normal Distribution ', ylabel='Frequency')
plt.show()

# Normal distribution in the range 90-110

def get_truncated_norm(mean=0, sd=1, low=0, upp=10):
    return truncnorm((low-mean)/sd, (upp-mean)/sd, loc=mean, scale=sd)

data = get_truncated_norm(mean=95, sd=1, low=90, upp=110).rvs(1000)
ax = sns.distplot(data,
                  bins=100,
                  kde=True,
                  color='skyblue',
                  hist_kws={"linewidth": 15,'alpha':1})
ax.set(xlabel='Normal Distribution between 90-110', ylabel='Frequency')
plt.show()