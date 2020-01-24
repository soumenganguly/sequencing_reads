# Create a random distribution, uniform distribution and normal distribution 
import numpy as np
from scipy.stats import uniform
import seaborn as sns
import matplotlib.pyplot as plt

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
plt.plot(list_of_uniform_random_numbers, uniform.pdf(list_of_uniform_random_numbers))
plt.show()