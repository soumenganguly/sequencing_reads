# Create a random distribution, uniform distribution and normal distribution 
import random
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns

# Random distribution
list_of_rand_numbers = []
for i in range(0,100):
    list_of_rand_numbers.append(random.randint(90,110))

plt.hist(list_of_rand_numbers, 21, density=1)
plt.xlabel('Random numbers')
plt.ylabel('Probability density')
plt.title('Probability distribution of random numbers in the range 90-110')
plt.show()

