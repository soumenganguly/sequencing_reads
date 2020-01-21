import re
import matplotlib.pyplot as plt

# Read the Fasta file containing long reads
file = open('dataset/m54238_180901_011437.Q20.fastq' , 'r')
file_contents = file.read()

elements = file_contents.split('\n')

file.close()

# Store the ID and Reads in a dictionary
id_read_dict = {}
for i in range(len(elements)):
    if re.match(r'^@', elements[i]):
        id_read_dict[elements[i]] = []
        if re.match(r'[A-Z]', elements[i+1]):
            id_read_dict[elements[i]].append(elements[i+1])
    else:
        continue

# Print the ID's along with reads
#for i in id_read_dict:
#    print(i, id_read_dict[i])

# Plot the distribution of the long read lengths in the fasta file
read_length_list = []
for i in id_read_dict:
    if len(id_read_dict[i]) == 0:
        read_length_list.append(0)
    else:
        read_length_list.append(len(id_read_dict[i][0]))

read_length_counter = {}
for i in read_length_list:
    read_length_counter[i] = read_length_list.count(i)

plt.plot(list(read_length_counter.keys()), list(read_length_counter.values()))
plt.title('Distribution of sequence lenghts over all sequences')
plt.xlabel('Length of long reads')
plt.ylabel('Frequency')
plt.show()


