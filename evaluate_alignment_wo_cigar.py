import time
import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt

start = time.time()

lreads_df = pd.read_csv('long_reads_alignment.csv')

random_sreads = pd.read_csv('short_reads_alignment.csv')

random_sreads = random_sreads.sample(n=100000)

histogram_x_limit = 15000

outlier_start_tf = []
outlier_end_tf = []
outlier_overlap_tf = []
outlier_read_type_start = []
outlier_read_type_end = []
outlier_read_type_overlap = []

outlier_start_len = []
outlier_end_len = []
outlier_overlap_len = []

normal_start_tf = []
normal_end_tf = []
normal_overlap_tf = []

overlap = {}

shoot = []

counter = 0

for i in random_sreads.index:
    temp = []
    x = random_sreads['qname'][i].split('/')
    lread_name = x[0]+'/'+x[1]+'/ccs'
    start_position = x[2].split('_')[1]
    read_type = x[2].split('_')[2]
    
    y = lreads_df.loc[lreads_df['qname'] == lread_name]
    

    if random_sreads['ref_start'][i] < 0 or np.isnan(random_sreads['ref_end'][i]):
        shoot.append(random_sreads['qname'][i])
        continue
    else:
        expected_start = y.ref_start.iat[0] + np.int(start_position)
        actual_start = random_sreads['ref_start'][i]
        
            
        expected_end = expected_start + np.int(random_sreads['qlength'][i])
        actual_end = random_sreads['ref_end'][i]
        
        start_deviation = abs(expected_start - actual_start)
        end_deviation = abs(expected_end - actual_end)
    
        percentage_overlap = abs(((random_sreads['qlength'][i] - start_deviation)/(random_sreads['qlength'][i])))*100
    
        ########################################### Plot data ##########################################
        
        if start_deviation > histogram_x_limit:
            temp.append(histogram_x_limit + 1)  # Set deviation to a number which is 1 greater than the short read length
        else:
            temp.append(start_deviation)
            
        if end_deviation > histogram_x_limit:
            temp.append(histogram_x_limit + 1) # Set deviation to a number which is 1 greater than the short read length
        else:
            temp.append(end_deviation)
        
        if percentage_overlap > 100:
            temp.append(0) # No overlap
        else:
            temp.append(percentage_overlap)
    
        overlap[random_sreads['qname'][i]] = temp
        
        ################################# Handling Outliers #####################################
        
        if start_deviation > random_sreads['qlength'][i]:
            outlier_start_tf.append(random_sreads['read_alignment'][i])
            outlier_read_type_start.append(read_type)
            outlier_start_len.append(start_deviation)
        else:
            normal_start_tf.append(random_sreads['read_alignment'][i])
        
        if end_deviation > random_sreads['qlength'][i]:
            outlier_end_tf.append(random_sreads['read_alignment'][i])
            outlier_read_type_end.append(read_type)
            outlier_end_len.append(end_deviation)
        else:
            normal_end_tf.append(random_sreads['read_alignment'][i])
        
        if percentage_overlap > 100:
            outlier_overlap_tf.append(random_sreads['read_alignment'][i])
            outlier_read_type_overlap.append(read_type)
        else:
            normal_overlap_tf.append(random_sreads['read_alignment'][i])
        
df = pd.DataFrame.from_dict(overlap, orient='index', columns=['start_deviation','end_deviation','percentage_overlap'])
axarr = df.hist(bins=100, color='#86bf91', grid=True, layout=(3,1), figsize=(15,15))


fig = axarr.flatten()

fig[0].set_xlabel('Deviation amount (in bp\'s)')
fig[0].set_ylabel('No. of reads')
fig[0].set_title('End deviation')

fig[1].set_xlabel('Percentage (%)')
fig[1].set_ylabel('No. of reads')
fig[1].set_title('Percentage Overlap')

fig[2].set_xlabel('Deviation amount (in bp\'s)')
fig[2].set_ylabel('No. of reads')
fig[2].set_title('Start deviation')

plt.tight_layout()
plt.savefig('histogram.pdf')

end = time.time()
print('It took: '+str(end-start)+' secs')