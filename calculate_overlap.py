import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt



sreads_df = pd.read_csv('short_reads_alignment.csv')
lreads_df = pd.read_csv('long_reads_alignment.csv')

overlap = {}

start = time.time()

for count in range(100):
    temp = []
    x = sreads_df['qname'][count].split('/')
    lread_name = x[0]+'/'+x[1]+'/ccs'
    start_position = x[2].split('_')[1]
    
    y = lreads_df.loc[lreads_df['qname'] == lread_name]
    

    if sreads_df['ref_start'][count] < 0 or np.isnan(sreads_df['ref_end'][count]):
        shoot.append(sreads_df['qname'][count])
        continue
    else:
        expected_start = y.ref_start.iat[0] + np.int(start_position)
        actual_start = sreads_df['ref_start'][count]
        
            
        expected_end = expected_start + np.int(sreads_df['qlength'][count])
        actual_end = sreads_df['ref_end'][count]
        
        start_deviation = abs(expected_start - actual_start)
        end_deviation = abs(expected_end - actual_end)
    
        percentage_overlap = abs(((sreads_df['qlength'][count] - (start_deviation+end_deviation))/(sreads_df['qlength'][count])))*100
    
        if start_deviation > sreads_df['qlength'][count]:
            temp.append(sreads_df['qlength'][count] + 1)  # Set deviation to a number which is 1 greater than the short read length
        else:
            temp.append(start_deviation)
        
        if end_deviation > sreads_df['qlength'][count]:
            temp.append(sreads_df['qlength'][count] + 1) # Set deviation to a number which is 1 greater than the short read length
        else:
            temp.append(end_deviation)
        
        if percentage_overlap > 100:
            temp.append(0)                               # No overlap 
        else:
            temp.append(percentage_overlap)
    
        overlap[sreads_df['qname'][count]] = temp
        
df = pd.DataFrame.from_dict(overlap, orient='index', columns=['start_deviation','end_deviation','percentage_overlap'])
plt.hist(df['start_deviation'], bins=50, color='#86bf91')
plt.show()

# Add subplots, title, labels

end = time.time()
print('It took: '+str(end-start)+' secs')
