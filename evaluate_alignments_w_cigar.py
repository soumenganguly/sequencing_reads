start = time.time()


short_samfile = pysam.AlignmentFile('aln-shortreads.sam', 'r' )
long_samfile = pysam.AlignmentFile('minimapalignment.sam', 'r' )

sreads=[]
for s in short_samfile.fetch():
    sreads.append(s)

lreads=[]
for l in long_samfile.fetch():
    lreads.append(l)


lreads_df = pd.read_csv('long_reads_alignment.csv')

random_sreads = pd.read_csv('short_reads_alignment.csv')

random_sreads = random_sreads.sample(n=1000)


outlier_start_tf = []
outlier_end_tf = []
outlier_overlap_tf = []
outlier_read_type = []

normal_start_tf = []
normal_end_tf = []
normal_overlap_tf = []

overlap = {}

shoot = []

counter = 0
reference_positions = [] # Store the reference positions of long reads for length 110 for overlap


for i in random_sreads.index:
    temp = []
    x = sreads[i].query_name.split('/')
    lread_name = x[0]+'/'+x[1]+'/ccs'
    start_position = x[2].split('_')[1]
    read_type = x[2].split('_')[2]
    
    y = lreads_df.loc[lreads_df['qname'] == lread_name]
    
    long_index = y.index[0]

    if sreads[i].reference_start < 0 or np.isnan(sreads[i].reference_start):
        shoot.append(sreads[i].query_name)
        continue
    else:
        ############################## Start deviation #####################################
        actual_start = sreads[i].reference_start
        expected_start_wo_cigar = lreads[long_index].reference_start + np.int(start_position)

        short_in_long_position = lreads[long_index].reference_start   # Initial long read reference start

        cigar_tuples = lreads[long_index].cigartuples # Cigar string for the long read

        for i in range(len(cigar_tuples)):
            if cigar_tuples[i][0] == 0 or cigar_tuples[i][0] == 2:   # Handle 'H' and 'S'
                counter+=cigar_tuples[i][1]
                if counter > int(start_position):
                    short_in_long_position = short_in_long_position + (cigar_tuples[i][1] - (counter - int(start_position)))
                    curr_cigar_length = cigar_tuples[i][1] - (cigar_tuples[i][1] - (counter - int(start_position))) # Clipped length of remaining operator
                    current_tuple_index = i # Store the index of the next cigartuple
                    break
                else:
                    short_in_long_position = short_in_long_position + cigar_tuples[i][1]
            elif cigar_tuples[i][0] == 1:
            # Skip for insertion
                continue
        
        expected_start_w_cigar = short_in_long_position
        
        start_deviation = abs(expected_start_w_cigar - actual_start)
        
        ############################### End deviation ############################################

        actual_end = sreads[i].reference_end
        expected_end_wo_cigar = lreads[long_index].reference_start + int(start_position) + sreads[i].query_length


        if curr_cigar_length <= sreads[i].query_length:
            traverse_length = short_in_long_position + curr_cigar_length
            reference_positions.extend(range(short_in_long_position, traverse_length))
            read_len_counter = curr_cigar_length
    
            for j in range(current_tuple_index+1,len(cigar_tuples)):
                if cigar_tuples[j][0] == 0 or cigar_tuples[j][0] == 2:   # Handle 'H' and 'S'
                    read_len_counter+=cigar_tuples[j][1]
                    if read_len_counter > int(sreads[i].query_length):
                        reference_positions.extend(range(traverse_length,traverse_length + (cigar_tuples[j][1] - (read_len_counter - int(sreads[i].query_length))) ))
                        traverse_length = traverse_length + (cigar_tuples[j][1] - (read_len_counter - int(sreads[i].query_length)))
                        break
                    else:
                        reference_positions.extend(range(traverse_length, traverse_length+cigar_tuples[j][1]))
                        traverse_length = traverse_length + cigar_tuples[j][1]
                elif cigar_tuples[j][0] == 1:
                    # Skip for insertion
                    continue       
        else:
            reference_positions.extend(range(short_in_long_position, short_in_long_position + int(sreads[i].query_length)))
            traverse_length = short_in_long_position + int(sreads[i].query_length)
        
        expected_end_w_cigar = traverse_length
        end_deviation = abs(expected_end_w_cigar - actual_end)
        
        
        ########################################## Overlap #########################################
        
        percentage_overlap = len(list(set(reference_positions) & set(sreads[i].get_reference_positions()))) * (sreads[i].query_length/100)
    
        ########################################## Plotting ########################################
        
        if start_deviation > sreads[i].query_length:
            temp.append(sreads[i].query_length + 1)  # Set deviation to a number which is 1 greater than the short read length
            outlier_start_tf.append(sreads[i].is_reverse)
            outlier_read_type.append(read_type)
        else:
            temp.append(start_deviation)
            normal_start_tf.append(sreads[i].is_reverse)
        
        if end_deviation > sreads[i].query_length:
            temp.append(sreads[i].query_length + 1) # Set deviation to a number which is 1 greater than the short read length
            outlier_end_tf.append(sreads[i].is_reverse)
        else:
            temp.append(end_deviation)
            normal_end_tf.append(sreads[i].is_reverse)
        
        if percentage_overlap > 100:
            temp.append(0) # No overlap
            outlier_overlap_tf.append(sreads[i].is_reverse)
        else:
            temp.append(percentage_overlap)
            normal_overlap_tf.append(sreads[i].is_reverse)
    
        overlap[sreads[i].query_name] = temp
        

df = pd.DataFrame.from_dict(overlap, orient='index', columns=['start_deviation','end_deviation','percentage_overlap'])
axarr = df.hist(bins=100, color='#86bf91', grid=True, layout=(3,1), figsize=(15,15))


fig = axarr.flatten()

fig[0].set_xlabel('Deviation amount')
fig[0].set_ylabel('No. of reads')
fig[0].set_title('End deviation')

fig[1].set_xlabel('Percentage (%)')
fig[1].set_ylabel('No. of reads')
fig[1].set_title('Percentage Overlap')

fig[2].set_xlabel('Deviation amount')
fig[2].set_ylabel('No. of reads')
fig[2].set_title('Start deviation')

plt.tight_layout()
plt.savefig('histogram.pdf')

end = time.time()
print('It took: '+str(end-start)+' secs')