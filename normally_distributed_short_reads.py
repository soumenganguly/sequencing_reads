from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import numpy as np
from scipy.stats import uniform, norm
from scipy.stats import truncnorm
import re

sequences = []
read_counter = 0

all_sequences = list(SeqIO.parse('data.fastq', 'fastq'))


def get_truncated_norm(m=0, sd=1, low=0, upp=10):
    return truncnorm((low - m) / sd, (upp - m) / sd, loc=m, scale=sd)


def sample_short_reads(mean, standard_dev, sreadLen):
    temp = 0
    for record in all_sequences:
        longReadLength = len(record.seq)
        counter = 0
        readCounter = 0
        boundary = 0
        qual_score = record.letter_annotations['phred_quality']

        while counter <= 0:
            insertSize = get_truncated_norm(m=mean, sd=standard_dev, low=200, upp=300).rvs()
            if (readCounter + 2 * sreadLen + int(insertSize)) <= longReadLength:
                r1 = record.seq[
                     readCounter:readCounter + sreadLen]  # The 0 is used to index the element stored as [[]] in r.
                q1 = qual_score[readCounter:readCounter + sreadLen]
                idr1 = record.id + '_' + str(readCounter) + '_' + 'r1'
                sequences.append(SeqRecord(r1, id=idr1, letter_annotations={'phred_quality': q1}))

                r2 = record.seq[
                     readCounter + sreadLen + int(insertSize): readCounter + sreadLen + int(insertSize) + sreadLen]
                q2 = qual_score[
                     readCounter + sreadLen + int(insertSize): readCounter + sreadLen + int(insertSize) + sreadLen]
                idr2 = record.id + '_' + str(readCounter + sreadLen + int(insertSize)) + '_' + 'r2'
                sequences.append(SeqRecord(r2, id=idr2, letter_annotations={'phred_quality': q2}))
                boundary = readCounter + sreadLen + int(insertSize) + sreadLen
                readCounter += sreadLen
                continue
            else:
                counter = 1

        temp += 1
        print(temp)
        # print(longReadLength, boundary)
        # print(longReadLength-boundary)
    return sequences


s = sample_short_reads(250, 1, 110)

# Store the contents into a fastq file

with open('normallydistributedreads.fastq', 'w') as fastafile:
    SeqIO.write(s, fastafile, 'fastq-illumina')