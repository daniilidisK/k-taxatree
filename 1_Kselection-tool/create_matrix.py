# libraries
import pandas as pd
import os
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pylab as plt
import shutil
from sklearn.feature_extraction import DictVectorizer
from collections import Counter

# Specify k to create matrix
dna_bases_set = set('ACGT')
k = 6


def create_matrix(file,k):
    
    # open fasta file
    fastq_filehandle = open(file, "r")
    
    # Initialization
    counts = {}
    seq = 0

    # Loop over each line in the file
    line = ''
    seq_id="ID-{0}".format(seq)
    rownames = []
    
    # count kmers
    for row in fastq_filehandle:
        
        # Keep the rows with data
        if row[0] != '>':
          line += row.strip()
        
        else:
          rownames.append(row.strip()[1:])
         
          if seq != 0:
              counts =count_kmers(seq_id,line,k,counts)
          
          seq += 1
          seq_id="ID-{0}".format(seq)
          counts[seq_id]={}
          line = ''
    
    counts = count_kmers(seq_id,line,k,counts)
    
    # close fastq file
    fastq_filehandle.close

    restructured = [counts[key] for key in counts]
    
    dictvectorizer = DictVectorizer(sparse=False,  sort=False)
    md = dictvectorizer.fit_transform(restructured)
    
    #rownames = ['-'.join(['ID',str(i)]) for i in range(1,md.shape[0]+1)]
    
    md = pd.DataFrame(md, columns=dictvectorizer.get_feature_names(), index = rownames)

    return md

# Functions - 2
def count_kmers(id,read,k,counts):

    # Calculate how many kmers of length k there are	
    num_kmers = len(read) - k + 1
    kmers = [read[i:i+k] for i in range(num_kmers) if set(read[i:i+k]) <= dna_bases_set ]
    counts[id] = dict(Counter(kmers))
    
    # Return the final counts
    return counts


############# MAIN ############################
if __name__== '__main__':

    input_directory = 'input'

    # listing files in 'input directory
    files = os.listdir(input_directory)

    # fasta file
    fasta_file = next(x for x in files if x.endswith('.fasta'))
    fasta_file_path = input_directory + '/' + fasta_file
    
    print('Creating kmerMatrix...')
    kmerMatrix = create_matrix(fasta_file_path,int(k))
    
    plots_folder = 'output'
    file_name = fasta_file[:-6] + '_matrix_k_' + str(k) + '.csv'
    file_name = 'kmerMatrix.csv'
    file_name = plots_folder + '/' + file_name
    kmerMatrix.to_csv(file_name,  header = True, index = True)
    print('done')

