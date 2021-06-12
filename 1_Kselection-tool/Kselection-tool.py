# libraries
import pandas as pd
import os
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pylab as plt
import shutil
from sklearn.feature_extraction import DictVectorizer
from collections import Counter
from sklearn.metrics import pairwise_distances

############## PARAMETERS ####################
dna_bases_set = set('ACGT')

# specify k-values
kvals = list(range(4,16))

# specify taxonomic levels
tax_levels = ['kingdom', 'phylum', 'class', 'order']

############## FUNCTIONS ######################
def create_matrix(file,k):
    
    # open fasta file
    fastq_filehandle = open(file, "r")
    
    # Initialization
    counts = {}
    seq = 0

    # Loop over each line in the file
    line = ''
    seq_id="ID-{0}".format(seq)

    # count kmers
    for row in fastq_filehandle:
        
        # Keep the rows with data
        if row[0] != '>':
            line += row.strip()
        
        else:
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
    md = pd.DataFrame(md, columns=dictvectorizer.get_feature_names())

    return md

# Functions - 2
def count_kmers(id,read,k,counts):

    # Calculate how many kmers of length k there are	
    num_kmers = len(read) - k + 1
    kmers = [read[i:i+k] for i in range(num_kmers) if set(read[i:i+k]) <= dna_bases_set ]
    counts[id] = dict(Counter(kmers))
    
    # Return the final counts
    return counts

def exclude_unassigned(targets, kmerMatrix):
    
    to_keep = [index for index,value in enumerate(targets) if value != "Unclassified"]
    
    return to_keep


############# MAIN ############################
if __name__== '__main__':

    input_directory = 'input'

    # listing files in 'input directory
    files = os.listdir(input_directory)

    # fasta file
    fasta_file = next(x for x in files if x.endswith('.fasta'))
    fasta_file_path = input_directory + '/' + fasta_file

    # taxonomies files
    taxonomies_file = next(x for x in files if x.endswith('.csv') or x.endswith('.tsv'))
    taxonomies_file = input_directory + '/' + taxonomies_file

    # read taxonomy table
    taxonomy = pd.read_csv(taxonomies_file, header = 0, index_col = None)

    # initialization of silhouette coefficient dictionary
    silh_coef_dict = {}
    for level in tax_levels:
        silh_coef_dict[level] = []


    # calculate silhouette coefficient for multiple k-values
    for k in kvals:
        
        # printing 
        print('Analysis for k = {0}'.format(k))
        
        # constructing kmerMatrix
        kmerMatrix = create_matrix(fasta_file_path,k)
        print("kmermatrix created")
        
        
        # compute pairwise distances
        D = pairwise_distances(kmerMatrix, metric='hamming', n_jobs = 4)
        print("done")
        
        for level in tax_levels:
            
            # specify targets
            targets = taxonomy[level]
            
            # exclude unassigned
            to_keep = exclude_unassigned(targets, kmerMatrix)
            y = targets[to_keep]
            d = D[to_keep]
            d = d[:,to_keep]
            print("Done2")
            
            # appending silhouette score
            silh_coef_dict[level].append(silhouette_score(d, y, metric='precomputed', n_jobs = 4))
            print("Done3")
        
        del D,d
    
    print('\n')
    
    # Delete and recreate plots folder
    plots_folder = 'output'
    if os.path.isdir(plots_folder):
        shutil.rmtree(plots_folder)
    
    os.makedirs(plots_folder)
    
    file1 = open("out.txt","w")
    # find kvalues that give max silhouette 
    for level in tax_levels:
        
        silh_list = silh_coef_dict[level]
        #print(silh_list)

        # max value and index
        silh_max = max(silh_list)
        index_of_max = silh_list.index(silh_max)
        k_optimal = kvals[index_of_max]
        
        printing_statemet = 'Level: {0}, max_silhouette = {1}, k_value = {2} \n'.format(level, silh_max, k_optimal)
        file1.write(printing_statemet)
    
    
    # Printing statement
    full_output_directory = os.getcwd() + '/' + plots_folder + '/'    
    printing_statemet = 'Based on the printed information above and the plots generated inside {0} directory, please select k-value to generate kmerMatrix and run create_matrix.py script'.format(full_output_directory)
    file1.write('\n')
    file1.write(printing_statemet)
    file1.close()

    # Create the figures
    
    # Find max
    max_y = max([max(silh_coef_dict[level]) for level in tax_levels])

    fig = plt.figure(figsize=(15,10), tight_layout=True)
    plt.ylim(0, max_y + 0.1)
    plt.title('Silhouette coefficient: level: {0}'.format(level))
    plt.xlabel('K values')
    plt.ylabel('Score')

    # Construct the figures
    
    for level in tax_levels:

        plt.plot(kvals,silh_coef_dict[level])
        
    
    plt.legend(tax_levels)
    plt.savefig(plots_folder + '/silhouette.png', bbox_inches='tight', dpi=150)
    plt.close()
    

