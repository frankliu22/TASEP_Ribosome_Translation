import pandas as pd
import numpy as np
import ast
import xlsxwriter
from codon_utils import *

def harvest_rho(bedgraph_file_path, tsv_file_path):
    """Given the file paths to the bedGraph and TSV files, for every gene in the human genome write into an
    Excel file the codon sequence for that gene, as well as the raw RefSeq values."""
    bedgraph = pd.read_csv(bedgraph_file_path, sep='\t', header = None, names = ['Chromosome', 'Start', 'End', 'Riboseq'])
    tsv = pd.read_csv(tsv_file_path, sep='\t', header = None, \
        names = ['Gene Name','Transcript Name','Chromosome','+/-',"5'-UTR","Exon Seq","3'-UTR","5'-UTR Coords","Coding Coords","3'-UTR Coords"])

    # Below is a dictionary of chromosome data frames.
    chromosome_dfs = dict()
    for i in range(1,23):
        chr_string = "chr" + str(i)
        chromosome_dfs[chr_string] = bedgraph.loc[bedgraph['Chromosome'] == chr_string]
    chromosome_dfs['chrX'] = bedgraph.loc[bedgraph['Chromosome'] == 'chrX']
    chromosome_dfs['chrY'] = bedgraph.loc[bedgraph['Chromosome'] == 'chrY']

    num_genes = len(tsv.index)

    codon_str_sequences = np.asarray(tsv['Exon Seq'].values)   # np array of all the codon sequences, each as a string
    codon_num_sequences = [None] * num_genes   # list of np arrays, where each array is numeric codon sequence
    gene_names = np.asarray(tsv['Gene Name'].values)   # np array of all the gene names, each as strings
    mapped_coordinates = [ast.literal_eval(s) for s in tsv['Coding Coords'].values]   # list of lists, each holding coding coordinates
    senses = np.asarray(tsv['+/-'].values)   # np array of +/- values, depending on which strand of DNA gene is on
    chromosomes = np.asarray(tsv['Chromosome'].values)   # np array of chromosome locations of the genes
    raw_rho_vectors = [None] * num_genes

    # First, get all the codon sequences
    for g in range(num_genes):
        codon_num_sequences[g] = nucleotide_seq_to_codon_nums(codon_str_sequences[g])

    # Next, assemble together all of the raw rho vectors
    for g in range(num_genes):
        if senses[g] == '-':
            continue   # skip those genes for now
        seq_length = len(codon_str_sequences[g])   # number of RNA nucleotides inside coding transcript
        raw_rho_vector = np.zeros(seq_length // 3)   # want integer number of indices
        contiguous_coordinates = []   # will hold a list, in order, of all the coordinates for the coding sequence
        coordinate_list = mapped_coordinates[g]
        for coordinate_str in coordinate_list:
            elems = coordinate_str.split(" ")
            begin, end = int(elems[2]), int(elems[3])
            contiguous_coordinates.extend(list(range(begin, end)))
        for i in range(seq_length // 3):
            raw_rho_vector[i] = np.mean([get_riboseq_val(chromosome_dfs[chromosomes[g]], contiguous_coordinates[3*i+c]) for c in range(0,3)])
        raw_rho_vectors[g] = raw_rho_vector

    # Write outputs to Excel file
    workbook = xlsxwriter.Workbook('Human Genome Input.xlsx')
    sheet = workbook.add_worksheet()
    for g in range(num_genes):
        if senses[g] == '-':
            continue   # skip negative genes for now
        row_index = 3 * g
        sheet.write(row_index, 0, gene_names[g])
        for column, value in enumerate(codon_num_sequences[g].tolist()):
            sheet.write(row_index+1, column, value)
        for column, value in enumerate(raw_rho_vectors[g].tolist()):
            sheet.write(row_index+2, column, value)
    workbook.close()

def get_riboseq_val(chromosome_df, coord):
    """Given the chromosome data frame and the integer coordinate, return the Riboseq
    value associated with that coordinate value, or 0 if it doesn't exist inside the dataframe."""
    rval = chromosome_df[(chromosome_df['Start'] <= coord) & (chromosome_df['End'] > coord)]
    if rval.empty:
        return 0
    return rval['Riboseq'].values[0]

def nucleotide_seq_to_codon_nums(nucleotide_seq):
    """Takes in a string consisting of the nucleotide sequence, and returns np array of codon numbers."""
    codon_seq_length = len(nucleotide_seq) // 3
    codon_num_list = np.zeros(codon_seq_length, dtype = int)
    for i in range(codon_seq_length):
        codon_num_list[i] = get_codon_num(nucleotide_seq[(3*i) : (3*i+3)])
    return codon_num_list

harvest_rho("GM18502.bedGraph", "test.tsv")
print("code reached completion")
