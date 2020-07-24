#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script takes an optional command-line argument which can be specified as the target location where the data should be downloaded and saved.
By default, all files are downloaded in the present working directory. 
"""
#import OS module to use os methods and functions
import os
import subprocess
import sys
import re
import pandas as pd
from Bio import SeqIO
import pdb
from Bio import Entrez
import csv

Entrez.email="sofia.stamouli@folkhalsomyndigheten.se"

name="accession_list.txt"

#set present working directory

cwd = sys.argv[1]
name = sys.argv[2]

if len(sys.argv) < 1:
    print ("Fatal: You forgot to include the directory name and the list with the contaminants")
    sys.exit(1)


#get the current working directory
os.chdir(cwd)
print(os.getcwd())
#function to process ftp url file that is created from assembly files
def process_url_file(inputurlfile):
    url_file=open(inputurlfile,'r')
    #The r means that the string is to be treated as a raw string, which means all escape codes to be ignored
    file_suffix=r'genomic.fna.gz'
    for line in url_file:
        url=line.rstrip('\n').split(',')
        ftp_url= url[0]+'/'+url[1]+'_'+url[2]+'_'+file_suffix
        print("Downloading"+ ftp_url)
        subprocess.call("wget "+ftp_url,shell=True)
        subprocess.call("gunzip *.gz",shell=True)
    subprocess.call("cat *.fna > contaminants.fna", shell=True)
    return

def download_contaminants(outfile='outfile_contaminants.txt'):
    assembly_summary_file="ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    #assembly_summary_file="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt"
    if os.path.exists('assembly_summary_refseq.txt'):
        os.remove('assembly_summary_refseq.txt')
    #Download the file using wget sysyem call
    subprocess.call("wget "+assembly_summary_file, shell=True)
    #Reformat the file to pandas-friendly format
    #Remove the first line of the file
    subprocess.call("sed -i '1d' assembly_summary_refseq.txt",shell=True)
    subprocess.call("sed -i 's/^# //' assembly_summary_refseq.txt", shell=True)
    #Read the file as a dataframe - using read_table
    #Use read_table if the column separator is tab
    df = pd.read_table('assembly_summary_refseq.txt',dtype='unicode')
    with open(name) as accession:
        accessionRows = [l.strip("\n") for l in accession]
    print(accessionRows)
    contaminants_list = df.organism_name.isin(accessionRows)
    print(contaminants_list)
    filtered_df = df[contaminants_list]
    my_df=filtered_df[['ftp_path','assembly_accession','asm_name']]
    my_df.to_csv(outfile,mode='w',index=False,header=None)
    process_url_file(outfile)
    return



download_contaminants('contaminants_complete_genome_url.txt')








#from Bio import Entrez


#Entrez.email="sofia.stamouli@folkhalsomyndigheten.se"

#name="accession_list.txt"

#def searchNCBI(identifier):

#        handle = Entrez.efetch(db='nucleotide', id=identifier, rettype='fasta')

#        name = identifier + ".fasta"
#        record = open(name, 'w')
#        record.write(handle.read())
#        handle.close()
#        record.close()


#with open(name, "r") as f:
#    for line in f:
#        ID = line.rstrip('\n')
#        print ("Getting fasta file for " +ID)
#        searchNCBI(ID)



