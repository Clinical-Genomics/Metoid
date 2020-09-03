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
#name = sys.argv[2]

if len(sys.argv) < 1:
    print ("Fatal: You forgot to include the directory name")
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
    return

#function to download bacterial sequences
def download_bacterial_genomes(outfile='outfile.txt'):
    assembly_summary_file=r'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
    if os.path.exists('assembly_summary.txt'):
       os.remove('assembly_summary.txt')
#    #Download the file using wget sysyem call
    subprocess.call("wget "+assembly_summary_file, shell=True)
#    #Reformat the file to pandas-friendly format
    subprocess.call("sed -i '1d' assembly_summary.txt",shell=True)
    subprocess.call("sed -i 's/^# //' assembly_summary.txt", shell=True)
    #Read the file as a dataframe - using read_table
    #Use read_table if the column separator is tab
    assembly_sum = pd.read_table('assembly_summary.txt',dtype='unicode')
    #filter the dataframe and save the URLs of the complete genomes in a new file
    my_df=assembly_sum[(assembly_sum['version_status'] == 'latest') &
                   (assembly_sum['assembly_level']=='Complete Genome') 
                  ]
    my_df=my_df[['ftp_path','assembly_accession','asm_name']]
    #output_file.write
    my_df.to_csv(outfile,mode='w',index=False,header=None)
    process_url_file(outfile)
    return
     
#function to download reference genomes 
#this function downloads latest version human reference genome by default 

def download_refseq_genome(taxid=9606,outfile='refseq_genome.txt'):
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
    assembly_sum = pd.read_table('assembly_summary_refseq.txt',dtype='unicode')
    #print('assembly_sum'+assembly_sum)
    my_df=assembly_sum[(assembly_sum['taxid'] == '9606') &
                       ((assembly_sum['refseq_category'] == 'reference genome') |
                        (assembly_sum['refseq_category'] == 'representative genome')
                       )]
    print('my_df'+my_df)
    my_df=my_df[['ftp_path','assembly_accession','asm_name']]
    #print('test'+my_df)
    #Process the newly created file and download genomes from NCBI website
    my_df.to_csv(outfile,mode='w',index=False,header=None)
    process_url_file(outfile)
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


#function to download viral sequences
def download_viral_genomes(outfile='outfile.txt'):
    #assembly_summary_file=r'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt'
    assembly_summary_file=r'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt'
    if os.path.exists('assembly_summary.txt'):
       os.remove('assembly_summary.txt')
#    #Download the file using wget sysyem call
    subprocess.call("wget "+assembly_summary_file, shell=True)
#    #Reformat the file to pandas-friendly format
    subprocess.call("sed -i '1d' assembly_summary.txt",shell=True)
    subprocess.call("sed -i 's/^# //' assembly_summary.txt", shell=True)
    #Read the file as a dataframe - using read_table
    #Use read_table if the column separator is tab
    assembly_sum = pd.read_table('assembly_summary.txt',dtype='unicode')
    #filter the dataframe and save the URLs of the complete genomes in a new file
    my_df=assembly_sum[(assembly_sum['version_status'] == 'latest') &
                   (assembly_sum['assembly_level']=='Complete Genome')
                  ]
    my_df=my_df[['ftp_path','assembly_accession','asm_name']]
    #output_file.write
    my_df.to_csv(outfile,mode='w',index=False,header=None)
    process_url_file(outfile)
    return


print('Downloading human genome'+'\n')
#change argument in the following function if you want to download other reference genomes
#taxonomy ID 9606 (human) should be replaced with taxonomy ID of genome of interest
download_refseq_genome(9606,'human_genome_url.txt')


print('Downloading bacterial genomes'+'\n')
download_bacterial_genomes('bacterial_complete_genome_url.txt')

download_viral_genomes('viral_complete_genome_url.txt')

#download_contaminants('contaminants_complete_genome_url.txt')

#Set name for the krakendb directory
#os.chdir("/srv/rs6/sofia/Metoid/Metoid/Test")
krakendb='HumanViralBacteria'
#kraken_dir="/srv/rs6/sofia/Metoid/Metoid/results/databases"
#human_db="Human"
#human_path=os.path.join(kraken_dir, human_db)
#os.mkdir(human_path)
subprocess.call('kraken2-build --download-taxonomy --use-ftp --db '+krakendb+' --threads 4', shell=True)
print('Running Kraken DB build for '+krakendb+'\n')
print('This might take a while '+'\n')
#for fna_file in os.listdir(cwd):
for fna_file in os.listdir(os.curdir):
    print ('fna ' + fna_file)
    #For all the databases
    #if fna_file.endswith('.fna') or fna_file.endswith('.fasta'):
    if fna_file.endswith('.fna'):
        print (fna_file)
        subprocess.call('kraken2-build --add-to-library '+fna_file +' --db '+krakendb+ ' --threads 4',shell=True)
subprocess.call('kraken2-build --build --db '+krakendb+' --threads 4' ,shell=True)

