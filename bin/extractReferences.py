#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Filters results from taxonomic classifiers and extracts taxonomic IDs for filtered hits"""
"""To do: Add argparse and logging"""
import os
import sys
import subprocess
import pandas as pd

def main():
    # Input
    classifier = sys.argv[1]
    params = sys.argv[2]
    hits_report = sys.argv[3]
    ntc = None
    if len(sys.argv) > 4:
        ntc = sys.argv[4]
    reference_summary_refseq = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
    reference_summary_genbank = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
    #reference_summary_refseq = "rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
    #reference_summary_genbank = "rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

    # Get thresholds for filtering
    min_reads = ""
    min_reads_w_ntc = ""
    max_reads_ntc = ""
    min_ratio = ""
    min_reads_ratio = ""
    try:
        params = params.strip().split(",")
        if len(params) == 1:
            min_reads = int(params[0])
        elif len(params) == 4:
            min_reads_w_ntc = int(params[0])
            max_reads_ntc = int(params[1])
            min_ratio = float(params[2])
            min_reads_ratio = int(params[3])
        else:
            print("No filter parameters specified")
            sys.exit(1)
    except Exception as e:
        print("Failed parsing filtering options\n" + str(e))
        sys.exit(1)

    print("min_reads", min_reads)
    print("min_reads_w_ntc", min_reads_w_ntc)
    print("max_reads_ntc", max_reads_ntc)
    print("min_ratio", min_ratio)
    print("min_reads_ratio", min_reads_ratio)

    # Output
    summary_filtered = "{}_{}_filtered_taxonomic_hits.tsv".format(classifier,hits_report)
    summary_skipped = "{}_{}_skipped_taxonomic_hits.tsv".format(classifier,hits_report)

    # Parse hits
    all_hits = get_all_hits(hits_report, classifier)
    if ntc:
        ntc_hits = get_all_hits(ntc, classifier)
    # Filter hits
    filtered_hits = {}
    skipped_hits = {}
    for taxid, data in all_hits.items():
        reads = parse_reads(data)
        # Filter based on ntc sample
        if ntc:
            # Get number of reads for ntc
            if taxid in ntc_hits.keys():
                ntc_data = ntc_hits[taxid]
                ntc_reads = parse_reads(ntc_data)
            else:
                ntc_reads = 0
            # Filter based on difference in read counts between sample and ntc
            if (reads > min_reads_w_ntc) and (ntc_reads < max_reads_ntc):
                add_taxid(all_hits, taxid, filtered_hits, skipped_hits)
            elif ((ntc_reads / reads) < min_ratio) and (reads > min_reads_ratio):
                add_taxid(all_hits, taxid, filtered_hits, skipped_hits)
        # Filter based on number of reads
        else:
            if reads > min_reads:
                add_taxid(all_hits, taxid, filtered_hits, skipped_hits)
    sorted_hits = sorted(list(filtered_hits.items()), key=lambda k: k[1][0], reverse=True)
    print("Found {} species in the sample".format(len(filtered_hits)))
    for hit in filtered_hits:
      print (hit, filtered_hits[hit])

    # Download refseq assembly summary file from ncbi
    print("Downloading NCBI refseq reference summary file")
    download_ftp(reference_summary_refseq)

    # Get reference ids for filtered taxids
    taxid_accesion_map = {}
    # Search refseq for references
    missing_acc_refseq = parse_accession(reference_summary_refseq, filtered_hits.keys(), taxid_accesion_map)
    if len(missing_acc_refseq) > 0:
        # Download genbank assembly summary file from ncbi
        print("Downloading NCBI genbank reference summary file")
        download_ftp(reference_summary_genbank)
        # Search genbank for references that are not found in refseq
        missing_acc = parse_accession(reference_summary_genbank, missing_acc_refseq, taxid_accesion_map)
        if len(missing_acc) > 0:
            print("Failed to identify accession number of reference for taxons: {}".format(", ".join(missing_acc)))

    # Create summary files
    with open(summary_filtered, "w") as out:
        out.write("taxonomic_id\treads\tspecies\taccession_id\n")
        for taxid, info in sorted_hits:
            out.write("\t".join([taxid, info[0], info[1], taxid_accesion_map[str(taxid)]]) + "\n")
    with open(summary_skipped, "w") as out:
        out.write("taxonomic_id\tspecies\n")
        for taxid, name in skipped_hits.items():
            out.write(str(taxid) + "\t" + name + "\n")

def get_all_hits(results, classifier):
    """Parse all hits in classifier report for species level"""
    all_hits = {}
    with open(results, 'r') as f:
        res = f.readlines()[1:]
        for hit in res:
            hit = hit.strip().split("\t")
            if classifier == "kaiju":
                all_hits[hit[3]] = [hit[2], hit[4]]
            elif classifier in ["kraken2", "centrifuge"]:
                reads = hit[1]
                level = hit[3]
                name = hit[5].strip().rstrip()
                if level == "S":
                    taxid = hit[4]
                    all_hits[taxid] = [reads, name]
                elif level.startswith("S") or level == "-":
                    if name not in ["root","cellular organisms","Opisthokonta","Eumetazoa","Bilateria","Deuterostomia","Craniata","Vertebrata","Gnathostomata","Teleostomi","Euteleostomi","Sarcopterygii"," Eutheria","Dipnotetrapodomorpha","Tetrapoda","Amniota","Theria","Eutheria","Boreoeutheria","Euarchontoglires","Haplorrhini","Simiiformes","Catarrhini","Hominoidea","Homininae","Pseudomonas fluorescens group"]:
                         all_hits[taxid].append([reads, name, hit[4]])
            else:
                print("No classifier specified. Supported classifiers are kraken2, centrifuge and kaiju.")
                sys.exit(1)
    return all_hits

def parse_reads(data):
    reads = int(data[0])
    return reads

def add_taxid(all_hits, taxid, output, skipped):
    """Add taxonomic ID of species or one of the strains of the species to the output"""
    data = all_hits[taxid]
    # Check if substrain of taxon should be used
    if len(data) > 2:
        # Compare reads with substrains for species
        reads = parse_reads(data)
        for sub in data[2:]:
            if int(sub[0])/reads > 0.75:
                taxid = sub[2]
                data = sub
    name = data[1]
    # Filter out unwanted species
    if "phage" in name.lower() or taxid == "9606" or "group" in name.lower() or "unclassified" in name.lower():
        skipped[taxid] = name
        return
    #elif contaminant
    output[taxid] = [data[0], name]

def download_ftp(path):
    """Download file"""
    try:
        command = "wget {}".format(path)
        proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = proc.communicate()
        print("output", output)
        print("error", error)
    except Exception as e:
        print("Could not download file: {}\n".format(path) + e)
        sys.exit(1)

def parse_accession(referencefile, taxids, taxid_accesion_map):
    """Parses NCBI summary file and chooses one reference accession number for each species in input list, returns list of missing taxons"""
    missing_acc = []
    with open(os.path.basename(referencefile), "r") as f:
        references = pd.read_csv(f, sep='\t', skiprows=1, index_col=0, dtype=str)
        for taxid in taxids:
            refs = references.copy()
            refs = refs.loc[refs["taxid"] == taxid]
            if len(refs) == 1:
                acc = refs.index.tolist()[0]
            elif "reference genome" in list(refs["refseq_category"]):
                acc = refs.index[refs["refseq_category"] == "reference genome"].tolist()[0]
            elif "representative genome" in list(refs["refseq_category"]):
                acc = refs.index[refs["refseq_category"] == "representative genome"].tolist()[0]
            elif "Complete Genome" in list(refs["assembly_level"]):
                acc = refs.index[refs["assembly_level"] == "Complete Genome"].tolist()[0]
            elif len(refs) > 0:
                acc = refs.index.tolist()[0]
            else:
                acc = None
                missing_acc.append(taxid)
            # Add identified references
            taxid_accesion_map[taxid] = acc
    return missing_acc

main()

