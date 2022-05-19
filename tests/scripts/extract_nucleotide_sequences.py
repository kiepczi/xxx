"""This scipt was used to extract nucleotide sequences for all SCO needed for backtranslation. 
"""

#SetUp
from Bio import SeqIO
from pathlib import Path
from Bio import SeqIO
from shutil import copy
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import pandas as pd


#Get a dictionary with OG as a key, and the genome_accession-protein_id as a value
datadir = Path("SCO_MSA")
filenames = sorted(datadir.glob("OG*"))

#Create a nested dictionary, with separate dictionaries for each OG, containning genome accession as a key and protein id as a value
SCO_genome_protein_ids = {}
for filename in filenames:
    if str(filename).split('.')[1] == 'fa':
        OG = str(filename).split('/')[-1].split('.')[0]
        SCO_genome_protein_ids[OG] = {}
        records = SeqIO.parse(filename, "fasta")
        for record in records:
            prot_seq = record.seq.replace('-', '')
            description = "GCF_"+record.description.split('_')[1]
            prot = "WP_"+record.description.split('_')[-1]
            SCO_genome_protein_ids[OG][description] = prot




def extract_cds_for_backtranslation(file, protein_id):
    """Extract CDS sequnces and their translated sequences for backtranslation for T-coffee. 
    """


    records = list(SeqIO.parse(file, "fasta")) #Load data
    #Extract the sequence of inetrest by checking if the protein_id in a description
    cds_of_interest = [record for record in records if protein_id in record.description]
    #Create new seq description
    genome_accession = "GCF_" + str(file).split('/')[1].split('_')[1]
    replicon = ''.join([record.description.split('|')[1].split('_cds')[0] for record in cds_of_interest])


    #Trim the sequence to exclude the stop codon 
    stop_codons = ["TAG", "TAA", "TGA"]
    no_stop_codon_seq = []
    translated_cds = []


    for record in cds_of_interest:
        replicon = record.description.split('|')[1].split('_cds')[0]
        locus_tag = record.description.split('[locus_tag=')[1].split(']')[0]
        new_description = genome_accession +' '+ replicon+' '+locus_tag
        if record.seq[-3:] in stop_codons:
            cds_record = SeqRecord(
                record.seq[:-3], 
                description = new_description,
                id=protein_id, #assign the protein id as sequence id
                name=new_description #assign the genome acession as the sequence description

            )
            no_stop_codon_seq.append(cds_record)

            prot_record = SeqRecord(
                record.seq[:-3].translate(),
                description = new_description,
                id = protein_id,
                name = new_description

            )
            translated_cds.append(prot_record)

        else:
            cds_record = SeqRecord(
                record.seq, 
                description = new_description,
                id=protein_id, #assign the protein id as sequence id
                name=new_description #assign the genome acession as the sequence description

            )
            no_stop_codon_seq.append(cds_record)

            prot_record = SeqRecord(
                record.seq.translate(),
                description = new_description,
                id = protein_id,
                name = new_description

            )
            translated_cds.append(prot_record)





    return no_stop_codon_seq, translated_cds




datadir = Path("CDS_files")
filenames = sorted(datadir.glob("GCF*"))

cds_per_OG = defaultdict(list)
prot_per_OG = defaultdict(list)

for filename in filenames:
    file = 'GCF_' + str(filename).split('/')[1].split('_')[1].split('.')[0]
    for k, v in SCO_genome_protein_ids.items():
        prot = SCO_genome_protein_ids[k][file]
        cds, translated = extract_cds_for_backtranslation(filename, prot)
        cds_per_OG[k].extend(cds)
        prot_per_OG[k].extend(translated)


for k, v in cds_per_OG.items():
    SeqIO.write(v, f'SCO_extracted_cds/{k}_cds.fasta', 'fasta')

for k, v in prot_per_OG.items():
    SeqIO.write(v, f'SCO_translated_cds/{k}_translated_cds.fasta', 'fasta')