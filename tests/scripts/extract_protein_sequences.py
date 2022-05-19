#SetUp
import pandas as pd
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict




#Load data
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
            description = "GCF_"+record.description.split('_')[1]
            prot = "WP_"+record.description.split('_')[-1]
            SCO_genome_protein_ids[OG][description] = prot


def extract_SCO_protein_seq(file, protein_id):
    """Return protein sequnces from GenBank files with aived protein_id.
    """

    records = SeqIO.parse(file,"gb") #Load genbank
    accession = 'GCF_' + str(file).split('/')[-1].split('_')[1] #Get genome accession from the file name

    protein_records = []


    for record in records:
        for feature in record.features:
            for key, value in feature.qualifiers.items():
                if protein_id in value:
                    protein_seq = SeqRecord(
                        Seq(feature.qualifiers['translation'][0]), 
                        description = accession + ' '+ record.id +' '+feature.qualifiers['locus_tag'][0],
                        id=protein_id, #assign the protein id as sequence id
                        name=accession + ' '+ record.id +' '+feature.qualifiers['locus_tag'][0] #assign the genome acession as the sequence description
                    )
                    protein_records.append(protein_seq)



    return protein_records


all_OG_records = defaultdict(list)

datadir = Path("GenBank_files")
filenames = sorted(datadir.glob("GCF*"))

for filename in filenames:
    file = 'GCF_' + str(filename).split('/')[1].split('_')[1].split('.')[0]
    for k, v in SCO_genome_protein_ids.items():
        prot = SCO_genome_protein_ids[k][file]
        prot_record = extract_SCO_protein_seq(filename, prot)
        all_OG_records[k].extend(prot_record)



for k, v in all_OG_records.items():
    SeqIO.write(v, f'SCO_genbank_prot_seq/{k}_prot.fasta', 'fasta')