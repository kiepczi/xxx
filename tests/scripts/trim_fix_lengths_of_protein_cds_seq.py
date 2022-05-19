#SetUp
from Bio.Seq import Seq
from Bio import SeqIO
from pathlib import Path

# data = "SCO_extracted_cds/OG0002136_cds.fasta"

def trim_cds_protein_sequences(file):
    """Return protein and cds sequences for each OG
    with correct lengths for backtranslation. 
    """

    records = list(SeqIO.parse(file, "fasta")) #Load data
    OG = str(file).split('/')[-1].split('_')[0]
    
    seen = []
    cds_sequences = []
    protein_sequences = []
    for cds_record in records:
        
        for protein_record in SeqIO.parse(f"test_input/SCO_genbank_prot_seq/{OG}_prot.fasta", "fasta"):
            if cds_record.description.split(' ')[1] in protein_record.description and cds_record.description.split(' ')[1] not in seen:
                seen.append(cds_record.description.split(' ')[1])
                if len(cds_record.seq) - len(protein_record.seq)*3 < 0:
                    prot_idx = len(cds_record.seq) - len(protein_record.seq)*3
                    protein_record.seq = cds_record.seq.translate()[:prot_idx]
                    protein_sequences.append(protein_record)
                    cds_idx = len(protein_record.seq)*3 - len(cds_record.seq)
                    cds_record.seq = cds_record.seq[:cds_idx]
                    cds_sequences.append(cds_record)
                elif len(cds_record.seq) - len(protein_record.seq)*3 >= 1:
                    cds_idx = len(protein_record.seq)*3 - len(cds_record.seq)
                    cds_record.seq = cds_record.seq[:cds_idx]
                    cds_sequences.append(cds_record)
                    protein_record.seq = cds_record.seq.translate()
                    protein_sequences.append(protein_record)
                else:
                    cds_sequences.append(cds_record)
                    protein_sequences.append(protein_record)

    # SeqIO.write(cds_sequences, f'SCO_cds_trimmed/{OG}_cds_trimmed.fasta', "fasta")
    # SeqIO.write(protein_sequences, f'SCO_protein_trimmed/{OG}_prot_trimmed.fasta', "fasta")

    return cds_sequences, protein_sequences

#Path
datadir_cds = Path("SCO_extracted_cds")
cds_filenames = sorted(datadir_cds.glob("OG*"))

for filename in cds_filenames:
    trim_cds_protein_sequences(filename)