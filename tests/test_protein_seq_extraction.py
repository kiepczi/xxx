from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path
from scripts.extract_protein_sequences import extract_SCO_protein_seq
import pytest


import zipfile
with zipfile.ZipFile("test_input/GCF_000092385.1_ASM9238v1_genomic.gbff.zip","r") as zip_ref:
    zip_ref.extractall("test_input")


"""During extraction of protein sequences from GenBank files, we found that most SCO identified by OrthoFinder
were found to as a single copy genes in their corresponding genomes, however in some cases there were multiple copies found.

The following 4 tests, were written to check if the extract_SCO_protein_seq function used to extract
protein sequences from GenBank files for a given protein id returns a correct protein sequence. All known sequences
used in these tests were downloaded directly from NCBI.  
"""



@pytest.fixture
def known_WP_003948652_1():
    """Known WP_003948652.1 sequence"""
    known = [record.seq for record in SeqIO.parse("fixtures/known_WP_003948652_1.fasta", "fasta")]
    return known


"""test_1_protein_extraction_from_GenBank - In this test, the WP_003948652.1 sequence for GCF_000092385.1 genome was 
extracted from the GenBank file using the extract_SCO_protein_seq function. In this case there should only be a single 
copy of WP_003948652.1 in GCF_000092385.1 genome. When coparing extracted protein sequence and known sequence, the test 
should pass. NOTE: That in here we compare a list of sequences to a list of sequences, and not strings!"""

def test_1_protein_extraction_from_GenBank(known_WP_003948652_1):
    """Test if the extract_SCO_protein_seq function extracts correct protein sequence from GenBank"""

    extracted_WP_003948652_1 = [record.seq for record in extract_SCO_protein_seq("test_input/GCF_000092385.1_ASM9238v1_genomic.gbff", "WP_003948652.1")]

    assert extracted_WP_003948652_1 == known_WP_003948652_1


@pytest.fixture
def known_WP_030387191_1():
    """Known WP_003948652.1 sequence"""
    known = [record.seq for record in  SeqIO.parse("fixtures/known_WP_030387191_1.fasta", "fasta")]
    return known


"""test_2_protein_extraction_from_GenBank - In this test, the WP_030387191.1 sequence for GCF_002028385.1 genome was 
extracted from the GenBank file using the extract_SCO_protein_seq function. In this instance, there should be 2 copies 
of WP_030387191.1 in GCF_002028385.1 genome. When comparing extracted protein sequences (list of 2 sequences) 
to a known sequence (list of a single sequence), the test should fail. """

def test_2_protein_extraction_from_GenBank(known_WP_030387191_1):
    """Test extract_SCO_protein_seq function. Test should fail. """

    extracted_WP_030387191_1 = [record.seq for record in extract_SCO_protein_seq("test_input/GCF_002028385.1_ASM202838v1_genomic.gbff", "WP_030387191.1")]    
    assert extracted_WP_030387191_1 == known_WP_030387191_1


"""test_3_protein_extraction_from_GenBank - From NCBI's IPGs (Identical Protein Groups), we known that there are more
than one instance of WP_030387191.1 in GCF_002028385.1, and these copies are identical. We can now check if the 
extract_SCO_protein_seq function, returns two identical sequences in this case. The test should pass. """

def test_3_protein_extraction_from_GenBank():
    """Test if the extract_SCO_protein_seq function extracts correct protein sequence from GenBank"""
    extracted_WP_030387191_1 = [record.seq for record in extract_SCO_protein_seq("test_input/GCF_002028385.1_ASM202838v1_genomic.gbff", "WP_030387191.1")]  

    assert extracted_WP_030387191_1[0] == extracted_WP_030387191_1[1]


"""test_4_protein_extraction_from_GenBank - In this test, we again compare the extracted WP_030387191.1 sequence form GCF_002028385.1 
genome. As, there are two copies of WP_030387191.1 in GCF_002028385.1 genome, and we compare list of sequences to list of sequences, 
we need the same number of sequences in both lists to find out if the correct sequences are extracted. In here, we duplicated the known 
sequence, and the test should pass as both copies of the extracted protein sequences are identical. """

def test_4_protein_extraction_from_GenBank(known_WP_030387191_1):
    """Test the extract_SCO_protein_seq function. Test should pass"""

    extracted_WP_030387191_1 = [record.seq for record in extract_SCO_protein_seq("test_input/GCF_002028385.1_ASM202838v1_genomic.gbff", "WP_030387191.1")] 
    assert extracted_WP_030387191_1 == known_WP_030387191_1 * 2

"""test_5_protein_extraction_from_GenBank = In this test, the WP_003948652_1 sequence for GCF_002028385.1 genome will be extracted,
and compared to a known known_WP_030387191_1 sequence. As the sequences are diffrent, the test should fail."""

def test_5_protein_extraction_from_GenBank(known_WP_030387191_1):
    """Test extract_SCO_protein_seq. Test should fail. """
    extracted_WP_003948652_1 = [record.seq for record in extract_SCO_protein_seq("test_input/GCF_000092385.1_ASM9238v1_genomic.gbff", "WP_003948652.1")]

    assert extracted_WP_003948652_1 == known_WP_030387191_1


