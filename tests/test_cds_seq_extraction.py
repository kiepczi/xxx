from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path
from scripts.extract_nucleotide_sequences import extract_cds_for_backtranslation
import pytest



"""OrthoFinder was used to find SCO, so the evolutionary history could be extimated
based on concatenated sequences. To build tree, the protein alignments had to be backtranslated. 

The following tests were written to test the extract_cds_for_backtranslation function, and check
if the correct sequences were retrived from FASTA files. 
"""

@pytest.fixture
def extracted_cds_WP_003948652_1():
    """Extract and return WP_003948652.1 cds sequence"""
    

    return extracted_cds_seq

@pytest.fixture
def known_cds_WP_003948652_1():
    """Known WP_003948652.1 cds sequence"""
    stop_codons = ['TAA', 'TAG', 'TGA']
    known = [record.seq[:-3] if record.seq[-3:] in stop_codons else record.seq for record in SeqIO.parse("fixtures/known_WP_003948652_1_cds.fasta", "fasta")]
    return known


"""test_1_cds_extraction - In this test, we again used the WP_003948652.1 sequence for GCF_000092385.1 genome (cds sequence) was 
extracted from the FASTA file using the extract_cds_for_backtranslation function. As we found in 
test_1_protein_extraction_from_GenBank, there was only 1 copy of that sequence in GenBank, therefore the test should pass. """

def test_1_cds_extraction(known_cds_WP_003948652_1):
    """Test extract_cds_for_backtranslation function. Test should pass. """
    extracted_cds, ectracted_tranlstaed_cds = extract_cds_for_backtranslation("test_input/GCF_000092385.1_ASM9238v1_cds_from_genomic.fna", "WP_003948652.1")
    extracted_WP_003948652_1_cds_seq = [record.seq for record in extracted_cds]

    assert extracted_WP_003948652_1_cds_seq == known_cds_WP_003948652_1


@pytest.fixture
def extracted_cds_WP_030387191_1():
    """Extract and return WP_030387191.1 cd sequence"""
    extracted_cds, extracted_translated_cds = extract_cds_for_backtranslation("test_input/GCF_002028385.1_ASM202838v1_cds_from_genomic.fna", "WP_030387191.1")
    extracted_cds_seq = [record.seq for record in extracted_cds]
    return extracted_cds_seq

@pytest.fixture
def known_cds_WP_030387191_1():
    """Known WP_030387191.1 cds sequence"""
    stop_codons = ['TAA', 'TAG', 'TGA']
    known = [record.seq[:-3] if record.seq[-3:] in stop_codons else record.seq for record in SeqIO.parse("fixtures/known_WP_030387191_1_cds.fasta", "fasta")]
    return known


"""test_2_cds_extraction - In this test, WP_030387191.1 cds sequence is extracted from GCF_002028385.1, and compared againsts a
known WP_030387191.1. This test should fail, as again the GCF_002028385.1 genome has two copies of WP_030387191.1, and 
we are comparing list of sequences to list of sequences!"""

def test_2_cds_extracted(known_cds_WP_030387191_1):
    """Test extract_cds_for_backtranslation function. Test should fail."""
    extracted_cds, extracted_translated_cds = extract_cds_for_backtranslation("test_input/GCF_002028385.1_ASM202838v1_cds_from_genomic.fna", "WP_030387191.1")
    extracted_cds_WP_030387191_1 = [record.seq for record in extracted_cds]
    assert known_cds_WP_030387191_1 == extracted_cds_WP_030387191_1


"""test_3_cds_extraction - From NCBI's IPGs (Identical Protein Groups), we known that there are more
than one instance of WP_030387191.1 in GCF_002028385.1, and these copies are identical. We can now check if the 
function indeed returns two identical sequences, and the test should pass. """

def test_3_cds_extracted():
    """Test if the extract_cds_for_backtranslation function returns two 
    identical copies of duplicated sequences. Test should pass."""
    extracted_cds, extracted_translated_cds = extract_cds_for_backtranslation("test_input/GCF_002028385.1_ASM202838v1_cds_from_genomic.fna", "WP_030387191.1")
    extracted_cds_WP_030387191_1 = [record.seq for record in extracted_cds]
    assert extracted_cds_WP_030387191_1[0] == extracted_cds_WP_030387191_1[1]


"""test_4_cds_extraction - In this test, we again compare the cds extracted WP_030387191.1 sequence form GCF_002028385.1 
genome. As, there are two copies of WP_030387191.1 in GCF_002028385.1 genome, and we compare list of sequences to list of sequences, 
we need the same number of sequences in both lists to find out if the correct sequences are extracted. In here, we duplicated the known 
sequence, and the test should pass as both copies of the extracted protein sequences are identical. """

def test_4_cds_extracted(known_cds_WP_030387191_1):
    """Test extract_cds_for_backtranslation funtion. Test should pass."""
    extracted_cds, extracted_translated_cds = extract_cds_for_backtranslation("test_input/GCF_002028385.1_ASM202838v1_cds_from_genomic.fna", "WP_030387191.1")
    extracted_cds_WP_030387191_1 = [record.seq for record in extracted_cds]

    assert extracted_cds_WP_030387191_1 == known_cds_WP_030387191_1 *2


"""test_5_cds_extraction - In this test, we compare two non-identical cds sequences, therefore the test should fail. """

def test_5_cds_extracted(known_cds_WP_030387191_1):
    """Test extract_cds_for_backtranslation funtion. Test should fail"""
    extracted_cds, ectracted_tranlstaed_cds = extract_cds_for_backtranslation("test_input/GCF_000092385.1_ASM9238v1_cds_from_genomic.fna", "WP_003948652.1")
    extracted_cds_WP_003948652_1  = [record.seq for record in extracted_cds]
    
    assert extracted_cds_WP_003948652_1 == known_cds_WP_030387191_1
