from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path
from scripts.trim_fix_lengths_of_protein_cds_seq import trim_cds_protein_sequences
import pytest


import zipfile
with zipfile.ZipFile("test_input/GCF_000092385.1_ASM9238v1_genomic.gbff.zip","r") as zip_ref:
    zip_ref.extractall("test_input")




"""T-coffee was used for backtranslation of the protein alignments. However, some cds sequences lengths 
did not match the corresponding protein sequences x 3. There were 3 behaviours that were spotted druning the analysis:
Behaviour 1 - the length of cds sequence matched the length of protein sequence x 3. 
Behaviour 2 - the cds sequences were longer than the length of the correspoding sequence x 3.
Behaviour 3 - the protein sequence x 3 was longer than the length ofthe correspoding cds sequence.

To solve this problem, a trimming function was written to deal with the problem appropriately. 

The following test were written to test if the sequences are appropriately trimmed:
"""
@pytest.fixture
def known_cds_WP_003948652_1_len():
    """Known length of cds sequence for WP_003948652.1."""
    known_cds_WP_003948652_1_seq = [len(record.seq) for record in SeqIO.parse('fixtures/SCO_extracted_cds/OG0002119_cds.fasta', 'fasta')][0]

    return known_cds_WP_003948652_1_seq

@pytest.fixture
def known_prot_WP_003948652_1_len():
    """Known length of protein sequence x 3 for WP_003948652.1"""
    knwon_prot_WP_003948652_1_seq = [len(record.seq) for record in SeqIO.parse("fixtures/SCO_genbank_prot_seq/OG0002119_prot.fasta", "fasta")][0]*3

    return knwon_prot_WP_003948652_1_seq


"""test_1_WP_003948652_1_len - Known/original length of cds sequence WP_003948652_1 and length of it's 
corresponding protein sequence x3 was checked. As this sequence was observed to show behaviour 1, 
where the cds sequence length and protein sequence length x3 matched, the test should pass."""

def test_1_WP_003948652_1_len(known_prot_WP_003948652_1_len, known_cds_WP_003948652_1_len):
    """Test if len(known_prot_record)*3 == len(knwon_cds_record). Test should pass."""

    assert known_prot_WP_003948652_1_len == known_cds_WP_003948652_1_len


"""test_2_WP_003948652_1_len - The trimming algorith was applied to WP_003948652_1 cds and protein sequences.
It was expected to see that the length of these sequences should match the following statement: 
len(trimmed_cds) == len(trimmed_protein)*3, and the test should pass."""

def test_2_WP_003948652_1_len():
    """Test if len(prot_record)*3 == len(cds_record) after trimming function is applied. Test should pass."""
    trim_WP_003948652_1 = trim_cds_protein_sequences('test_input/SCO_extracted_cds/OG0002119_cds.fasta')
    trimmed_WP_003948652_1_cds_record_len = [len(record.seq) for record in trim_WP_003948652_1[0]][0]
    trimmed_WP_003948652_1_prot_record_len = [len(record.seq) for record in trim_WP_003948652_1[1]][0]
    assert trimmed_WP_003948652_1_cds_record_len == trimmed_WP_003948652_1_prot_record_len*3

"""test_3_WP_003948652_1_len - Here, the length of original WP_003948652_1 cds sequence is checked against trimmed cds sequence.
As we known that the sequences should not be trimmed by the written function, the test should pass."""

def test_3_WP_003948652_1_len(known_cds_WP_003948652_1_len):
    """Test if the length of known WP_003948652_1 cds is the same as the length of
     WP_003948652_1 after trimming function is applied.
     Test should pass, as the sequence should not be trimmed!
    """
    trim_WP_003948652_1 = trim_cds_protein_sequences('test_input/SCO_extracted_cds/OG0002119_cds.fasta')
    trimmed_WP_003948652_1_cds_record_len = [len(record.seq) for record in trim_WP_003948652_1[0]][0]

    assert trimmed_WP_003948652_1_cds_record_len == known_cds_WP_003948652_1_len


"""test_4_WP_003948652_1_len - Here, the length of original WP_003948652_1 protein sequence is checked against trimmed protein sequence.
The lengths of sequences should be matched, as no trimming should be performed, the test shoul pass. """

def test_4_WP_003948652_1_len(known_prot_WP_003948652_1_len):
    """Test if the length of known WP_003948652_1 protein is the same as the length of
     WP_003948652_1 after trimming function is applied.
     Test should pass, as the sequence should not be trimmed!
    """
    trim_WP_003948652_1 = trim_cds_protein_sequences('test_input/SCO_extracted_cds/OG0002119_cds.fasta')
    trimmed_WP_003948652_1_prot_record_len = [len(record.seq) for record in trim_WP_003948652_1[1]][0]
    assert trimmed_WP_003948652_1_prot_record_len*3 == known_prot_WP_003948652_1_len



@pytest.fixture
def known_cds_WP_200743075_1_len():
    """Length of the known cds sequence for WP_200743075.1"""
    known_cds_WP_200743075_1_seq = [len(record.seq) for record in SeqIO.parse('fixtures/SCO_extracted_cds/OG0002124_cds.fasta', 'fasta')][0]

    return known_cds_WP_200743075_1_seq

@pytest.fixture
def known_prot_WP_200743075_1_len():
    """Length of the known protein WP_200743075.1 sequence x 3."""
    known_prot_WP_200743075_1_seq = [len(record.seq) for record in SeqIO.parse("fixtures/SCO_genbank_prot_seq/OG0002124_prot.fasta", "fasta")][0]*3

    return known_prot_WP_200743075_1_seq



"""test_1_WP_200743075_1_len - Known/original length of cds sequence WP_200743075_1 and length of it's corresponding
protein sequence x3 was checked. As the WP_200743075_1 sequence was observed to show behaviour 1, where the cds sequence (295bp long)
was longer than the protein sequence x 3 (294bp long). The test should fail. """

def test_1_WP_200743075_1_len(known_prot_WP_200743075_1_len, known_cds_WP_200743075_1_len):
    """Test if len(known_prot_record)*3 == len(knwon_cds_record). Test should fail."""

    assert known_prot_WP_200743075_1_len == known_cds_WP_200743075_1_len


"""test_2_WP_200743075_1_len - The trimmed length of cds sequence WP_200743075_1, and the length of it's corresponding 
protein sequence x3 was checked. As the function was written to adequately trim sequences to match for backtranslation, we expect
the test to pass. """
def test_2_WP_200743075_1_len():
    """Test if len(trimmed_cds) == len(trimmed_prot)*3. As we expect, the function
    to solve the problem with the length of sequences for bactranslation, we extect the test to pass."""
    trimmed_record_WP_200743075_1 = trim_cds_protein_sequences('fixtures/SCO_extracted_cds/OG0002124_cds.fasta')
    trimmed_WP_200743075_1_cds_record_len = [len(record.seq) for record in trimmed_record_WP_200743075_1[0]][0]
    trimmed_WP_200743075_1_prot_record_len = [len(record.seq) for record in trimmed_record_WP_200743075_1[1]][0]
    assert trimmed_WP_200743075_1_cds_record_len == trimmed_WP_200743075_1_prot_record_len*3


"""test_3_WP_200743075_1_len - The trimmed cds and original sequences WP_200743075_1 lengths are checked. As in this case the length
of the cds sequence was longer than protein sequences x3, the cds sequence should be trimmed and the lenths should not match.
The test should fail. """

def test_3_WP_200743075_1_len(known_cds_WP_200743075_1_len):
    """Test if the length of known WP_200743075_1_len cds is the same as the length of
     WP_200743075_1_len after trimming function is applied.
     Test should fail, as the cds sequence should be trimmed!
    """
    trimmed_record_WP_200743075_1 = trim_cds_protein_sequences('fixtures/SCO_extracted_cds/OG0002124_cds.fasta')
    trimmed_WP_200743075_1_cds_record_len = [len(record.seq) for record in trimmed_record_WP_200743075_1[0]][0]
    assert trimmed_WP_200743075_1_cds_record_len == known_cds_WP_200743075_1_len


"""test_4_WP_200743075_1_len - The trimmed protein and original sequences for WP_200743075_1 lengths are checked. The test should pass,
as we did not expect or want to trim them. """

def test_4_WP_200743075_1_len(known_prot_WP_200743075_1_len):
    """Test if the length of known WP_200743075_1_len protein is the same as the length of
     WP_200743075_1_len after trimming function is applied.
     Test should pass, as the protein sequence should not be trimmed, as opposed to cds sequences.
    """
    trimmed_record_WP_200743075_1 = trim_cds_protein_sequences('fixtures/SCO_extracted_cds/OG0002124_cds.fasta')
    trimmed_WP_200743075_1_prot_record_len = [len(record.seq) for record in trimmed_record_WP_200743075_1[1]][0]
    assert trimmed_WP_200743075_1_prot_record_len*3 == known_prot_WP_200743075_1_len



@pytest.fixture
def known_cds_WP_086721276_1_len():
    """Length of the known cds sequence for WP_086721276.1"""
    known_cds_WP_086721276_1_seq = [len(record.seq) for record in SeqIO.parse('fixtures/SCO_extracted_cds/OG0002164_cds.fasta', 'fasta')][0]

    return known_cds_WP_086721276_1_seq

@pytest.fixture
def known_prot_WP_086721276_1_len():
    """Length of the known protein WP_086721276.1 sequence x 3."""
    known_prot_WP_086721276_1_seq = [len(record.seq) for record in SeqIO.parse("fixtures/SCO_genbank_prot_seq/OG0002164_prot.fasta", "fasta")][0]*3

    return known_prot_WP_086721276_1_seq



"""test_1_WP_086721276_1_len - The WP_086721276.1 sequence displayed the behaviour 3, where the protein sequence x 3 was longer than the corresponding
cds sequence. Here we check, the length of cds WP_086721276.1 sequence (2287) against the length of protein sequence x3 (2298). The test should
fail.  """

def test_1_WP_086721276_1_len(known_prot_WP_086721276_1_len, known_cds_WP_086721276_1_len):
    """Test if len(known_prot_record)*3 == len(knwon_cds_record). Test should fail."""

    assert known_prot_WP_086721276_1_len == known_cds_WP_086721276_1_len


"""test_2_WP_086721276_1_len - The trimmed protein sequence legth x3 and the trimmed cds sequence were checked. We would expect the written function
to deal with the problem adequately, and solve the problem so that they should be used for backtranslation. Test should pass."""
def test_2_WP_086721276_1_len():
    """Test if len(trimmed_cds) == len(trimmed_prot)*3. As we expect, the function
    to solve the problem with the length of sequences for bactranslation, we extect the test to pass."""

    trimmed_record_WP_086721276_1 = trim_cds_protein_sequences('fixtures/SCO_extracted_cds/OG0002164_cds.fasta')
    trimmed_cds_WP_086721276_1_record_len = [len(record.seq) for record in trimmed_record_WP_086721276_1[0]][0]
    trimmed_prot_WP_086721276_1_record_len = [len(record.seq) for record in trimmed_record_WP_086721276_1[1]][0]
    assert trimmed_cds_WP_086721276_1_record_len == trimmed_prot_WP_086721276_1_record_len*3


"""test_3_WP_086721276_1_len - The length of known/original WP_086721276.1 cds and the length of trimmed cds sequence were checked. We would expect the
written function to trim the sequnece, and the test should fail. """
def test_3_WP_086721276_1_len(known_cds_WP_086721276_1_len):
    """Test if the length of known WP_086721276.1 cds is the same as the length of
     WP_086721276.1 after trimming function is applied.
     Test should fail, as the cds sequence should be trimmed!
    """
    trimmed_record_WP_086721276_1 = trim_cds_protein_sequences('fixtures/SCO_extracted_cds/OG0002164_cds.fasta')
    trimmed_cds_WP_086721276_1_record_len = [len(record.seq) for record in trimmed_record_WP_086721276_1[0]][0]

    assert trimmed_cds_WP_086721276_1_record_len  == known_cds_WP_086721276_1_len


"""test_4_WP_086721276_1_len - The length of known/original WP_086721276.1 x 3 and the length of trimmed protein sequence x 3 were checked. 
As we would expect the function to also trim the protein sequence, thetest should fail."""

def test_4_WP_086721276_1_len(known_prot_WP_086721276_1_len):
    """Test if the length of known WP_086721276.1 protein is the same as the length of
     WP_086721276.1 after trimming function is applied.
     Test should fail, as the protein sequence should be trimmed!
    """
    trimmed_record_WP_086721276_1 = trim_cds_protein_sequences('fixtures/SCO_extracted_cds/OG0002164_cds.fasta')
    trimmed_prot_WP_086721276_1_record_len = [len(record.seq) for record in trimmed_record_WP_086721276_1[1]][0]
    assert trimmed_prot_WP_086721276_1_record_len*3 == known_prot_WP_086721276_1_len