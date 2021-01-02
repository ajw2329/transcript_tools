import cds_tools.seq_methods
from cds_tools.transcript_classes import Transcript, Exon
import pytest
import pyranges as pr

@pytest.fixture
def toy_transcript_plus():
	'''
	Creates a transcript object on the '+' strand of a single chromosome
	'''

	seq_dict = {"chr1": "ATAGCAGTAGACAGAGAGGCAGACCCACGACAAGGACCCAGATTAGCAGACCCAGACAGGGGCTATTTACCGA"}

	tx = Transcript(
			seq_dict = seq_dict, 
			exon_tuple = (("chr1", 1, 10, "+", None), ("chr1", 30, 40, "+", None)))

	tx.complete_transcript()

	return tx

@pytest.fixture
def toy_transcript_minus():
	'''
	Creates a transcript object on the '+' strand of a single chromosome
	'''

	seq_dict = {"chr1": "ATAGCAGTAGACAGAGAGGCAGACCCACGACAAGGACCCAGATTAGCAGACCCAGACAGGGGCTATTTACCGA"}

	tx = Transcript(
			seq_dict = seq_dict, 
			exon_tuple = (("chr1", 1, 10, "-", None), ("chr1", 30, 40, "-", None)))

	tx.complete_transcript()

	return tx


def test_transcipt_complete_sorting_plus(toy_transcript_plus):

	assert toy_transcript_plus.exons[0].start == 1 and toy_transcript_plus.exons[1].start == 30


def test_transcipt_complete_sorting_minus(toy_transcript_minus):

	assert toy_transcript_minus.exons[0].start == 30 and toy_transcript_minus.exons[1].start == 1


def test_transcript_sequence_plus(toy_transcript_plus):

	assert toy_transcript_plus.sequence == 'ATAGCAGTAGACAAGGACCCA'


def test_transcript_sequence_minus(toy_transcript_minus):

	assert toy_transcript_minus.sequence == 'TGGGTCCTTGTCTACTGCTAT'

def test_transcript_ranges_plus(toy_transcript_plus):
	'''
	Tests if pyranges created by Transcript.complete_transcript
	is equivalent to manually specified pyranges.  To do this,
	the PyRanges.as_df method must be called to convert to a
	pandas DataFrame, after which the pandas.DataFrame.equals
	method can be called.  It does not look like the 
	__eq__ method was specified in PyRanges
	'''

	assert toy_transcript_plus.ranges.as_df().equals(pr.PyRanges(
			chromosomes = ["chr1", "chr1"],
			starts = [1, 30],
			ends = [10, 40],
			strands = ["+", "+"]).as_df())