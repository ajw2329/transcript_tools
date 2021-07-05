import transcript_tools.seq_methods
from transcript_tools.transcript_classes import Transcript, Exon
import pytest
import pyranges as pr


@pytest.fixture
def toy_exon_plus():
	'''
	Creates an exon in the '+' strand of a single toy chromosome
	'''

	ex = Exon("chr1", 30, 40, "+")

	return ex

@pytest.fixture
def toy_exon_minus():
	'''
	Creates an exon in the '+' strand of a single toy chromosome
	'''

	ex = Exon("chr1", 30, 40, "-")

	return ex


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

def test_exon_plus_basic_attributes(toy_exon_plus):

	assert toy_exon_plus.chrom == "chr1"
	assert toy_exon_plus.start == 30
	assert toy_exon_plus.end == 40
	assert toy_exon_plus.strand == "+"

def test_exon_plus_length(toy_exon_plus):

	assert toy_exon_plus.length == 11

def test_exon_plus_sequence(toy_exon_plus):

	seq_dict = {"chr1": "ATAGCAGTAGACAGAGAGGCAGACCCACGACAAGGACCCAGATTAGCAGACCCAGACAGGGGCTATTTACCGA"}

	toy_exon_plus.get_seq(seq_dict)

	assert toy_exon_plus.sequence == "ACAAGGACCCA"



def test_exon_plus_boundary_alteration(toy_exon_plus):

	seq_dict = {"chr1": "ATAGCAGTAGACAGAGAGGCAGACCCACGACAAGGACCCAGATTAGCAGACCCAGACAGGGGCTATTTACCGA"}

	toy_exon_plus.get_seq(seq_dict)

	toy_exon_plus.alter_boundaries(32,38)

	assert toy_exon_plus.chrom = "chr1"
	assert toy_exon_plus.start = 32
	assert toy_exon_plus.end = 38
	assert toy_exon_plus.strand = "+"
	assert toy_exon_plus.sequence = "AGGACC"


def test_exon_minus_basic_attributes(toy_exon_minus):

	assert toy_exon_minus.chrom == "chr1"
	assert toy_exon_minus.start == 30
	assert toy_exon_minus.end == 40
	assert toy_exon_minus.strand == "-"

def test_exon_minus_length(toy_exon_minus):

	assert toy_exon_minus.length == 11

def test_exon_minus_sequence(toy_exon_minus):

	seq_dict = {"chr1": "ATAGCAGTAGACAGAGAGGCAGACCCACGACAAGGACCCAGATTAGCAGACCCAGACAGGGGCTATTTACCGA"}

	toy_exon_minus.get_seq(seq_dict)

	assert toy_exon_minus.sequence == "TGGGTCCTTGT"



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