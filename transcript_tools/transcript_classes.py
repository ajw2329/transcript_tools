
## TODO 
### implement conversion of genomic to txomic coords
### Add check for stop codon from downstream start codon
### Plan out how to implement Transcriptome class that allows
###		rapid intersection of start codons with exons using
### 	pyranges

### decide if class should take seq dict as input
### when an exon is sliced, its sequence should also get sliced - add this code and test it

import pyranges as pr
import transcript_tools.seq_methods as seq_methods


def multi_index(query_list, search_term_list, retain_none = False):
	'''
	Wraps index() such that a list of objects to identify the
	indices of can be passed instead of a single object. Further,
	if the object is not contained in the query list, the index
	will either not be returned, or will optionally be returned
	as None.


	Parameters
	----------
	
	query_list : list
		list to be searched

	search_term_list : list
		list of objects to identify the indices of
	
	retain_none : bool
		Boolean indicating whether indices of search terms not found
		in the query list will be returned as None or will simply
		be omitted.

	Returns
	-------

	index_list : list
		list of indices for the search terms. If retain_none is set
		to True, the returned indices will be in the same order as 
		the search term list, with a NoneType value indicating 
		that the search term was not found in the query list. If
		retain_none is Set to False (the default), the order of
		the returned list is not necessarily meaningful. This is
		useful if, for instance, one is searching for the earliest
		index for all search terms.

	>>> multi_index(['ATT', 'AGG', 'GGG', 'TAG', 'AAA', 'CCC', 'TAG'], ["TAA", "TAG"])
	[3]

	>>> multi_index(['ATT', 'AGG', 'GGG', 'TAG', 'AAA', 'CCC', 'TAG'], ["TAA", "TAG"], retain_none = True)
	[None, 3]	
	'''

	indices = [None]*len(search_term_list)

	for index, term in enumerate(search_term_list):

		try:
			search_index = query_list.index(term)
			indices[index] = search_index
		except ValueError:
			indices[index] = None

	if not retain_none:

		indices = [i for i in indices if i is not None]

	return indices


class TranscriptError(Exception):
    """Base class for exceptions in this module."""
    pass

class TranscriptCompletionError(TranscriptError):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class Transcript():
	'''
	Transcript - a collection of transcribed exons spliced together,
	including chromosomal coordinates and strand, as well as
	spliced sequence.  Transcript objects may also contain
	CodingTranscript objects, which contain information about
	a putative coding sequence and associated sequelae, such as
	the presence of a putative PTC.  Note that a given Transcript 
	object may contain multiple associated CodingTranscript objects,
	due to the fact that a given set of spliced exons may in 
	principle by translated from more than one start codon.

	Parameters
	-----------
	exon_tuple : tuple
		Tuple of tuples containing information about component exons. 
		Can be added manually after instantiation with Transcript.add_exon() 
	require_common_strand : bool
		Logical indicating whether transcript is required to be composed of 
		exons from the same strand 
	require_common_chrom : bool
		Logical indicating whether transcript is required to be composed of 
		exons from the same chromosome
	require_non_overlap : bool
		Logical indicating whether chromosomal coordinates of exons are 
		required to be non-overlapping
	allow_nonstop : bool
		Logical indicating whether candidate coding sequences added to 
		the transcript object should include any putative nonstop
		scenarios (i.e. situations in which translation from a specified
		start codon can not result in an in-frame stop codon. If such a 
		scenario were to occur, the transcript might be targeted by
		surveillance processes such as nonstop decay)

	Attributes
	----------
	exon_complete : bool
		Logical indicating whether transcript contains all constituent exons. If false,
		more exons need to be added, followed by subsequent transcript object processing
		steps.
	exons : tuple
		Tuple of Exon objects characterizing the constituent exons
	junctions : tuple
		Tuple of strings characterizing splice junctions between exons
	chrom : str
		Chromosome on which the transcript is encoded (relevant only when exons
		are confined to a single chromosome)
	strand : str
		Strand ('+' or '-') on which the transcript is encoded (relevant only 
		when exons are confined to a single strand)
	cumulative_exon_lengths : tuple
		Values containing the cumulative lengths of all exons in the transcript.

	seq_dict = seq_methods.import_fasta("../data/GRCh38.primary_assembly.genome.fa.gz") 
	x = Transcript(seq_dict, (("chr1", 1000000, 1000100, "+"), ("chr1", 1000300, 1000400, "+")))
	'''

	def __init__(
		self, 
		seq_dict = None,
		exon_tuple = None, 
		allow_nonstop = False,
		allowable_start_codons = ["AUG"],
		allowable_stop_codons = ["TAG", "TGA", "TAA"]):

		if exon_tuple is not None:

			if seq_dict is None:

				raise ValueError("A seq_dict must be supplied if exon information is supplied.")

		self.allow_nonstop = allow_nonstop
		self.allowable_start_codons = allowable_start_codons
		self.allowable_stop_codons = allowable_stop_codons
		self.exon_complete = False
		self.exon_sorted = False
		self.exons = ()

		for exon in exon_tuple:

			self.add_exon(exon[0], exon[1], exon[2], exon[3], exon[4], seq_dict)



	def add_exon(self, chrom, start, end, strand, position_in_tx, seq_dict):
		'''
		Add another exon to a Transcript object

		Parameters
		----------
		chrom : str
			Chromosome on which the exon is located
		start : int
			1-based starting coordinate of the exon
		end : int
			1-based end coordinate of the exon
		strand : str
			Strand of the chromosome on which the exon is encoded. 
			Can be either "+" or "-".
		position_in_tx 
			Relative 1-based position of the exon within the transcript
			(e.g. 2 if the exon is the second exon from the 5'-end). 
			None can be passed if the transcript is assumed to consist
			of consecutive exons on the same strand of a single chromosome.
		seq_dict : dict
			Dictionary containing chromosome sequences indexed by 
			chromosome names. 

		Raises
		------
		TranscriptCompletionError
			If a Transcript is marked as complete new exons cannot be 
			added. If this method is called on a complete Transcript
			object (i.e. Transcript.exon_complete == True), this
			error will be raised.

		'''

		if not self.exon_complete:

			new_exon = Exon(chrom, start, end, strand, position_in_tx)

			new_exon.get_seq(seq_dict)

			self.exons += tuple([new_exon])

		else:

			raise TranscriptCompletionError(
				"Can't add exon to exon-complete transcript.")



	def check_exon_positions(self):
		'''
		Check whether transcript exons were provided with their
		position in the transcript. This is not necessary  fortranscripts
		consisting of sequential, non-overlapping exons from the
		same strand of the same chromosome, but is necessary for
		unconventional transcripts such as trans-spliced transcripts.

		Raises
		------

		TranscriptCompletionError
			If some but not all of the exons were provided with their
			position within the transcript, the final exon order
			cannot be determined, so an exception will be raised.
		'''

		exon_positions = set([exon.position_in_tx for exon in self.exons if exon.position_in_tx is not None])

		if len(exon_positions) == 0:

			self.use_exon_positions = False

		elif len(exon_positions) == len(self.exons):

			self.use_exon_positions = True

		else:

			raise TranscriptCompletionError("Some but not all exons were provided transcript positions.")



	def sort_exons(self):
		'''
		Sorts exons according to their position in the transcript,
		from 5'-> 3'. Attempts to use Exon.position_in_tx, if 
		available. If unavailable, the exons will be sorted 
		according to chromosomal coordinates (from lowest to
		highest in the case of exons on the '+' strand, and 
		in reverse in the case of exons on the '-' strand).

		'''

		if self.use_exon_positions:

			self.exons = sorted(self.exons, key = lambda x: x.position_in_tx)

		else:

			self.exons = sorted(self.exons, key = lambda x: x.start)

		self.exon_sorted = True


	def check_common_chrom(self):
		'''
		Checks whether all provided exons are encoded on 
		the same chromosome.

		Raises
		------
		TranscriptCompletionError
			If provided exons are encoded on more than one chromosome, 
			an exception is raised. This method is only called when
			multiple chromosomes are disallowed.
		'''

		chromosomes = set([exon.chrom for exon in self.exons])

		if len(chromosomes) > 1:

			raise TranscriptCompletionError("Exons from multiple chromosomes provided when disallowed.")

		else:

			self.chrom = chromosomes.pop()

	def check_common_strand(self):
		'''
		Checks whether all provided exons are encoded on 
		the same strand.

		Raises
		------
		TranscriptCompletionError
			If provided exons are encoded on more than one strand 
			an exception is raised. This method is only called when
			multiple strands are disallowed.
		'''

		strands = set([exon.strand for exon in self.exons])

		if len(strands) > 1:

			raise TranscriptCompletionError("Exons from multiple strands provided when disallowed.")

		else:

			self.strand = strands.pop()

	def check_exon_overlap(self):
		'''
		Checks whether exons on the same chromosome overlap
		'''

		for index, exon in enumerate(self.exons):

			if index < len(self.exons) - 1:

				if exon.end >= self.exons[index + 1].start:

					raise TranscriptError("Coordinates of exons overlap when disallowed.")


	def get_splice_junctions(self):
		'''
		Determines splice junction identities from the complete and ordered
		transcript exon list.
		'''

		if not self.exon_complete:

			raise TranscriptCompletionError(
				"Can't get splice junctions if the transcript isn't exon-complete.")

		self.junctions = [None]*(len(self.exons) - 1)

		for index, exon in enumerate(self.exons):

			if index < len(self.exons) - 1:

				nindex = index + 1

				junction = "_".join((
					exon.chrom, 
					exon.strand, 
					str(exon.end), 
					self.exons[nindex].chrom,
					self.exons[nindex].strand,
					str(self.exons[nindex].start)))

				self.junctions[index] = junction

		self.junctions = tuple(self.junctions)


	def get_transcript_length(self):
		'''
		Determine the total length of the 
		transcript by summing the lengths of the exons
		'''

		if not self.exon_complete:

			raise TranscriptCompletionError(
				"Can't get meaningful transcript length before the transcript is complete")

		self.transcript_length = sum(exon.length for exon in self.exons)


	def get_transcript_boundaries(self):
		'''
		Assign transcript genomic start and end as separate atrributes. 
		"start" and "end" refer to the "leftmost" and "rightmost" 
		genomic coordinates respectively, not necessarily the 5' and 3'
		ends of the RNA transcript itself.
		'''

		if not self.exon_complete:

			raise TranscriptCompletionError(
				"Can't get transcript boundaries before the transcript is complete")

		self.start = self.exons[0].start
		self.end = self.exons[-1].end


	def get_transcript_sequence(self):
		'''
		Concatenate exon sequences into transcript sequence
		'''

		seqs = [exon.sequence for exon in self.exons]

		self.sequence = "".join(seqs)

		assert len(self.sequence) == self.transcript_length


	def create_pyranges(self):
		'''
		Create pyranges object with transcript exons.
		'''

		self.ranges = pr.PyRanges(
			chromosomes = [exon.chrom for exon in self.exons],
			starts = [exon.start for exon in self.exons],
			ends = [exon.end for exon in self.exons],
			strands = [exon.strand for exon in self.exons])


	def assert_exon_complete(self):
		'''
		Sets exon_complete attribute to `True`, asserting
		that all constituent exons have been added to the
		transcript.
		'''

		if len(self.exons) == 0:

			raise TranscriptCompletionError('Transcript with zero exons cannot be complete.')

		self.exon_complete = True


	def complete_transcript(self):
		'''
		Runs methods that flesh out transcript features after
		all exons have been added to the object.
		'''

		self.assert_exon_complete()

		self.check_common_chrom()

		self.check_common_strand()

		self.check_exon_overlap()

		## check presence of exon positions

		self.check_exon_positions()

		## sort exons

		self.sort_exons()

		## calculate transcript length

		self.get_transcript_length()

		## get transcript gx start and end

		self.get_transcript_boundaries()

		## get transcript sequences

		self.get_transcript_sequence()

		## add junctions

		self.get_splice_junctions()

		## add pyranges

		self.create_pyranges()

		## calculate cumulative exon lengths

		self.calc_cumulative_exon_lengths()


	def assert_exon_incomplete(self):
		'''
		Resets exon_complete attribute to `False`.
		'''

		self.exon_complete = False


	def add_cds(self, chrom, start, strand):
		'''
		Attempts to add a CDS to a transcript by translating in silico
		from the provided start codon position
		'''

		if not self.exon_complete:

			raise TranscriptCompletionError(
				"CDS cannot be added to a transcript until Transcript.complete_transcript has been called")

		else:

			if self.position_in_exons(chrom, start, strand):

				cds_tx_start = self.convert_gx_to_tx_coords(chrom, start, strand, "GT")

				unbounded_candidate_seq = self.sequence[cds_tx_start - 1:]

				candidate_codons = seq_methods.get_non_overlapping_codon_triplets(unbounded_candidate_seq)


	def calc_cumulative_exon_lengths(self):
		'''
		Calculates the cumulative lengths of the ordered exons in the transcript. A
		the cumulative lengths are stored as a tuple, with a 0 prepended to the beginning
		of the tuple for convenience in coordinate calculations. The value at a given 
		index can thereby be interpreted as 'the length up until the start of the exon at
		the same index'. For instance, for a transcript with 5 exons of length 100,
		self.cumulative_exon_lengths would adopt the value (0, 100, 200, 300, 400, 500).
		'''

		if not self.exon_complete:

			raise TranscriptCompletionError(
				"CDS cannot be added to a transcript until Transcript.complete_transcript has been called")

		cumulative_exon_lengths = [0]

		for i, exon in enumerate(self.exons):

			cumulative_exon_lengths.append(exon.length + cumulative_exon_lengths[i])

		self.cumulative_exon_lengths = tuple(cumulative_exon_lengths)



	def position_in_exons(self, chrom, start, strand):
		'''
		Takes a genomic position and determines which transcript exon contains it,
		if a transcript exon does indeed contain it.

		Parameters
		----------

		chrom : str
			Chromosome of the putatively overlapping position
		start : int
			Start position of the putatively overlapping position
		strand : str
			Strand ('+' or '-') of the putatively overlapping position

		Returns
		-------

		int
			Index of the containing exon.

		'''


		for i, exon in enumerate(self.exons):

			if chrom == exon.chrom and strand == exon.strand and exon.start <= start <= exon.end:

				return i

		else:

			return None


	def slice_exons(self, start = None, end = None):
		'''
		Slices exons bounded by a set of genomic coordinates. Set bounding coordinates must
		occur within transcript exon boundaries, however, if either start or end is set to
		none, the bounding occurs only at the provided coordinate.  For instance, if start
		is provided by end is not, the returned exons will include genomic space from the
		provided start coordinate until the end of the input exons. At least one bounding
		coordinate must be provided. Chromosome and strand of the bounding coordinates are
		assumed to match. Exons and exon segments falling outside the boundaries are removed 
		from the returned exon list.

		Parameters
		----------

		start : int
			Leftmost of the bounding coordinates
		end : int
			Rightmost of the bounding coordinates
		unbounded : str
			Indicates if the 

		Returns
		-------

		list
			List of contained exons, with outer exons 'shaved' if the bounding coordinates
			do not include their entirety.
		'''

		if start is None and end is None:

			raise ValueError(
				"At least one of 'start' and 'end' must be provided an integer value. Both cannot be 'None'.")

		if start is None:

			start_contained = 0,
			start = 0

		if end is None:

			end_contained = len(self.exons)
			end = exon_list[-1].end

		if start_contained is not None and end_contained is not None:

			exon_list = self.exons[start_contained, end_contained+1]

			exon_list[0].alter_boundaries(new_start = start)
			exon_list[-1].alter_boundaries(new_end = end)

		return exon_list
		


	def convert_gx_to_tx_coords(self, chrom, position, strand, direction = "GT"):
		'''
		Takes a genomic coordinate and converts it to a transcriptomic coordinate,
		or vice versa, provided it overlaps an exon in the transcript.
		'''

		if direction not in ["GT", "TG"]:

			raise ValueError("direction must be GT (genome-to-transcript) or TG (transcript-to-genome)")

		if not self.exon_complete:

			raise TranscriptCompletionError(
				"CDS cannot be added to a transcript until Transcript.complete_transcript has been called")

		if direction == "GT":

			exon_index = self.position_in_exons(chrom, start, strand)

			if exon_index is None:

				raise ValueError("Coordinate cannot be converted because it does not overlap any exons")

			## first get position from the most upstream position with reference to the genome
			### (i.e. the 'leftmost' position)

			overhang = position - self.exons[exon_index].start
			position_from_left = overhang + self.cumulative_exon_lengths[exon_index]

			## then adjust if the transcript is in the "-" strand

			if strand == "-":

				position = self.transcript_length - position_from_left

			else:

				position = position_from_left

			return position

		else:

			pass


	def __repr__(self):

		rep_l = ["Chromosome: " + self.chrom,
				 "Start: " + str(self.start),
				 "End: " + str(self.end),
				 "Strand: " + self.strand]

		return("\n".join(rep_l))


class CodingTranscript(Transcript):
	'''
	Coding transcript, with one CDS per object
	'''

	def __init__(
		self, 
		transcript, 
		cds_g_start,
		cds_g_end,
		cds_tx_start, 
		cds_tx_end,
		nmd_threshold = 55):

		self.cds_tx_start = cds_tx_start
		self.cds_tx_end = cds_tx_end

		pass


	def get_five_utr_exons(self):
		'''
		Recovers the partial and complete exons that
		constitute the 5'-UTR. 
		'''

		self.five_utr_exons = self.slice_exons(end = cds_g_start)


	def get_three_utr_exons(self):
		'''
		Recovers the partial and complete exons that
		constitute the 3'-UTR. 
		'''

		self.three_utr_exons = self.slice_exons(start = cds_g_end)


	def get_coding_exons(self):
		'''
		Recovers the partial and complete exons that
		constitute the coding region.
		'''

		self.coding_exons = self.slice_exons(start = cds_g_start, end = cds_g_end)


	def get_five_utr_seq(self):
		'''
		Recovers the sequence of the 5'-UTR.
		'''

		self.five_utr_seq = "".join([i.sequence for i in self.five_utr_exons])

	def get_three_utr_seq(self):
		'''
		Recovers the sequence of the 3'-UTR
		'''

		self.three_utr_seq = "".join([i.sequence for i in self.three_utr_exons])

	def get_coding_seq(self):
		'''
		Recovers the sequence of the coding region.
		'''

		self.coding_seq = "".join([i.sequence for i in self.coding_exons])

	def translate(self):
		'''
		Recovers the amino acid sequence of the putative protein.
		'''

		pass

	def has_ptc(self):
		'''
		Determines whether the stop codon is substantially
		upstream of the 3'-most splice site.
		'''

		pass

	def find_uorf(self):
		'''
		Identifies potential open reading frames within the 5'-UTR.
		'''

		pass

	def find_dorf(self):
		'''
		Identifies potential open reading frames within the 3'-UTR
		'''

		pass



class Exon():
	'''Perceptron classifier

	Parameters
	-----------
	chrom : str
		Chromosome on which the exon is found
	start : int
		1-based left-most chromosomal coordinate of the exon
	end : int
		1-based right-most chromosomal coordinate of the exon
	strand : str
		Chromosome strand encoding the exon - can be "+" or "-" only
	position_in_tx : int
		1-based position of the exon relative to other exons in the transcript.
		Is inferred from chromosomal coordinates and strand if not supplied.
		Explicit provision of this value can account for unconventional 
		transcripts. Default is None.


	Attributes
	----------
	chrom
	start
	end
	strand
	position_in_txt
	'''

	def __init__(self, chrom, start, end, strand, position_in_tx = None):

		if strand not in ("+", "-"):

			raise ValueError("Exons only allowed to have strand '+' or '-', not " + strand)

		self.chrom = chrom
		self.start = start
		self.end = end
		self.length = end - start + 1
		self.strand = strand
		self.position_in_tx = position_in_tx
		self.sequence = None

	def get_seq(self, seq_dict):
		'''
		Recovers sequence of exon.
		'''

		seq = seq_methods.slice_seq(
			self.chrom, self.start, self.end, seq_dict)

		if self.strand == "+":

			self.sequence = seq

		elif self.strand == "-":

			self.sequence = seq_methods.gen_rev_comp(seq)


	def alter_boundaries(self, new_start = None, new_end = None):
		'''
		Alters the start or end coordinate of the exon, and alters the sequence,
		if available, appropriately. One or both boundaries can be replaced,
		but they must be within the oroginal boundaries.

		Parameters
		----------

		new_start : int
			New start coordinate -- must be >= original start and <= original end
		new_end : int
			New end coordinate -- must be >= original start and  <= original end
		'''

		alter = False

		if new_start is not None:

			if self.start <= new_start <= self.end:

				alter = True

				slice_start_coord = new_start - self.start

			else:

				raise ValueError("New boundaries must be within old boundaries.")

		else:

			new_start = self.start
			slice_start_coord = 0


		if new_end is not None:

			if self.start <= new_end <= self.end:

				alter = True

				slice_end_coord = new_end - self.start

			else:

				raise ValueError("New boundaries must be within old boundaries.")

		else:

			new_end = self.end
			slice_end_coord = None


		if alter:

			self.start = new_start
			self.end = new_end


			if self.sequence is not None:

				self.sequence = self.sequence[slice_start_coord:slice_end_coord + 1]



	def __repr__(self):

		rep_l = ["Chromosome: " + self.chrom,
				 "Start: " + str(self.start),
				 "End: " + str(self.end),
				 "Strand: " + self.strand]

		return("\n".join(rep_l))


class Transcriptome:
	'''
	A dictionary of Transcript objects, along with a single
	pyranges object containing all exons in the transcriptome.
	A transcriptome can be instantiated by passing a path to
	a GTF file, or alternatively by passing a list of 
	transcript objects.
	'''

	def __init__(self):

		pass

	def add_transcript(self, gtf_path = None, transcript_list = None):

		pass

	def create_pyranges(self):

		pass

	def summarize(self):

		pass

	def __repr__(self):

		pass

