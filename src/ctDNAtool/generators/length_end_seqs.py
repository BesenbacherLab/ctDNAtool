import numpy as np
from scipy.sparse import dok_matrix, csr_matrix
import logging

from .bed import load_bed_file
from .bam import BAM
from .utils import seq_to_index, fetch_seqs, reverse_complement
from ..data import Data
from ..py2bit_context import Py2bitContext

logger = logging.getLogger()


def length_end_seqs(
    bam_file, bed_file, ref_genome_file, output_file, max_length=500, flank=1, mapq=20
):
    """Create a tensor where the first dim. represents a region from the bed file,
    the second dim. represent read lengths from 1 to max_length and the third dim.
    represents the sequences at the fragment ends, taken from the reference genome, 
    endcoded as an index. Then length of an end sequence is 2 times the flank parameter.

    :param bam_file: File path to the bam sample file
    :type bam_file: str
    :param bed_file: File path to the bed file, which can be compiled by the preprocessing function
    :type bed_file: str
    :param output_file: File path to the output file
    :type output_file: str
    :param ref_genome_file: File path to the reference genome given as a 2bit file.
    :type ref_genome_file: str
    :param max_length: Maximum read length to be counted
    :type max_length: int > 0
    :param flank: Determines how many base pairs are examined in each end.
    type flank: Int
    :param mapq: map quality. Ignores all reads below the threshold.
    :type mapq: Int
    :returns:  None
    """
    region_lst = load_bed_file(bed_file)
    bam = BAM(bam_file)
    id_lst = list()
    tensor = np.empty((len(region_lst),), dtype=np.object)
    N_seqs = 4 ** (2 * flank)

    with Py2bitContext(ref_genome_file) as tb:
        chroms_lengths = tb.chroms()
        for i, region in enumerate(region_lst):
            matrix = dok_matrix((max_length, N_seqs), dtype=np.uint32)
            for read in bam.pair_generator(
                region.chrom, region.start, region.end, mapq
            ):
                length = read.length
                if (
                    length <= max_length
                    and chroms_lengths[region.chrom] >= (read.end + flank)
                    and 0 <= (read.start - flank)
                ):
                    start_seq, end_seq = fetch_seqs(tb, region.chrom, read.start, read.end, flank)
                    if not 'N' in start_seq:
                        matrix[length - 1, seq_to_index(start_seq)] += 1
                    if not 'N' in end_seq:
                        end_seq = reverse_complement(end_seq)
                        matrix[length - 1, seq_to_index(end_seq)] += 1
            tensor[i] = csr_matrix(matrix)
            id_lst.append(region.region_id)
        logger.info(str(bam))
        Data.write(Data(tensor, id_lst, bam.report), output_file)
