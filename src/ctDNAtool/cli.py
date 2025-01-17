import click

from . import generators
from . import preprocessors
from . import manipulations
from .utils import tsv_reader
from .preprocessors.bin_genome import Chromosomes
from . import cli_common


@click.group()
@cli_common.quiet
@cli_common.debug
def cli(quiet, debug):
    cli_common.setup_debugger(quiet, debug)


@cli.command()
@click.argument("annotation_input_file")
@click.option(
    "-k",
    "--region-size",
    default=10000,
    type=click.IntRange(min=0),
    help="Size of region around tss",
)
@click.option("--bed-output-file", default="transcription_start_sites.bed")
@click.option("--tss-output-file", default="transcription_start_sites.tsv")
def find_tss(annotation_input_file, region_size, bed_output_file, tss_output_file):
    """Finds all Transcription Start Sites given an annotation file."""
    preprocessors.find_tss(
        annotation_input_file, region_size, bed_output_file, tss_output_file
    )


@cli.command()
@click.argument("genome_ref_file")
@click.option("-o", "--output-file")
@cli_common.mbp
@cli_common.include_x
def bin_genome(genome_ref_file, output_file, mbp, include_x):
    """Creates a bed file, splitting the genome in bins of the given mbp size"""
    chromosomes = Chromosomes.AUTOSOMES_X if include_x else Chromosomes.AUTOSOMES
    preprocessors.bin_genome_Mbp(genome_ref_file, output_file, mbp, chromosomes)


@cli.command()
@click.argument("genome-ref-file")
@click.option("-o", "--output_file")
@cli_common.include_x
def bin_genome_chromosome(genome_ref_file, output_file, include_x):
    """Creates a bed file, splitting the genome into bins, where each bin corresponds to a chromosome"""
    chromosomes = Chromosomes.AUTOSOMES_X if include_x else Chromosomes.AUTOSOMES
    preprocessors.bin_genome_chromosome(genome_ref_file, output_file, chromosomes)


@cli.command()
@click.argument("bam_file")
@click.argument("bed_file")
@click.option("-o", "--output-file", default="length_matrix.pickle")
@cli_common.max_length
@cli_common.map_quality
def generate_length(bam_file, bed_file, output_file, max_length, map_quality):
    """Creates a tensor with fragment length data"""
    generators.length_matrix(bam_file, bed_file, output_file, max_length, map_quality)


@cli.command()
@click.argument("bam_file")
@click.argument("bed_file")
@click.argument("reference_genome")
@click.option("-o", "--output-file", default="length_seq_matrix.pickle")
@cli_common.max_length
@cli_common.flank
@cli_common.map_quality
@cli_common.end_type
def generate_length_end_seq(
    bam_file,
    bed_file,
    reference_genome,
    output_file,
    max_length,
    flank,
    map_quality,
    end_type,
):
    """Creates a tensor with length and end sequence data"""
    generators.length_end_seqs(
        bam_file,
        bed_file,
        reference_genome,
        output_file,
        max_length,
        flank,
        map_quality,
        end_type,
    )


@cli.command()
@click.argument("bam_file")
@click.argument("bed_file")
@click.argument("reference_genome")
@click.option("-o", "--output-file", default="length_seq_matrix.pickle")
@cli_common.max_length
@cli_common.flank
@cli_common.map_quality
def generate_length_end_seq_marginal(
    bam_file, bed_file, reference_genome, output_file, max_length, flank, map_quality
):
    """Creates a tensor with length and marginal end sequence data"""
    generators.length_end_seqs_marginal(
        bam_file,
        bed_file,
        reference_genome,
        output_file,
        max_length,
        flank,
        map_quality,
    )


@cli.command()
@click.argument("bam_file")
@click.argument("bed_file")
@click.argument("reference_genome")
@click.option("-o", "--output-file", default="length_matrix.pickle")
@cli_common.max_length
@cli_common.flank
@cli_common.map_quality
def generate_mate_length_end_seq(
    bam_file, bed_file, reference_genome, output_file, max_length, flank, map_quality
):
    """Create a tensor with length and end sequence data, where the first dimension represents whether a read came from the first or the second mate"""
    generators.mate_length_end_seqs(
        bam_file,
        bed_file,
        reference_genome,
        output_file,
        max_length,
        flank,
        map_quality,
    )


@cli.command()
@click.argument("sample_files", nargs=-1)
@cli_common.file_of_files
@click.option("-o", "--output-file", default="collapsed_samples.pickle")
@click.option("--uint32", is_flag=True)
def sample_sum(sample_files, file_of_files, output_file, uint32):
    """Collapses the samples value by value"""
    if file_of_files:
        files = cli_common.get_files_from_file(file_of_files)
        manipulations.sample_sum(files, output_file, uint32)
    else:
        if len(sample_files) > 0:
            manipulations.sample_sum(sample_files, output_file, uint32)


@cli.command()
@click.argument("sample_file")
@click.option("-o", "--output-file", default="collapsed_sample.pickle")
def region_sum(sample_file, output_file):
    """Sums the regions of the sample file"""
    manipulations.region_sum(sample_file, output_file)


@cli.command()
@click.argument("input_file")
@click.option("-o", "--output-file", default="tsv_length_matrix.csv")
@cli_common.min_length
@cli_common.max_length
def convert_to_tsv_length(input_file, output_file, min_length, max_length):
    """Converts a .pickle file containing length data to a .tsv file"""
    manipulations.convert_to_tsv_length(input_file, output_file, min_length, max_length)


@cli.command()
@click.argument("input_file")
@click.option("-o", "--output-file", default="tsv_length_seq.csv")
def convert_to_tsv_length_seq_sum(input_file, output_file):
    """Converts a .pickle file containing length seq data to a .tsv file"""
    manipulations.convert_to_tsv_length_seq_sum(input_file, output_file)


@cli.command()
@click.argument("input_file")
@click.option("-o", "--output-file", default="tsv_length_seq.csv")
def convert_to_tsv_length_seq_sum_wide(input_file, output_file):
    """Converts a .pickle file containing length seq data to a .tsv file"""
    manipulations.convert_to_tsv_length_seq_sum_wide(input_file, output_file)


@cli.command()
@click.option("-o", "--output-file", default="combined_data.pickle")
@cli_common.file_of_files
@click.argument("input_files", nargs=-1)
def combine_data(output_file, file_of_files, input_files):
    """Combines multiple .pickle files with Data objects into one .pickle file. Can optionally take a file of files as input"""

    if file_of_files:
        files = cli_common.get_files_from_file(file_of_files)
        manipulations.combine_data(output_file, files)
    else:
        manipulations.combine_data(output_file, input_files)


@cli.command()
@click.argument("input_sample")
@click.argument("ids_file")
@click.option("-o", "--output-file", default="subset_sample.pickle")
def pick_subset(input_sample, ids_file, output_file):
    """Creates a subset of the sample given a list of row identifiers"""
    ids = list()
    with open(ids_file) as fp:
        for line in tsv_reader(fp):
            if line[0].startswith("#"):
                continue
            ids.append(line[0])
    manipulations.pick_subset(input_sample, output_file, ids)


@cli.command()
@click.argument("input_matrix")
@click.option("-o", "--output-file", default="binned_matrix.pickle")
@click.option(
    "-b", "--bin-size", default=1, type=click.IntRange(min=1)
)  # TODO: Could this be @mbp as well?
@cli_common.stride
def binning(input_matrix, output_file, bin_size, stride):
    """Bins the sample in the first axis."""
    manipulations.binning(input_matrix, output_file, bin_size, stride)
