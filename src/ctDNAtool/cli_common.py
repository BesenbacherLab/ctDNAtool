import click
import logging


def include_x(function):
    function = click.option(
        "-x",
        "--include-X",
        is_flag=True,
        help="Option for specifying whether the X chromosome should be included",
    )(function)

    return function


def min_length(function):
    function = click.option(
        "--min-length",
        default=1,
        type=click.IntRange(min=1),
        help="Minimum read length to be counted",
    )(function)

    return function


def max_length(function):
    function = click.option(
        "--max-length",
        default=500,
        type=click.IntRange(min=1),
        help="Maximum read length to be counted",
    )(function)

    return function


def map_quality(function):
    function = click.option(
        "-q",
        "--map-quality",
        default=20,
        type=click.IntRange(min=0),
        help="Minimum map quality",
    )(function)

    return function


def pickle_output(function):
    function = click.option(
        "-p",
        "--pickle-output",
        is_flag=True,
        help="Specifies if output should be to .pickle",
    )(function)

    return function


def mbp(function):
    function = click.option(
        "-m",
        "--mbp",
        default=1.0,
        type=click.FloatRange(min=1e-06),
        help="Binning size given in mega Base Pairs",
    )(function)

    return function


def flank(function):
    function = click.option(
        "-f",
        "--flank",
        default=1,
        type=click.IntRange(min=1),
        help="Number of base pairs to include measured from fragment end",
    )(function)

    return function

def end_type(function):
    function = click.option(
        "--end_type",
        default='both',
        type=click.Choice(['start', 'end', 'both']),
        help='Type of ends to look at. When ends are included the reverse complement of their k-mer will be counted.',
    )(function)

    return function


def stride(function):
    function = click.option(
        "-s",
        "--stride",
        default=1,
        type=click.IntRange(min=1),
        help="Number of base pairs to slide, when binning with a sliding window",
    )(function)

    return function


def file_of_files(function):
    function = click.option(
        "-f", "--file-of-files", help="File containing files to be combined"
    )(function)

    return function


def quiet(function):
    function = click.option(
        "-q", "--quiet", is_flag=True, help="Suppress logging output"
    )(function)

    return function


def debug(function):
    function = click.option(
        "-d", "--debug", is_flag=True, help="Prints debugging info"
    )(function)

    return function


def setup_debugger(quiet_flag, debug_flag):
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()
    if quiet_flag:
        logger.setLevel(logging.CRITICAL)
    if debug_flag:
        logger.setLevel(logging.DEBUG)


def get_files_from_file(fp):
    with open(fp, "r") as f:
        files = f.read().splitlines()
        files = [file for file in files if file != "" and not file.isspace()]

    return files
