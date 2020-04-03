import logging
import argparse
from argparse import Namespace

import click

import bbgmod.stats as st
import bbgmod.logistic_errors as le

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

# ------------------------------------------------------------------------------
# ARGPARSER
# ------------------------------------------------------------------------------

@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    '--debug', help="Show more progress details", is_flag=True)
def cli(debug):
    logging_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=logging_level)

    if not debug:
        # Hide bgdata messages
        logging.getLogger('bgdata').setLevel(logging.WARNING)




@cli.command(short_help='Obtain list of reads that passed q-score')
@click.option(
    '-sp', '--signif-pos', default='', 
    help='Significant positions or motif to restrict the error comparison'
)
@click.option(
    '-te', '--treat-errors', default='', 
    help='Frequency file of the treated sample'
)
@click.option(
    '-ue', '--untreat-errors', default='', 
    help='Frequency file of the untreated sample'
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def error_analysis(treat_errors, untreat_errors, signif_pos, output):
    le.main(treat_errors, untreat_errors, signif_pos, output)



@cli.command(short_help='Level sample compare in house')
@click.option(
    '-f5', '--fast5-basedirs', default='', 
    help='Directory containing fast5 reads of the control'
)
@click.option(
    '-af5', '--alternate-fast5-basedirs', default='',  
    help='Directory containing fast5 reads of the treated'
)
@click.option(
    '-cg', '--corrected-group', default='RawGenomeCorrected_000',  
    help='Corrected group within fast5 after resquiggle command is run'
)
@click.option(
    '-bs', '--basecall-subgroups', default='BaseCalled_template',  
    help='Basecalled group withing fast5'
)
@click.option(
    '-cpu', '--cpus', default=1, help='number of processes to be used, default 1'
)
@click.option(
    '-o', '--output', default='', help='Folder to ouput files'
)
def detect_modifications(fast5_basedirs, alternate_fast5_basedirs, cpus, 
    corrected_group, basecall_subgroups, output):
    st.main(fast5_basedirs, alternate_fast5_basedirs, cpus, 
        corrected_group, basecall_subgroups, output)


if __name__ == '__main__':
    cli()