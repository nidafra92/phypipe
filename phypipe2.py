#!/usr/bin/env python3
version = "2.0" #parallel features
description = \
"""
    PhyPipe2: Automated pipeline for phylogenetic reconstruction. \n
    Version {}
    coded by Nicolás D. Franco-Sierra

    Expected inputs:
    (1) Nucleotide sequences of a desired marker from several taxa in
     FASTA format (more than one file is supported) (-i/--input_file).
    (2) Configuration file listing all options for all software in the
     routine (-c/--config_file).

    Output: final phylogenetic reconstructed by the desired method.

    What it does? it runs a typical phylogenetic reconstruction routine
    combining all the required software for all the involved steps (i.e.
    alignment, model selection, phylogenetic reconstruction, topological
    tests)

""".format(version)

epilog = \
"""
    If only one FASTA file is specified as input, "single-locus routine" will
    be executed. If multiple FASTA files are given PhyPipe will be executed
    in its traditional multi-locus mode.
"""

import sys, os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from Bio import SeqIO
from Bio.Alphabet import IUPAC

exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))

dirpath = os.getcwd()

parser = ArgumentParser(prog="PhyPipe",
                        formatter_class=RawDescriptionHelpFormatter,
                        description=description, epilog=epilog)
parser.add_argument("-i", "--input_file", dest="input_file",
                  metavar="FASTA_FILE",
                  nargs='+',
                  type=str,
                  help="Name and path of input file(s) in FASTA format. ")
parser.add_argument("-c", "--config_file", dest="config_file",
                  metavar="CONFIG_FILE",
                  type=str,
                  help="Name and path of config file to excecute. ")
parser.add_argument("--result_file", dest="result_file",
                  metavar="RESULT_FILE",
                  type=str,
                  default="final.phylo_results",
                  help="Name and path of where the files should end up. ")
parser.add_argument("-d", "--output_directory", dest="output_directory",
                  metavar="OUTPUT_PATH",
                  type=str,
                  default="phypipe_results",
                  help="Name and path of output directory where calculations"
                            "should take place. ")

parser.add_argument('--version', action='version', version='%(prog)s {}'.format(version))
# verbosity
parser.add_argument("-v", "--verbose", dest = "verbose",
                  action="count", default=0,
                  help="Print more detailed messages for each additional verbose level."
                       " E.g. run_parallel_blast --verbose --verbose --verbose ... (or -vvv)")
#   pipeline
#
parser.add_argument("-j", "--jobs", dest="jobs",
                  default=1,
                  metavar="THREADS",
                  type=int,
                  help="Specifies the number of jobs (operations) to run \
                  in parallel.")
parser.add_argument("--flowchart", dest="flowchart",
                  metavar="FILE",
                  type=str,
                  help="Print flowchart of the pipeline to FILE. Flowchart "
                       "format depends on extension. Alternatives include "
                       "('.dot', '.jpg', '*.svg', '*.png' etc). Formats "
                       "other than '.dot' require "
                       "the dot program to be installed (http://www.graphviz.org/).")
parser.add_argument("-n", "--just_print", dest="just_print",
                    action="store_true", default=False,
                    help="Only print a trace (description) of the pipeline. "
                         " The level of detail is set by --verbose.")

options = parser.parse_args()


if not options.flowchart:
    if not options.config_file:
        parser.error("\n\n\tMissing parameter --config_file FILE\n\n")
    if not options.input_file:
        parser.error("\n\n\tMissing parameter --input_file FILE(s)\n\n")

from ruffus import *
import subprocess

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Functions
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

def run_cmd(cmd_str):
    """
    Throw exception if run command fails
    """
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str, process.returncode))
    return stdout_str.decode(), stderr_str.decode()

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Logger
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

import logging
logger = logging.getLogger("run_phypipe")
#
# We are interesting in all messages
#
if options.verbose:
    logger.setLevel(logging.DEBUG)
    stderrhandler = logging.StreamHandler(sys.stderr)
    stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
    stderrhandler.setLevel(logging.DEBUG)
    logger.addHandler(stderrhandler)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Pipeline tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

input_file = options.input_file
config_file = options.config_file
working_dir = options.output_directory
result_file = options.result_file

if input_file:
    num_files = len(input_file)
    print("Number of files detected as input: {}".format(num_files))
    for filename in input_file :
        print(filename)
else:
    num_files = 0

#print(input_file)
#print(type(input_file))

phypipe_single_locus = Pipeline(name = "Single-locus analysis")
phypipe_multi_locus = Pipeline(name = "Multi-locus analysis")

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Print list of tasks

def task_finished():
    print("Task finished")

def align(input_seq, out_files):
    """
    """
    print("\n [Step: Alignment] \n MAFFT is going to be used as selected aligner.\n")
    run_cmd("mafft --auto {} > {}".format(input_seq, out_files))
    print("File written in: {}".format(out_files))

def modeltest(input_seq, outreports):
    """
    """
    output, error = run_cmd("jModelTest.jar -d {} -s 5 -f -i -g 4 -BIC -p -tr 4".format(input_seq))
    best_model_stats = output.partition("::Best Models::")[2].partition("---\n")[2].split("\t")
    criterion, best_model = best_model_stats[0:2]
    with open(outreports,"a") as model_file:
        model_file.write("{}\t{}".format(criterion, best_model))
    print(criterion, best_model)

def fasta2nexus(input_fasta, output_nexus):
    """
    """
    with open(input_fasta) as fasta, open(output_nexus,"a") as nexus:
        sequences = SeqIO.parse(fasta, "fasta", alphabet=IUPAC.ambiguous_dna)
        SeqIO.write(sequences, nexus, "nexus")

###################################

# Single-locus phypipe
#phypipe_single_locus.mkdir(working_dir)
phypipe_single_locus.transform(task_func = align,
                                input    = input_file,
                                filter   = suffix('.fasta'),
                                output   = r'\1_aligned.fasta',
                                output_dir = working_dir)\
                    .posttask(task_finished)\
                    .mkdir(working_dir)
phypipe_single_locus.transform(task_func = modeltest,
                                input = output_from("align"),
                                filter = suffix("_aligned.fasta"),
                                output = r'\1_best_model.txt',
                                output_dir = working_dir)\
                    .posttask(task_finished)
phypipe_single_locus.transform(task_func = fasta2nexus,
                                input = output_from("align"),
                                filter = suffix('.fasta'),
                                output =  r'\1.nexus',
                                output_dir = working_dir)\
                    .posttask(task_finished)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

if num_files == 1:
    phypipe = phypipe_single_locus
    print("\nPhyPipe is going to be executed in single-locus mode.\n")
else:
    phypipe = phypipe_multi_locus
    print("\nPhyPipe is going to be executed in multi-locus mode.\n")
    print("Not yet fully implemented :/")

print("Results are going to be written in:\n {}".format(dirpath + "/" + working_dir))

if options.just_print:
    phypipe.printout(sys.stdout, verbose=options.verbose)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print flowchart
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
elif options.flowchart:
    # use file extension for output format
    output_format = os.path.splitext(options.flowchart)[1][1:]
    phypipe.printout_graph (open(options.flowchart, "w"),
                             output_format,
                             no_key_legend = True)
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Run Pipeline
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
else:
    phypipe.run(multiprocess = options.jobs,
                        logger = logger, verbose=options.verbose)
