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

################################################################################
#   Imports
################################################################################
import sys, os
import configparser
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from ruffus import *
import subprocess
import logging

exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))

dirpath = os.getcwd()
################################################################################
#   Argument parser definitions
################################################################################
parser = ArgumentParser(prog="PhyPipe",
                        formatter_class=RawDescriptionHelpFormatter,
                        description=description, epilog=epilog)
parser.add_argument("-i", "--input_file", dest="input_file",
                  metavar="FASTA_FILE", required = True,
                  nargs='+',
                  type=str,
                  help="Name and path of input file(s) in FASTA format. ")
parser.add_argument("-c", "--config_file", dest="config_file", required = True,
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
parser.add_argument("-v", "--verbose", dest = "verbose",
                  action="count", default=1,
                  help="Print more detailed messages for each additional verbose level."
                       " E.g. run_parallel_blast --verbose --verbose --verbose ... (or -vvv)")
parser.add_argument("-t", "--threads", dest="jobs",
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
# get basic parameters from command line
config_file = options.config_file
input_file = options.input_file
working_dir = options.output_directory
result_file = options.result_file
################################################################################
#   Helper functions
################################################################################
def run_cmd(cmd_str):
    """
    This helper function excecute a bash command as a separate process
    and returns an exception if the command fails.

    It also returns the standard output and the standard error of the command as
    str objects.

    """
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str, process.returncode))
    return stdout_str.decode(), stderr_str.decode()

def model_parser2garli(model_string):
    """
    Parses a model string to the equivalent command block for Garli
    e.g. converts "GTR+I+G" to: "ratematrix = 6rate
                                 statefrequencies = estimate
                                 ratehetmodel = gamma
                                 numratecats = 4
                                 invariantsites = estimate"

    """
    #define equivalencies in dictionaries
    model_dict = {"GTR":"ratematrix = 6rate\nstatefrequencies = estimate\n",
                  "SYM":"ratematrix = 6rate\nstatefrequencies = equal\n",
                  "HKY":"ratematrix = 2rate\nstatefrequencies = estimate\n",
                  "K80":"ratematrix = 2rate\nstatefrequencies = equal\n",
                  "F81":"ratematrix = 1rate\nstatefrequencies = estimate\n",
                  "JC":"ratematrix = 1rate\nstatefrequencies = equal\n"}
    gamma_dict = {True:"ratehetmodel = gamma\nnumratecats = 4\n",
                    False:"ratehetmodel = none\nnumratecats = 1\n"}
    inv_dict = {True:"invariantsites = estimate\n",
                False:"invariantsites = none"}
    # parse model name
    model_info = model_string.split("+")
    model = model_info[0]
    inv_status = True if "I" in model_info else False
    gamma_status = True if "G" in model_info else False
    # compose Garli block
    model_string = model_dict[model] if model in model_dict else "GTR"
    composed_string = model_string + gamma_dict[gamma_status] + \
                      inv_dict[inv_status]
    return composed_string

################################################################################
#   logging
################################################################################
logger = logging.getLogger("run_phypipe")

if options.verbose:
    logger.setLevel(logging.DEBUG)
    stderrhandler = logging.StreamHandler(sys.stderr)
    stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
    stderrhandler.setLevel(logging.DEBUG)
    logger.addHandler(stderrhandler)

################################################################################
#   Get run parameters from config file
################################################################################
config_handle = configparser.ConfigParser()
config_handle.read(config_file)

# provisional parameters
garli_params = config_handle['Garli run parameters']
searchreps = garli_params['searchreps']
bootstrapreps = garli_params['bootstrapreps']

print(searchreps,bootstrapreps)
################################################################################
#   General pipelines definitions
################################################################################
phypipe_single_locus = Pipeline(name = "Single-locus analysis")
phypipe_multi_locus = Pipeline(name = "Multi-locus analysis")
################################################################################
#   Tasks definitions
################################################################################
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

def garliconf_gen(inputs, output, searchreps, bootstrapreps):
    """
    """
    base_conf = "[general]\ndatafname = {0}\nconstraintfile = none\nstreefname = stepwise\nattachmentspertaxon = 50\nofprefix = {1}\nrandseed = -1\navailablememory = 512\nlogevery = 10\nsaveevery = 100\nrefinestart = 1\noutputeachbettertopology = 0\noutputcurrentbesttopology = 0\nenforcetermconditions = 1\ngenthreshfortopoterm = 20000\nscorethreshforterm = 0.05\nsignificanttopochange = 0.01\noutputphyliptree = 0\noutputmostlyuselessfiles = 0\nwritecheckpoints = 0\nrestart = 0\noutgroup = 1\nresampleproportion = 1.0\ninferinternalstateprobs = 0\noutputsitelikelihoods = 0\noptimizeinputonly = 0\ncollapsebranches = 1\n\nsearchreps = {2}\nbootstrapreps = {3}\n\n[model1]\n{4}\n\n"
    master = "[master]\nnindivs = 4\nholdover = 1\nselectionintensity = 0.5\nholdoverpenalty = 0\nstopgen = 5000000\nstoptime = 5000000\n\nstartoptprec = 0.5\nminoptprec = 0.01\nnumberofprecreductions = 10\ntreerejectionthreshold = 50.0\ntopoweight = 1.0\nmodweight = 0.05\nbrlenweight = 0.2\nrandnniweight = 0.1\nrandsprweight = 0.3\nlimsprweight =  0.6\nintervallength = 100\nintervalstostore = 5\nlimsprrange = 6\nmeanbrlenmuts = 5\ngammashapebrlen = 1000\ngammashapemodel = 1000\nuniqueswapbias = 0.1\ndistanceswapbias = 1.0"
    print(inputs)
    nexus_name, input_model = inputs
    with open(output, "a") as conf_file, open(input_model) as model_file:
        model_block = model_parser2garli(model_file.read().split("\t")[1])
        prefix = output.partition("_garli.conf")[0]
        conf_string = base_conf.format(nexus_name, prefix, searchreps, bootstrapreps, model_block) + \
                master
        conf_file.write(conf_string)

def run_garli(conf_file, outfile):
    run_cmd("Garli {}".format(conf_file))
    open(outfile)

def sumtrees(input_trees, sum_tree):
    best_tree, boot_trees = input_trees
    run_cmd("sumtrees.py -p -d0 -o {} -t {} {}".format(sum_tree, best_tree, boot_trees))
    open(sum_tree)

################################################################################
# Pipeline definition for single-locus phypipe
################################################################################
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
phypipe_single_locus.merge(task_func = garliconf_gen,
                            name = "ML_garliconf_gen",
                            input = [output_from("fasta2nexus"), output_from("modeltest")],
                            output = working_dir + "/ML_garli.conf",
                            extras = [searchreps, "0"])
phypipe_single_locus.merge(task_func = garliconf_gen,
                            name = "BS_garliconf_gen",
                            input = [output_from("fasta2nexus"), output_from("modeltest")],
                            output =  working_dir + "/BS_garli.conf",
                            extras = ["1", bootstrapreps])
phypipe_single_locus.transform(task_func = run_garli,
                                name = "garli_ML",
                                input = output_from("ML_garliconf_gen"),
                                filter = suffix("_garli.conf"),
                                output = ".best.tre",
                                output_dir = working_dir)
phypipe_single_locus.transform(task_func = run_garli,
                                name = "garli_BS",
                                input = output_from("BS_garliconf_gen"),
                                filter = suffix("_garli.conf"),
                                output = ".boot.tre",
                                output_dir = working_dir)
phypipe_single_locus.merge(task_func=sumtrees,
                            input = [output_from("garli_ML"), output_from("garli_BS")],
                            output = "ML_w_bootstrap.tre")

################################################################################
#   Pipeline selection according number of files found in input (-i FILE(s))
################################################################################
if input_file:
    num_files = len(input_file)
    print("Number of files detected as input: {}".format(num_files))
    for filename in input_file :
        print(filename)
else:
    num_files = 0

if num_files == 1:
    phypipe = phypipe_single_locus
    print("\nPhyPipe is going to be executed in single-locus mode.\n")
else:
    phypipe = phypipe_multi_locus
    print("\nPhyPipe is going to be executed in multi-locus mode.\n")
    print("Not yet fully implemented :/")

print("Results are going to be written in:\n {}".format(dirpath + "/" + working_dir))

################################################################################
#   print list of tasks for the selected workflow (only if -n option is given)
################################################################################
if options.just_print:
    phypipe.printout(sys.stdout, verbose=options.verbose)
################################################################################
#   print flowchart for selected pipeline (if --flowchart FILE is given)
################################################################################
elif options.flowchart:
    output_format = os.path.splitext(options.flowchart)[1][1:]
    phypipe.printout_graph(stream = options.flowchart,
                             output_format = output_format,
                             no_key_legend = True)
################################################################################
#   Run selected pipeline (if -n is not present)
################################################################################
else:
    phypipe.run(multiprocess = options.jobs,
                        logger = logger, verbose=options.verbose)
