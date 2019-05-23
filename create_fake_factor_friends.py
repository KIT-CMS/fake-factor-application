#!/usr/bin/env python
# -*- coding: utf-8 -*-

from python.calculate_fake_factors import *

import os
import yaml
import ROOT
import numpy
import copy
from array import array
from multiprocessing import Pool

import argparse
import logging
logger = logging.getLogger()


def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate fake factors and create friend trees.")
    parser.add_argument(
        "-i",
        "--input-directory",
        required=True,
        type=str,
        help="Directory with Artus outputs.")
    parser.add_argument(
        "-o",
        "--output-directory",
        required=True,
        type=str,
        help="Directory to write outputs.")
    parser.add_argument(
        "--et-friend-directories",
        nargs='+',
        default=[],
        type=str,
        help=
        "Directory arranged as Artus output and containing a friend tree for et."
    )
    parser.add_argument(
        "--mt-friend-directories",
        nargs='+',
        default=[],
        type=str,
        help=
        "Directory arranged as Artus output and containing a friend tree for mt."
    )
    parser.add_argument(
        "--tt-friend-directories",
        nargs='+',
        default=[],
        type=str,
        help=
        "Directory arranged as Artus output and containing a friend tree for tt."
    )
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        type=str,
        help="Key of desired config in yaml config file.")
    parser.add_argument(
        "-C",
        "--configpath",
        default="fake-factor-application/config.yaml",
        type=str,
        help="Path to yaml config file.")
    parser.add_argument(
        "--et-fake-factor-directory",
        required=True,
        type=str,
        help="Directory containing et fake factor inputs.")
    parser.add_argument(
        "--mt-fake-factor-directory",
        required=True,
        type=str,
        help="Directory containing mt fake factor inputs.")
    parser.add_argument(
        "--tt-fake-factor-directory",
        required=True,
        type=str,
        help="Directory containing tt fake factor inputs.")
    parser.add_argument(
        "--era", type=str, required=True, help="Experiment era.")
    parser.add_argument(
        "--num-threads",
        default=32,
        type=int,
        help="Number of threads to be used.")
    parser.add_argument(
        "--category-mode", type=str, help="Category mode. If 'inclusive' fake factors are calculated inclusively, otherwise depending on NN categories")
    parser.add_argument(
        "-w",
        "--fractions-from-worspace",
        action="store_true",
        default=False,
        help="Use fractions from workspace."
    )
    parser.add_argument(
        "--workspace", type=str, help="Path to workspace for fractions if -w option is set")
    return parser.parse_args()


def apply_fake_factors_arg_wrapper(config):
    apply_fake_factors(**config)

def main(args):
    if args.category_mode=="inclusive":
        logger.warning("Option to calculate fake factors inclusively has been set. No categorization applied!")
        
    fractions = load_fractions(args.configpath, args.config, args.fractions_from_worspace, args.workspace, "fake-factor-application/%s_ff_yields.root" % args.era, args.era)

    #find paths to data files the fake factors are appended to
    datafiles = []
    for entry in os.listdir(args.input_directory):
        if not("HToTauTau" in entry
                or "DoubleEG" in entry
                or "DoubleMuon" in entry
                or "MuonEG" in entry
                or "WJets" in entry
                or "W1Jets" in entry
                or "W2Jets" in entry
                or "W3Jets" in entry
                or "W4Jets" in entry):
            path = os.path.join(entry, "%s.root" % entry)

            #check whether expected files exist
            if not os.path.isfile(os.path.join(args.input_directory, path)):
                logger.critical(
                    "Expected file %s does not exist. Check --input-directory option!"
                )
                raise Exception
            
            datafiles.append(path)

    #create output directory
    if os.path.exists(args.output_directory):
        logger.critical(
            "Output directory %s already exists. I don't want to overwrite it!"
            % args.output_directory)
        raise Exception
    else:
        os.mkdir(args.output_directory)
        logger.info("Createt output directory")

    logger.info("Create friend trees...")

    pool = Pool(processes=args.num_threads)
    pool.map(apply_fake_factors_arg_wrapper, [{
            "datafile" : os.path.join(args.input_directory, datafile),
            "friendfilelists" : {
                "et" : [os.path.join(frienddir, datafile) for frienddir in args.et_friend_directories],
                "mt" : [os.path.join(frienddir, datafile) for frienddir in args.mt_friend_directories],
                "tt" : [os.path.join(frienddir, datafile) for frienddir in args.tt_friend_directories]
                },
            "outputfile" : os.path.join(args.output_directory, datafile),
            "category_mode" : args.category_mode,
            "fakefactordirectories" : {
                "et" : args.et_fake_factor_directory,
                "mt" : args.mt_fake_factor_directory,
                "tt" : args.tt_fake_factor_directory
                },
            "fractions" : fractions,
            "use_fractions_from_worspace" : args.fractions_from_worspace,
            "configpath" : args.configpath,
            "expression" : args.config,
            "era" : args.era
        } for datafile in datafiles])
    pool.close()
    pool.join()
    del pool


if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("{}_calculate_fake_factor_friends.log".format(args.era),
                  logging.INFO)
    main(args)
