#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 21 April 2014 20:48 PDT (-0700)
"""


from __future__ import absolute_import
import os
import sys
import yaml
import shutil
from curb import core
from curb.log import setup_logging

import pdb


def get_best_ML_tree(raxml, alignment):
    pass


def worker(work):
    args, raxml, orig_alignment, constraints = work
    # make a directory to hold the locus info
    orig_aln_full_name = os.path.basename(orig_alignment)
    orig_aln_name = os.path.splitext(orig_aln_full_name)[0]
    working_dir = os.path.join(args.output, orig_aln_name)
    os.makedirs(working_dir)
    # copy old aln to new directory
    shutil.copyfile(orig_alignment, working_dir)
    working_alignment = os.path.join(working_dir, orig_aln_full_name)
    # estimate the best ML tree for the data
    best_tree = get_best_ML_tree(args, raxml, working_alignment)
    # under each of the constraints, estimate the best constrained ML tree
    for each in constraints:
        constraint_trees = get_best_constraint_tree(
            args,
            raxml,
            working_alignment,
            constraints
        )
    # run SH test of unconstrained against constrained
    get_sh_test_results(args, raxml, working_alignment, best_tree, constraint_trees)
    # output per-site liklihoods


def main(args):
    # setup logging
    log, my_name = setup_logging(args)
    log.info("Getting alignments")
    # load yaml
    config = yaml.load(open(args.config))
    # get raxml
    raxml = core.which("raxmlHPC-SSE3")
    # get and check alignments for taxon membership
    alignments = core.get_alignments(args.alignments)
    valid_alignments = set([])
    for alignment in alignments:
        try:
            if core.satisfy_all_taxon_groups(alignment, config["orders"]):
                valid_alignments.add(alignment)
        except core.GroupError, e:
            if e.message == "Not all taxa present in Group":
                log.warn("Dropped {} due to missing taxa from '{}'".format(
                    e.alignment,
                    e.group
                ))
            else:
                raise core.GroupError(e)
    # convert remaining alignments to a list
    valid_alignments = list(valid_alignments)
    # package up alignments with other params
    work = [(args, raxml, alignment, config["constraints"]) for alignment in alignments]
    map(worker, work)
