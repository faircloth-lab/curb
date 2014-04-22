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
import re
import sys
import yaml
import glob
import shutil
import random
import dendropy
import subprocess
import multiprocessing

from curb import core
from curb.log import setup_logging

import pdb


def get_best_ML_tree(working_dir, raxml, alignment, orig_aln_name, searches=20):
    cmd = [
        raxml,
        "-m",
        "GTRGAMMA",
        "-p",
        str(random.randrange(1,100000000)),
        "-s",
        alignment,
        "-N",
        str(searches),
        "-n",
        "{}.BEST".format(orig_aln_name)
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    best_tree = os.path.join(working_dir, "RAxML_bestTree.{}.BEST".format(orig_aln_name))
    # cleanup temp files
    cleanup_raxml_temp_files(working_dir, orig_aln_name, "BEST")
    # return best tree name
    return best_tree


def cleanup_raxml_temp_files(working_dir, aln_name, postfix):
    for name_stub in [
            "RAxML_log.{}.{}.RUN.*".format(aln_name, postfix),
            "RAxML_parsimonyTree.{}.{}.RUN.*".format(aln_name, postfix),
            "RAxML_result.{}.{}.RUN.*".format(aln_name, postfix)
        ]:
        [os.remove(f) for f in glob.glob(os.path.join(working_dir, name_stub))]


def get_best_constraint_tree(working_dir, raxml, alignment, orig_aln_name, constraint, searches=20):
    constraint_name, constraint_tree_pth = constraint
    cmd = [
        raxml,
        "-g",
        constraint_tree_pth,
        "-m",
        "GTRGAMMA",
        "-p",
        str(random.randrange(1,100000000)),
        "-s",
        alignment,
        "-N",
        str(searches),
        "-n",
        "{}.{}.constraint.BEST".format(orig_aln_name, constraint_name)
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    best_constraint_tree = os.path.join(working_dir, "RAxML_bestTree.{}.{}.constraint.BEST".format(orig_aln_name, constraint_name))
    # cleanup temp files
    cleanup_raxml_temp_files(working_dir, orig_aln_name, "{}.constraint.BEST".format(constraint_name))
    # return best tree name
    return best_constraint_tree


def prune_and_normalize_constraint_tree(working_dir, alignment_taxa, orig_aln_name, constraint):
    constraint_name, constraint_string = constraint
    constraint_tree = dendropy.Tree.get_from_string(constraint_string, schema="newick",  preserve_underscores=True)
    # get taxa in constraint tree
    constraint_taxa = set(constraint_tree.infer_taxa().labels())
    # find which taxa are missing
    intersect = list(constraint_taxa.intersection(alignment_taxa))
    # prune constraint tree to remove missing taxa
    constraint_tree.retain_taxa_with_labels(intersect)
    # set_branch_lengths_equal
    for edge in constraint_tree.preorder_edge_iter():
        if edge.length:
            edge.length = 1.0
    constraint_tree_pth = os.path.join(
        working_dir,
        "{}.{}.constraint.tre".format(orig_aln_name, constraint_name)
    )
    constraint_tree.write_to_path(
            constraint_tree_pth,
            "newick"
        )
    return constraint_name, constraint_tree_pth


def get_merged_constraint_trees(working_dir, orig_aln_name, best_constraint_trees):
    merged_constraint_tree_pth = os.path.join(
        working_dir,
        "{}.MERGED.constraint.tre".format(orig_aln_name)
    )
    treelist = dendropy.TreeList()
    tree_map = {}
    for cnt, constraint_tree in enumerate(best_constraint_trees):
        treelist.read_from_path(constraint_tree, schema="newick", preserve_underscores=True)
        tree_map[cnt] = constraint_tree
    treelist.write_to_path(
        merged_constraint_tree_pth,
        "newick"
    )
    return merged_constraint_tree_pth, tree_map


def get_sh_test_results(working_dir, raxml, alignment, orig_aln_name, best_tree, merged_constraint):
    cmd = [
        raxml,
        "-f",
        "H",
        "-m",
        "GTRGAMMA",
        "-t",
        best_tree,
        "-z",
        merged_constraint,
        "-s",
        alignment,
        "-n",
        "{}.constraints.SHTEST".format(orig_aln_name)
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    sh_test_result = os.path.join(working_dir, "RAxML_info.{}.constraints.SHTEST".format(orig_aln_name))
    # return sh_test results
    return sh_test_result


def filter_sh_test_results(sh_test, tree_map, orig_aln_name):
    regex = re.compile("""
        Tree:\s(\d+)\s              # capture tree num
        Likelihood:\s(-\d+\.\d+)\s  # capture loglik
        \\D\(LH\):\s(-\d+\.\d+)\s   # capture loglike Delta
        SD:\s(\d+\.\d+)\s           # get SD
        Significantly\sWorse\:\s+   # capture test values
        (\w+)\s\((\d+)\%\),\s+(\w+)\s\((\d+)\%\),\s+(\w+)\s\((\d+)\%\)
        """, re.VERBOSE)
    with open(sh_test, 'rU') as infile:
        temp_results = re.findall(regex, infile.read())
    # add actual tree name to results
    results = [(tree_map[int(i[0])],) + i[1:] for i in temp_results]
    return orig_aln_name, results


def get_all_merged_trees(working_dir, orig_aln_name, tree_map, best, merged_constraint_trees):
    all_tree_pth = os.path.join(
        working_dir,
        "{}.ALL.MERGED.tre".format(orig_aln_name)
    )
    treelist = dendropy.TreeList()
    # join merged
    all_trees = [merged_constraint_trees, best]
    for cnt, tree in enumerate(all_trees):
        treelist.read_from_path(tree, schema="newick", preserve_underscores=True)
    # we're inserting the best tree last in the all-merged tree file
    # so add appropriate index
    tree_map[max(tree_map.keys()) + 1] = best
    # set branch lengths = 1
    for tree in treelist:
        for edge in tree.preorder_edge_iter():
            if edge.length != 1.0:
                edge.length = 1.0
    treelist.write_to_path(
        all_tree_pth,
        "newick"
    )
    return all_tree_pth, tree_map


def get_tree_puzzle(working_dir, raxml, alignment, orig_aln_name, all_tree_pth):
    pdb.set_trace()
    cmd = [
        raxml,
        "-f",
        "G",
        "-m",
        "GTRGAMMA",
        "-z",
        all_tree_pth,
        "-s",
        alignment,
        "-n",
        "{}.sitelh".format(orig_aln_name)
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    puzzle_result = os.path.join(working_dir, "RAxML_info.{}.puzzle.SITELH".format(orig_aln_name))
    # return sh_test results
    return puzzle_result


def worker(work):
    # get starting dir
    owd = os.getcwd()
    # unpack work
    args, raxml, orig_alignment, constraints = work
    # make a locus specific working dir
    orig_aln_full_name = os.path.basename(orig_alignment)
    orig_aln_name = os.path.splitext(orig_aln_full_name)[0]
    working_dir = os.path.join(args.output, orig_aln_name)
    os.makedirs(working_dir)
    # copy old aln to new directory
    working_alignment = os.path.join(working_dir, orig_aln_full_name)
    shutil.copyfile(orig_alignment, working_alignment)
    # change to new aln working dir
    os.chdir(working_dir)
    # estimate the best ML tree for the data
    best_tree = get_best_ML_tree(
        working_dir,
        raxml,
        working_alignment,
        orig_aln_name,
        args.searches
    )
    # under each of the constraints, estimate the best constrained ML tree
    taxa_present = core.get_taxa_in_alignment(working_alignment)
    best_constraint_trees = []
    for constraint in constraints.items():
        # pull out missing taxa from constraint tree
        constraint = prune_and_normalize_constraint_tree(
            working_dir,
            taxa_present,
            orig_aln_name,
            constraint
        )
        # feed constraint tree and alignment to raxml
        best_constraint_tree = get_best_constraint_tree(
            working_dir,
            raxml,
            working_alignment,
            orig_aln_name,
            constraint,
            args.searches
        )
        best_constraint_trees.append(best_constraint_tree)
    # get all constraint trees into one file
    merged_constraint_trees, merged_constraint_tree_map = get_merged_constraint_trees(
        working_dir,
        orig_aln_name,
        best_constraint_trees
    )
    # run SH test of unconstrained against constrained
    sh_test = get_sh_test_results(
        working_dir,
        raxml,
        working_alignment,
        orig_aln_name,
        best_tree,
        merged_constraint_trees
    )
    sh_test_results = filter_sh_test_results(
        sh_test,
        merged_constraint_tree_map,
        orig_aln_name
    )
    # prep a site likelihood file
    all_tree_pth, all_tree_map = get_all_merged_trees(
        working_dir,
        orig_aln_name,
        merged_constraint_tree_map,
        best_tree,
        merged_constraint_trees
    )
    get_tree_puzzle(
        working_dir,
        raxml,
        working_alignment,
        orig_aln_name,
        all_tree_pth
    )
    # output per-site likelihoods
    os.chdir(owd)
    # write some progress indicator
    sys.stdout.write(".")
    sys.stdout.flush()
    return sh_test_results


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
    try:
        assert len(valid_alignments) > 0
    except:
        raise IOError("There are not alignments to use.")
    # package up alignments with other params
    work = [(args, raxml, alignment, config["constraints"]) for alignment in alignments]
    # start run
    sys.stdout.write("Running")
    if args.cores > 1:
        assert args.cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
        pool = multiprocessing.Pool(args.cores)
        results = pool.map(worker, work)
        # close the pool.  so sad.
        pool.close()
    else:
        results = map(worker, work)
    pdb.set_trace()
    print ""
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))

