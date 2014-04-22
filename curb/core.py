#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 21 April 2014 20:54 PDT (-0700)
"""


import os
import sys
import glob
import shutil
import argparse
import subprocess

from Bio import AlignIO

import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class CreateDir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # get the full path
        d = os.path.abspath(os.path.expanduser(values))
        # check to see if directory exists
        if os.path.exists(d):
            answer = raw_input("[WARNING] Output directory exists, REMOVE [Y/n]? ")
            if answer == "Y":
                shutil.rmtree(d)
            else:
                print "[QUIT]"
                sys.exit()
        # create the new directory
        os.makedirs(d)
        # return the full path
        setattr(namespace, self.dest, d)


class GroupError(Exception):
    def __init__(self, message, group, alignment):
        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)
        # Now for your custom code...
        self.group = group
        self.alignment = alignment


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def which(prog):
    cmd = ["which", prog]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    if stderr:
        raise EnvironmentError("Program {} does not appear to be installed")
    else:
        return stdout.strip()


def get_alignments(alignment_dir):
    alignments = []
    for ftype in ('.phylip', '.phy'):
        alignments.extend(glob.glob(os.path.join(alignment_dir, "*{}".format(ftype))))
    return alignments


def satisfy_one_taxon_group(taxa_in_align, taxon_group):
    try:
        isinstance(taxon_group, list)
    except:
        raise AssertionError("Taxon group is not a list.")
    group_set = set(taxon_group)
    # ensure there is at least one member in each group
    if len(taxa_in_align.intersection(group_set)) >= 1:
        return True
    else:
        return False


def get_taxa_in_alignment(alignment):
    aln = AlignIO.read(alignment, "phylip-relaxed")
    taxa_in_align = set([taxon.id for taxon in aln])
    return taxa_in_align


def satisfy_all_taxon_groups(alignment, taxon_groups):
    """given an input alignment, see if any taxa in list are in file"""
    taxa_in_align = get_taxa_in_alignment(alignment)
    taxa_present = []
    for group_name, taxon_group in taxon_groups.iteritems():
        if satisfy_one_taxon_group(taxa_in_align, taxon_group):
            taxa_present.append(True)
        else:
            taxa_present.append(False)
    if all(taxa_present):
        return True
    else:
        raise GroupError(
            "Not all taxa present in Group",
            group_name,
            os.path.basename(alignment),
        )
