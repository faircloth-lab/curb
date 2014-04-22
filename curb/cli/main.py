#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 31 January 2014 11:38 PST (-0800)
"""

from __future__ import absolute_import
#import os
import argparse
from curb import core


def get_args():
    parser = argparse.ArgumentParser(
        description="Run constraint tests using RAxML",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--config',
        required=True,
        help='Path to a YAML config file with constraint info.',
        action=core.FullPaths
    )
    parser.add_argument(
        '--alignments',
        required=True,
        help='A directory of PHYLIP alignments to use for tree inference.',
        action=core.FullPaths,
        type=core.is_dir
    )
    parser.add_argument(
        '--output',
        required=True,
        help='The output directory for results.',
        action=core.CreateDir
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='The number of compute cores to use.',
    )
    parser.add_argument(
        '--searches',
        type=int,
        default=20,
        help='The number of RAxML search reps to use.',
    )
    return parser.parse_args()


def main():
    from curb.main import main
    args = get_args()
    main(args)
