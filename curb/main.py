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
import sys
from curb import core
from curb.log import setup_logging


def main(args):
    core.check_dependencies(args)
    # setup logging
    log, my_name = setup_logging(args)
