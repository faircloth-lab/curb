#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 21 April 2014 21:16 PDT (-0700)
"""

import os
import pytest
from sh_t import core

class TestAlignments:
    def test_correct_aligments(self):
        expected = set([
            'uce-10.phylip',
            'uce-13.phylip',
            'uce-19.phylip',
            'uce-3.phylip',
            'uce-7.phylip',
            'uce-8.phylip',
            'uce-508.phylip'
        ])
        test_alignments = os.path.join(
                os.path.dirname(__file__),
                "alignments"
            )
        alignments = core.get_alignments(test_alignments)
        observed = set([os.path.basename(i) for i in alignments])
        assert observed == expected

    def test_taxa_in_align(self):
        test_alignment = os.path.join(
                os.path.dirname(__file__),
                "alignments",
                "uce-10.phylip"
            )
        group_set = {"clupeiforms":[
            "thryssa_hamiltonii2",
            "chirocentrus_dorab2",
            "dorosoma_pentense"
        ]}
        observed = core.satisfy_all_taxon_groups(test_alignment, group_set)
        assert observed is True

    def test_taxa_not_in_align(self):
        test_alignment = os.path.join(
                os.path.dirname(__file__),
                "alignments",
                "uce-508.phylip"
            )
        group_set =  {"gonorynchiforms":[
            "gonorynchus_sp2"
            "chanos_chanos2"
            "parakneria_abbreviata"
        ]}
        with pytest.raises(core.GroupError):
            core.satisfy_all_taxon_groups(test_alignment, group_set)

