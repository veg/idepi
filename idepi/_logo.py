#
# idepi :: (IDentify EPItope) python libraries containing some useful machine
# learning interfaces for regression and discrete analysis (including
# cross-validation, grid-search, and maximum-relevance/mRMR feature selection)
# and utilities to help identify neutralizing antibody epitopes via machine
# learning.
#
# Copyright (C) 2011 N Lance Hepler <nlhepler@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from math import floor
from operator import itemgetter
from re import match

import PIL as pil

from _util import alph_to_base_26, base_n_to_10, BASE_ALPH

__all__ = ['AlphIndex', 'Logo']

__FONTPATH = '/Users/Lance/Library/Fonts/Inconsolata.otf'

__COLORS = {
# ASP,GLU bright red
    'D': 'rgb(230,10,10)',
    'E': 'rgb(230,10,10)',
# CYS,MET yellow
    'C': 'rgb(230,230,0)',
    'M': 'rgb(230,230,0)',
# LYS,ARG blue
    'K': 'rgb(20,90,255)',
    'R': 'rgb(20,90,255)',
# SER,THR orange
    'S': 'rgb(250,150,0)',
    'T': 'rgb(250,150,0)',
# PHE,TYR mid blue
    'F': 'rgb(50,50,170)',
    'Y': 'rgb(50,50,170)',
# ASN,GLN cyan
    'N': 'rgb(0,220,220)',
    'Q': 'rgb(0,220,220)',
# GLY light grey
    'G': 'rgb(235,235,235)',
# LEU,VAL,ILE green
    'L': 'rgb(15,130,15)',
    'V': 'rgb(15,130,15)',
    'I': 'rgb(15,130,15)',
# ALA dark grey
    'A': 'rgb(200,200,200)',
# TRP pink
    'W': 'rgb(180,90,180)',
# HIS pale blue
    'H': 'rgb(130,130,210)',
# PRO flesh
    'P': 'rgb(220,150,130)',
# UNK black
    'X': 'rgb(0,0,0)'
}


class AlphIndex(object):

    def __init__(self, k):
        m = match(r'^(\d+)([a-z]*)', k)
        if m:
            self.str = k
            self.n = int(m.group(1))
            self.a = base_n_to_10(alph_to_base_26(m.group(2)), BASE_ALPH)
        else:
            raise ValueError

    def __lt__(self, other):
        return self.n < other.n or (self.n == other.n and self.a < other.a)

    def __gt__(self, other):
        return self.n > other.n or (self.n == other.n and self.a > other.a)

    def __eq__(self, other):
        return self.n == other.n and self.a == other.a

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return self.str


class Logo(object):

    def __init__(self, refseq, statistics, folds):
        self.refseq = refseq
        self.statistics = statistics
        self.folds = folds
        self.seq = None
        self.margin = 40
        self.col_height = 0
        self.col_width = 0
        self.image = None
        self.draw = None
        self.font = None

    def __len__(self):
        len(Logo.sequence(self))

    def sequence(self):
        if self.seq is not None:
            return self.seq

        ret = dict([('%d' % (i+1), [(self.refseq[i], -1, -1)]) for i in xrange(len(self.refseq))])
        inserts = {}
        for pos, stats in self.statistics.items():
            m = match(r'^(\d+)([a-z]+)([A-Z])', pos)
            if m:
                # not -1 because insert() inserts things before a position
                ins = m.group(1) + m.group(2)
                if ins not in inserts:
                    ret[ins] = ['']
                ret[ins].append((m.group(3), sum([i for i in stats if i > 0]), sum([i for i in stats if i < 0])))
            else:
                m = match(r'^[A-Z](\d+)([A-Z])', pos)
                if m:
                    assert(m.group(1) in ret)
                    ret[m.group(1)].append((m.group(2), sum([i for i in stats if i > 0]), sum([i for i in stats if i < 0])))
        self.seq = ret
        return ret

    def draw_column(self, x_off, idx):
        pos = {}
        neg = {}
        ref_letter = None
        for letter, p, n in Logo.sequence(self)[idx]:
            letter = letter.upper()
            if letter not in __COLORS:
                letter = 'X'
            if p > 0:
                pos[letter] = p
            if n > 0:
                neg[letter] = n
            if p < 0 and n < 0:
                ref_letter = letter
        assert(ref_letter is not None)
        # draw positive side
        col_height = floor(2. * sum(pos.values()) / self.folds)
        y_off = self.col_height - col_height + self.margin
        for letter, val in sorted(pos, key=itemgetter(1)):
            width, height = self.font.getsize(letter)
            l = pil.Image.new('RGBA', (width, height))
            d = pil.ImageDraw.Draw(l)
            d.setfont(self.font)
            d.text((0, 0), letter, fill=__COLORS[letter])
            height = floor(2. * height * val / self.folds)
            l.resize((width, height), pil.Image.ANTIALIAS)
            self.image.paste(l, (x_off, y_off))
            y_off += height
        # draw central letter
        self.draw.text((x_off, 2 * self.col_height / 5 + self.margin), ref_letter, fill=__COLORS[ref_letter])
        y_off = 3 * col_height / 5 + self.margin
        for letter, val in sorted(pos, key=itemgetter(1), reverse=True):
            width, height = self.font.getsize(letter)
            l = pil.Image.new('RGBA', (width, height))
            d = pil.ImageDraw.Draw(l)
            d.setfont(self.font)
            d.text((0, 0), letter, fill=__COLORS[letter])
            height = floor(2. * height * val / self.folds)
            l.resize((width, height), pil.Image.ANTIALIAS)
            self.image.paste(l, (x_off, y_off))
            y_off += height

    def save(self, filename):
        font = pil.ImageFont.truetype(__FONTPATH, 20)
        char_width, char_height = font.getsize('A')
        # length of the reference sequence plus any insertions, times a single characters width (use monospace fonts), with 40px padding
        len_ = len(self)
        self.margin = 40
        self.col_width = char_width
        self.col_height = 5 * char_height
        size = len_ * self.col_width + 2 * self.margin, self.col_height + 2 * self.margin
        self.image = pil.Image.new('RGBA', size)
        self.draw = pil.ImageDraw.Draw(self.image)
        self.font = font
        indices = sorted(Logo.sequence(self).items(), key=AlphIndex)
        for i in xrange(len(indices)):
            Logo.draw_column(self, char_width * i + self.margin, indices[i])
        self.image.save(filename)
