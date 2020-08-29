#coding=utf-8
#!/usr/bin/env python

""" load 3D nPE-Map from guihong """

from ROOT import TFile, TH2D

import global_

def map_loader(object):

    def __init__(self, filename):
        self.filename = filename
