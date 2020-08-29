#coding=utf-8
#!/usr/bin/env python

""" prediciont nonlinearity mode module """
""" electron nonlinearity model prediction  """

import uproot as up
import numpy as np

import global_param as glo

class electron_nonlinearity(object):
    """ electron nonlinearity prediction """

    quenchNL = []
    cerNL = []
    
    def __init__(self):
        pass


    def read_quench(self):
        kB = glo.get_value("kB")
        filename = glo.get_value("QuenchFile")
        try:
            f = up.open(filename)
            hist_name = "kB%d"%kB
            hist = f[hist_name]
            cont, edge =  hist.numpy()
            self.quenchNL = cont.tolist()
        except IOError:
            print(" >>> Can not find electron quenchign file")

        #print(" Test: energy {:.2f} --> bin: {:d} and bincenter: {: .2f} and bincontent: {: .2f}".format(Etrue, idx, edge[idx], cont[idx]))


    def read_cerenkov(self):
        filename = glo.get_value("CerenkovFile")
        with open(filename) as f:
            for lines in f.readlines():
                line = lines.strip("\n")
                data = line.split(" ")
                self.cerNL.append(float(data[2]))

    
    def initialize(self):
        self.read_quench()
        self.read_cerenkov()


    def quench_NL(self, Etrue):
        """ linear interpolation"""
        idx1 = int(Etrue/0.001) if Etrue<=0.1 else int((Etrue-0.1)/0.01)+100
        idx2 = idx1+1
        frac = (Etrue-0.001*idx1)/0.001 if Etrue<=0.1 else (Etrue-0.1-0.01*(idx1-100))/0.01
        return self.quenchNL[idx1]*(1-frac) + self.quenchNL[idx2]*frac



    def cerenkov_NL(self, Etrue):
        """ linear interpolation"""
        idx1 = int(Etrue/0.01)
        idx2 = idx1+1
        frac = (Etrue-0.01*idx1)/0.01
        return (self.cerNL[idx1])/Etrue/glo.get_value("energy_scale")*(1-frac) + (self.cerNL[idx2])/Etrue/glo.get_value("energy_scale")*frac


    def total_NL(self, Etrue):
        kA = glo.get_value("kA")
        kC = glo.get_value("kC")
        return  kA*self.quench_NL(Etrue) + kC*self.cerenkov_NL(Etrue)

