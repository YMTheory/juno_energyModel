#coding=utf-8
#!/usr/bin/env python

""" prediciont mode module """
""" gamma model prediction  """

import math
import uproot as up
import random
import numpy as np

import global_param as glo
import electron_nonlinearity as en
import electron_resolution as er


class gamma_prediction(object):
    """ Gamma Resolution Prediction class"""

    # clsss variable: 

    """ option1: 2-layer sampling method """
    __mean_dist  = []
    __sigma_dist = []
    """ option2: average distribution method  """
    __prim_etru = []
    __prim_prob = []
    elec_nonl = en.electron_nonlinearity()

    def initialize(self):
        self.__mean_dist.clear()
        self.__sigma_dist.clear()

    def __init__(self, name, energy):
        """ gamma source name ane energy are required to offer """
        self.name = name
        self.energy = energy  # unit MeV

    def print_name(self):
        print(" Get Gamma Source Name : {} with gamma energy {: .3f} MeV". format(self.name, self.energy))


    def get_name(self):
        return self.name

    def get_energy(self):
        return self.energy


    # read primary electron distribution: option 1
    def read_primTXT(self):
        self.elec_nonl.initialize()
        filename = glo.get_dict_value("GammaPrimElecFile", self.name)
        with open(filename, 'r') as f:
            tmp_mean, tmp_sigma = 0, 0
            for lines in f.readlines():
                line = lines.strip("\n")
                data = line.split(" ")
                while '' in data:
                    data.remove('')
                data = list(map(float, data))
                tmp_mean = sum([self.elec_nonl.total_NL(etru)*etru for etru in data])
                self.__mean_dist.append(tmp_mean)
                tmp_sigma = sum([er.electron_sigma2_total_option2(etru) for etru in data])
                self.__sigma_dist.append(math.sqrt(tmp_sigma))


    def predict_option1(self):  # do a two-layer sampling 
        self.initialize()
        self.read_primTXT()
        sampleid = [ int( random.uniform(0, 5000) ) for i in range(50000)]  # first sampling 
        sample_npe = [ random.gauss(self.__mean_dist[idx], self.__sigma_dist[idx]) for idx in sampleid ] # second sampling
        evis = np.mean(np.array(sample_npe))
        width = np.std(np.array(sample_npe))
        #print(" {} prediction evis {:.4f} and resolution {: .4f}" .format(self.name, evis, width/evis) )
        return evis, width/evis



    # prediction: option 2
    def read_primROOT(self):
        filename = glo.get_value("GammaPrimElecDistFile")
        ff = up.open(filename)
        graph_name = "gamma"+self.name
        graph = ff[graph_name]
        for e, p in zip(graph._fX, graph._fY):
            self.__prim_etru.append(e)
            self.__prim_prob.append(p)

    
    def predict_option2(self):  # average distribution method
        self.elec_nonl.initialize()
        self.initialize()
        self.read_primROOT()
        part1, part2 = 0, 0
        for e, p in zip(self.__prim_etru, self.__prim_prob):
            part1 += e*p*self.elec_nonl.total_NL(e)
            part2 += p*e
        try:
            part1 / part2
        except ValueError:
            print("Error in nonlinearity")
        return part1/part2


