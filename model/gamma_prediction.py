#coding=utf-8
#!/usr/bin/env python

""" prediciont mode module """
""" gamma resolution model prediction  """

import math
import uproot as up

import global_param as glo

class gamma_resolution(object):
    """ Gamma Resolution Prediction class"""

    # clsss variable: 
    mean_dist  = []
    sigma_dist = []

    def __init__(self, name, energy):
        """ gamma source name ane energy are required to offer """
        self.name = name
        self.energy = energy  # unit MeV

    def print_name(self):
        print(" Get Gamma Source Name : {} with gamma energy {: .3f} MeV". format(self.name, self.energy))

    def read_primTXT(self):
        filename = glo.get_dict_value("GammaPrimElecFile", "Cs137")
        with open(filename, 'r') as f:
            for line in f.readlines():
                tmp_mean, tmp_sigma = 0, 0
                for data in line.split(" "):
                    pass
