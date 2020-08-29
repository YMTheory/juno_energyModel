#coding=utf-8
#!/usr/bin/env python

""" main function for energy model """

import numpy as np
import uproot as up

import global_param as glo
import param_config as pc
import plot_tpl as pl

import gamma_sources as gs
import simData as sd
import gamma_minimize as gmin


if __name__ == "__main__" :
    pc.param_config()  # parameter configuration firstly

    """ gamma prediction class """
    gamma_list = gs.gamma_sources()
    gamma_list.initialize()
    gamma_list.add_all_sources()
    gamma_list.read_all_sources_nonl()
    
    gamma_fitter = gmin.gamma_fitter(gamma_list)
    gamma_fitter.GetChiSquare()

    #pl.plot_func(gamma_list.get_etrue_list(), np.array(gamma_list.get_evis_list())/np.array(gamma_list.get_etrue_list()), "gamma nonlinearity", "o-", "blue")
    #pl.plot_err_func(sd.gamma_nonl_data()._fX, sd.gamma_nonl_data()._fY, "gamma nonlinearity data", sd.gamma_nonl_data()._fEX, sd.gamma_nonl_data()._fEY, "s", "hotpink")
    #pl.show(False, "Etrue/MeV", "Evis/Etrue")

    

