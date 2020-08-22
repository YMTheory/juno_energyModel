#coding=utf-8
#!/usr/bin/env python

""" LS photon smearing module """
""" parameterisation form: sigma2 = p0+p1*E+p2*E*E  """

import numpy as np
import uproot as up

import global_param as glo
import param_config as pc
import plot_tpl as pl
import prediction as pred
import simData as sd


if __name__ == "__main__" :
    pc.param_config()

    npe = [i for i in range(1000, 14000)]
    elec_resol_tot = [pred.electron_resolution_total(n, False, False) for n in npe]
    pl.plot_func([n/glo.get_value("energy_scale") for n in npe], elec_resol_tot, "LS NPE resolution")

    graph_elec = sd.electron_data( glo.get_value("simFile") )
    pl.plot_func(graph_elec._fX, graph_elec._fY, "electron data", "o")

    pl.show()

