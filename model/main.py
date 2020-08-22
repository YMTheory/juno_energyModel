#coding=utf-8
#!/usr/bin/env python

""" LS photon smearing module """
""" parameterisation form: sigma2 = p0+p1*E+p2*E*E  """

import numpy as np

import global_param as glo
import plot_tpl as pl
import LS_smear as LS
import uniform_smear as uni
import DN_smear as dn

def param_config():
    """ For current model, we must set up parameters configuration firstly! """
    glo._init()

    print("   ====> Configure All the Parameters in Current Model <====  ")
   
    # LS NPE sigma2
    glo.set_value("p0_LSPE",-3.69505e2 )
    glo.set_value("p1_LSPE", 1.30168)
    glo.set_value("p2_LSPE", 1.90109e-05)
    print(" In Global Configuration LSPE params : {: .2f}, {: .5f} and {: .5f}".format(glo.get_value("p0_LSPE"), glo.get_value("p1_LSPE"), glo.get_value("p2_LSPE")))

    # Residual Non-Uniformity sigma2
    glo.set_value("p0_RNU", 0)
    glo.set_value("p1_RNU", 0)
    glo.set_value("p2_RNU", 2.3e-5)
    print(" In Global Configuration RNU params : {: .1f}, {: .1f} and {: .5f}".format(glo.get_value("p0_RNU"), glo.get_value("p1_RNU"), glo.get_value("p2_RNU")))

    # DN sigma2
    glo.set_value("Nlpmt", 17612) # detector lpmt number config
    glo.set_value("window_length", 300e-9)  # unit s
    glo.set_value("dcr_lpmt", 27000) # unit: Hz
    print(" In Global Configuration NLpmt: {: d}".format(glo.get_value("Nlpmt")) )
    print(" In Global Configuration Window Length: {: .1f} s" .format(glo.get_value("window_length")))
    print(" In Global Configuration Lpmt DCR: {: .1f} Hz" .format(glo.get_value("dcr_lpmt")))

    # other config
    glo.set_value("energy_scale", 3350/2.22)  # energy scale
    print(" In Global Configuration Energy Scale: {: .2f} pe/MeV" .format(glo.get_value("energy_scale")))


if __name__ == "__main__" :
    param_config()

    npe = [i for i in range(10000)]
    npe_ls_sigma2 = [LS.LSPE_sigma2(n) for n in npe ]
    npe_uni_sigma2 = [uni.residual_nu_sigma2(n) for n in npe]
    npe_dc_sigma2 = [dn.DN_sigma2() for n in npe]
    npe_total_resol = [np.sqrt(a+b+c)/n for a,b,c,n in zip(npe_ls_sigma2, npe_uni_sigma2, npe_dc_sigma2, npe) ]
    pl.plot_func([n/glo.get_value("energy_scale") for n in npe], npe_total_resol, "total e- resolution")
    #pl.plot_func(npe, npe_ls_sigma2, "LS")
    #pl.plot_func(npe, npe_uni_sigma2, "RNU")
    #pl.plot_func(npe, npe_dc_sigma2, "DN")
    pl.show()

