#coding=utf-8
#!/usr/bin/env python

""" LS photon smearing module """
""" parameterisation form: sigma2 = p0+p1*E+p2*E*E"""

import global_param as glo

def LSPE_sigma2(NPE):
    """ LS NPE sigma2 calculator"""
    p0 = glo.get_value("p0_LSPE")
    p1 = glo.get_value("p1_LSPE")
    p2 = glo.get_value("p2_LSPE")
    val = p0 + p1*NPE + p2*NPE*NPE
    return val
