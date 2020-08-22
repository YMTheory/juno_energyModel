#coding=utf-8
#!/usr/bin/env python

""" Photon smearing from residual non-uniformity module """
""" use a constant model to describe currently """

import global_param as glo

def residual_nu_sigma2(NPE):
    """ residual non-uniformity induced extra sigma2 """
    p0 = glo.get_value("p0_RNU")
    p1 = glo.get_value("p1_RNU")
    p2 = glo.get_value("p2_RNU")
    val = p0 + p1*NPE + p2*NPE*NPE
    return val

