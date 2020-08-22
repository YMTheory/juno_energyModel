#coding=utf-8
#!/usr/bin/env python

""" prediciont mode module """
""" electron resolution model prediction  """

import math

import global_param as glo
import LS_smear as LS
import uniform_smear as uni
import DN_smear as dn

def electron_sigma2_LS(NPE):
    val = LS.LSPE_sigma2(NPE)
    return val

def electron_sigma2_RNU(NPE):
    val = uni.residual_nu_sigma2(NPE)
    return val

def electron_sigma2_DN(NPE):
    val = dn.DN_sigma2(NPE)
    return val


def electron_resolution_total(NPE, hasRSU=True, hasDN=True):
    val = LS.LSPE_sigma2(NPE)
    if hasRSU:
        val += uni.residual_nu_sigma2(NPE)
    if hasDN:
        val += dn.DN_sigma2(NPE)
    return math.sqrt(val)/NPE

