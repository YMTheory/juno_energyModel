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


def electron_sigma2_total(NPE, hasRSU=False, hasDN=False):
    val = LS.LSPE_sigma2(NPE)
    if hasRSU:
        val += uni.residual_nu_sigma2(NPE)
    if hasDN:
        val += dn.DN_sigma2(NPE)
    return val


def electron_resolution_total(NPE, hasRSU=False, hasDN=False):
    val = LS.LSPE_sigma2(NPE)
    if hasRSU:
        val += uni.residual_nu_sigma2(NPE)
    if hasDN:
        val += dn.DN_sigma2(NPE)
    return math.sqrt(val)/NPE


def electron_sigma2_total_option2(Etrue):
    p0 = 2.57576e-2
    p1 = 6.86431e-3
    p2 = 0.
    return (p0*p0/Etrue+p1*p1+p2*p2/Etrue/Etrue)*Etrue**2
