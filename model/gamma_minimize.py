#coding=utf-8
#!/usr/bin/env python

""" gamma fitter (minimize) """

import global_param as glo
import gamma_sources as gs

from ctypes import *
import numpy as np
from ROOT import TMinuit

class gamma_fitter(object):

    """ Minimization for gamma chi2 fitter """

    def __init__(self, cls):
        self.gamma_cls = cls

    def scan_parameter(self):
        pass

    def get_chi2(self):
        print("Current chi2: %.4f" %self.gamma_cls.chi2_nonl(2) )
        return  self.gamma_cls.chi2_nonl(2)

    def ChisqFCN(self, npar, grad, fval, par, flag):
        self.SetParameters(par)
        fval = self.get_chi2()

    def SetParameters(self, par):
        print("Parameters Config: {:.4f}, {:.4f}, {:.4f}" .format(par[0], par[1], par[2]) )
        glo.set_value("kA", par[0])
        glo.set_value("kB", par[1])
        glo.set_value("kC", par[2])

    def GetChiSquare(self):
        gammaNLMinuit = TMinuit()
        gammaNLMinuit.SetFCN(self.ChisqFCN)
        gammaNLMinuit.SetPrintLevel(1)

        arglist = np.zeros((10, 1))
        ierrflag = 1
        iPar = 0
        gammaNLMinuit.mnexcm("CLEAR", arglist, 0, c_int(ierrflag))

        # Configurate Parameters
        gammaNLMinuit.mnparm(iPar, "kA", 0.97, 0.1, 0.93, 0.98, c_int(ierrflag));      iPar+=1
        gammaNLMinuit.mnparm(iPar, "kB", 65, 1, 50, 80, c_int(ierrflag));                iPar+=1;
        gammaNLMinuit.mnparm(iPar, "kC", 1.0, 0.01, 0.97, 1.03, c_int(ierrflag));       iPar+=1;

        gammaNLMinuit.FixParameter(1)
        gammaNLMinuit.FixParameter(2)

        # Minimization strategy
        gammaNLMinuit.SetErrorDef(1)
        arglist[0] = 2
        gammaNLMinuit.mnexcm("SET STR", arglist, 1, c_int(ierrflag))

        arglist[0] = 5000 #maxCalls
        arglist[1] = 0.01 # tolerance
        gammaNLMinuit.mnexcm("MIGRAD", arglist, 1, c_int(ierrflag))

        m_min, m_edm, m_errdef = 0, 0, 0
        nvpar, nparx, icstat = 0, 0, 0
        gammaNLMinuit.mnstat(c_double(m_min), c_double(m_edm), c_double(m_errdef), c_int(nvpar), c_int(nparx), c_int(icstat) )

        print ("=======================")
        print (" minChi2: %.3f" %m_min  ) 
        print ("=======================")
