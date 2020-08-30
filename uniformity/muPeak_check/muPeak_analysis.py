#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Sun Aug 30 20:12:42 2020
# File Name: muPeak_analysis.py
"""

import uproot as up
from ROOT import TFile, TH1D
import sys

def read_root(name):
    hist_array = []
    try:
        ff = TFile(name, "read")
        hist =  TH1D(ff.Get("hComp_pos1") )
        print(hist)
        #for i in range(1):
        #    #hist_array[i] = ff.Get("hComp_pos%d"%i)
        #    hist = ff.Get("hComp_pos0")
        return hist
    except:
        print("No Such root file! ")
        sys.exit(0)




if __name__ == "__main__":

    hist_array = read_root("/junofs/users/huanggh/EnergyRec/GenCalibData/ACU_CLS_MAP/nPEMap_Truth/J20v1r1-Pre0/e+2MeVScaleCompFile.root")
    print(hist_array.GetNbinsX())
    
    """ find peak in whole range"""
    peak_value_arr, peak_pos_arr = [], []

    #for ipos in range(1):
    #    peak_value = -1
    #    peak_pos = -1
    #    for ibin in range(hist_array[ipos].GetNbinsX()):
    #        if hist_array[ipos].GetBinContent(ibin) > peak_value:
    #            peak_value = hist_array[ipos].GetBinContent(ibin+1)
    #            peak_pos = hist_array[ipos].GetBinCenter(ibin+1)
    #    if peak_value != 1:
    #        peak_value_arr.append(peak_value)
    #        peak_pos_arr.append(peak_pos)

    #plt.plot(peak_pos_arr, peak_value_arr, "o-")
    #plt.show()
