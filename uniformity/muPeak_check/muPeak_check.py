#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Sun Aug 30 00:48:22 2020
# File Name: muPeak_check.py
"""

import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import sys
from ROOT import TProfile, TCanvas, gPad, TGraph, TFile

plt.style.use("ggplot")

def readFile(name):
    try:
        ff = up.open(name)
        hittime = ff["evt"].array("hitTime").flatten()
        pmtid = ff["evt"].array("pmtID").flatten()
        return hittime, pmtid

    except:
        print("No such file!")
        sys.exit(0)


def pmt_pos():
    """read lpmt positions"""
    pmtx, pmty, pmtz = [], [], []
    with open("./PMTPos_Acrylic_with_chimney.csv") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            pmtx.append(float(data[1]))
            pmty.append(float(data[2]))
            pmtz.append(float(data[3]))
    return pmtx, pmty, pmtz



if __name__ == "__main__":

    pmtx, pmty, pmtz= pmt_pos()
    hittime, pmtid = readFile("../data/AmC_16580_0_2709.root")

    mu_pmt = [ 0 for i in range(17612) ]
    hittime_pmt = [ 0 for i in range(17612) ]
    for time, idx in zip(hittime[0::], pmtid[0::]):
        if idx >= 17612:
            continue  ## skip spmts
        mu_pmt[idx] += 1
        hittime_pmt[idx] += time

    source_x, source_y, source_z = 16580.14834, 0, 2709
    source_r = np.sqrt(source_x**2+source_y**2+source_z**2)
    pmt_theta = [ np.arccos((source_x*pmtx[i]+source_y*pmty[i]+source_z*pmtz[i])/np.sqrt(pmtx[i]**2+pmty[i]**2+pmtz[i]**2)/source_r)/np.pi*180. for i in range(17612)  ]
    

    prof1 = TProfile("prof1", "", 600, 0, 180, 0, 2)
    prof2 = TProfile("prof2", "", 600, 0, 180, 0, 2000)
    for i in range(17612):
        if mu_pmt[i] == 0 :
            continue
        prof1.Fill(pmt_theta[i], mu_pmt[i]/5000)
        prof2.Fill(pmt_theta[i], pmt_theta[i]/mu_pmt[i])
    graph1 = TGraph(); graph1.SetName("mu2theta")
    graph2 = TGraph(); graph2.SetName("time2theta")
    for i in range(prof1.GetNbinsX()):
        graph1.SetPoint(i, prof1.GetBinCenter(i+1), prof1.GetBinContent(i+1))
        graph2.SetPoint(i, prof2.GetBinCenter(i+1), prof2.GetBinContent(i+1))
    
    graph1.SetLineColor(34)
    graph2.SetLineColor(44)
    out = TFile("mu2theta.root", "recreate")
    graph1.Write()
    graph2.Write()
    out.Close()
    
