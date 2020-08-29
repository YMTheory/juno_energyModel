#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Thu Aug 27 14:03:38 2020
# File Name: calib_points_stat.py
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ROOT import TH2D, TCanvas, TProfile2D, gStyle, TFile
import uproot as up

r3max = 17.2**3
deltar = r3max / 20.

#data = np.loadtxt("./OrgCalibPos.txt", delimiter=" ")
data = np.loadtxt("./OrgCalibPos_Random2000.txt", delimiter=" ")
posx, posz = data[:, 0], data[:, 1]

print("Total CLS Points: %d" % len(posx) )

edepR = ( np.sqrt(posx**2 + posz**2)/1000 ) 
edepR3 = ( np.sqrt(posx**2 + posz**2)/1000 ) **3
costh = (posz/1000/edepR)

mtx = np.zeros((10, 20))
for r, t in zip(edepR3, costh):
    rbin = int(r/deltar)
    #print("%.2f, %.2f, %d" %(r, deltar, rbin) )
    thetabin = int( (t+1)/(2/10.) )
    #print("%.2f,  %d" %(thetabin, rbin) )
    if thetabin >= 10:
        thetabin = 9
    if rbin >= 20:
        continue;
        #rbin = 18
    mtx[thetabin][rbin] += 1
    
print(mtx)

hist2d_source_number = TH2D("source_number", "", 20, 0, r3max, 10, -1 ,1)
for i in range(20):
    for j in range(10):
        hist2d_source_number.SetBinContent( i+1, j+1, mtx[j][i] )


#infile = up.open("../../../uniformity/gamma/Ge68.root")
infile = up.open("../data/Ge68_CLS2000_new.root")
edepR = infile["petree"].array("edepR")
evis = infile["petree"].array("evis")
edepZ = infile["petree"].array("edepZ")

prof2d = TProfile2D("prof2d", "", 20, 0, r3max, 10, -1, 1)
for r, e, z in zip(edepR, evis, edepZ):
    r3 = (r/1000.)**3
    cos = z/r
    prof2d.Fill(r3, cos, e)
    


gStyle.SetOptStat(0)
out = TFile("CLS2000_new.root", "recreate")
#cc = TCanvas()
#prof2d.SetTitle("2000 CLS; R^3/m^3; cos_theta")
#prof2d.SetBarOffset(0.2)
#prof2d.Draw("COLZ")
#hist2d_source_number.SetBarOffset(-0.2)
#hist2d_source_number.Draw("TEXT SAME")
prof2d.Write()
hist2d_source_number.Write()
out.Close()
#cc.SaveAs("CLS2000.pdf")


#plt.plot(edepR**3, costh, "o", ms=3)
#plt.xlabel(r"R^3")
#plt.ylabel("cos_theta")
#plt.show()

#f, ax = plt.subplots(figsize=(8, 4))
#cmap = sns.diverging_palette(220, 10, as_cmap=True)
#sns.heatmap(matrix, cmap=cmap, vmax=.3, center=0,
#            square=True, linewidths=.5, cbar_kws={"shrink": .5})
#sns.set(font_scale=0.7)
#sns.heatmap(mtx, cmap="viridis", annot=True,fmt=".0f",  square=True, linewidths=.5, cbar_kws={"shrink": .5})
#ax.set_xticklabels([str(int(r*deltar)) for r in range(20)])
#ax.set_yticklabels(["%.2f"%((t*0.2-1)) for t in range(10)])
#ax.set_yticklabels(["200", "210", "220", "230", "240", "250", "260", "270", "280", "290", "300"])
#ax.set_yticklabels(["100", "130", "160", "200", "220", "250", "280", "310", "340", "370", "400"])
#ax.set_xlabel(r"$r^3/m^3$")
#ax.set_ylabel(r"$cos \theta$")
#ax.set_title("240 CLS", fontsize= 10);

#plt.savefig("slip1.pdf")
#plt.show()
