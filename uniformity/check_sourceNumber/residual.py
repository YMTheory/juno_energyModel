#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Sat Aug 29 19:02:32 2020
# File Name: residual.py
"""

import uproot as up
import matplotlib.pyplot as plt
import numpy as np
import sys

def readFiles(name):
    """ read evis after nPEMap corrections """
    try:
        ff = up.open(name)
        evis = ff["petree"].array("evis_lpmt")
        edepR = ff["petree"].array("edepR")
        edepZ = ff["petree"].array("edepZ")
        return edepZ, edepR, evis
    except:
        print("No such file: {}" .format(name) )
        sys.exit(0)

plt.style.use("bmh")


if __name__ == "__main__" :

    edepZ1, edepR1, evis1 = readFiles("../data/Ge68_CLS2000_new.root")
    edepZ2, edepR2, evis2 = readFiles("../data/Ge68_CLS240_new.root")
    
    deltaR = 17200**3/40.

    evis1_array = [ [] for i in range(40) ] # 40 radius bins
    for i, j in zip(edepR1, evis1):
        index = int(i*i*i/deltaR)
        if index >= 40:
            continue
        evis1_array[index].append(j)
    evis2r1 = [ np.array(evis1_array[m]).mean() for m in range(40) ]

    evis2_array = [ [] for i in range(40) ] # 40 radius bins
    for i, j in zip(edepR2, evis2):
        index = int(i*i*i/deltaR)
        if index >= 40:
            continue
        evis2_array[index].append(j)
    evis2r2 = [ np.array(evis2_array[m]).mean() for m in range(40) ]

    plt.plot([r*deltaR for r in range(40) ], evis2r1, "o-", label="2000 CLS" )
    plt.plot([r*deltaR for r in range(40) ], evis2r2, "o-", label="240 CLS" )
    plt.legend()
    plt.xlabel(r"R^{3/m^{3}}")
    plt.ylabel("evis")
    plt.show()

    

