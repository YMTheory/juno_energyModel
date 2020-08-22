import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile, TH1D

Nbins = 180
rsd_ratio = np.zeros((5, Nbins))
rsd_sigma = np.zeros((5, Nbins))

def read_corr(index):
    resol = ["0", "100", "230", "500"]
    global Nbins
    global rsd_ratio
    global rsd_sigma
    ff = TFile("./smear"+resol[index]+"mm.root", "read")
    for i in range(Nbins):
        h = ff.Get("sub%d" %i)
        rsd_ratio[index][i] = h.GetMean()
        rsd_sigma[index][i] = h.GetStdDev()

if __name__ == '__main__':
    resol = ["no smear", "smear100mm", "smear230mm", "smear500mm"]
    for idx in range(len(resol)):
        print("Processing " + str(idx))
        read_corr(idx);
        plt.plot([subid for subid in range(Nbins)], rsd_sigma[idx][:]/rsd_ratio[idx][:], "o-", ms=2, label=resol[idx])
    plt.legend(loc="upper left")
    plt.xlabel("sub-detector ID"); plt.ylabel("resolution")
    plt.title("Sub-Detector Resolution")
    plt.vlines(150, 0.02, 0.07, color="darkviolet")
    plt.text(100, 0.055, "total reflection \n region", color="darkviolet")
    #plt.show()
    plt.savefig("subDet-resol.pdf")
