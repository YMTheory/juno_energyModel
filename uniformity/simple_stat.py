import numpy as np
import random as rd
import matplotlib.pyplot as plt
from ROOT import TH1D
from rootpy.interactive import wait

h1 = TH1D("hist", "", 200, 700, 1300);
sigma = 50
for i in range(100000): # 200 bins to sample
    mu = rd.gauss(1000,50)
    sample = rd.gauss(mu, sigma)
    h1.Fill(sample)

print("%.2f, %.2f"%(h1.GetMean(), h1.GetStdDev())) 
h1.Draw()
wait()
