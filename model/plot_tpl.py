#coding=utf-8
#!/usr/bin/env python

""" Plot template scripts """
import matplotlib.pyplot as plt

def plot_func(x, y, ll):
    plt.plot(x, y, label=ll)


def show():
    plt.legend(loc='best')
    plt.show()
