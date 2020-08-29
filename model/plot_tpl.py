#coding=utf-8
#!/usr/bin/env python

""" Plot template scripts """
import matplotlib.pyplot as plt

def plot_func(x, y, ll, myfmt="-", mycolor=None):
    if mycolor == None:
        plt.plot(x, y, myfmt, label=ll)
    else:
        plt.plot(x, y, myfmt, color=mycolor, label=ll)


def plot_err_func(x, y, ll, x_err=None, y_err=None, myfmt='-', mycolor=None):
    if mycolor == None:
        plt.errorbar(x, y, xerr=x_err, yerr=y_err, fmt=myfmt, label=ll)
    else:
        plt.errorbar(x, y, xerr=x_err, yerr=y_err, fmt=myfmt, color=mycolor, label=ll)



def show(xlog=False, xlabel=None, ylabel=None):
    if xlog:
        plt.semilogx()
    plt.legend(loc='best')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
