#coding=utf-8
#!/usr/bin/env python

""" simData module """
""" read simulation samples data """

import uproot as up
import uproot_methods

import global_param as glo

def electron_resol_data(filename):
    try:
        file = up.open(filename)
        #print("Start to read {}" .format(filename) )
        elec_graph = file["elec"]
        return elec_graph
    except IOError:
        print ("Error: can not open {}" .format(filename)  )



def electron_nonl_data():
    filename = glo.get_value("simCenterNonlFile")
    try:
        file = up.open(filename)
        #print("Start to read {}" .format(filename) )
        elec_graph = file["elec"]
        return elec_graph
    except IOError:
        print ("Error: can not open {}" .format(filename)  )


def gamma_nonl_data():
    filename = glo.get_value("simCenterNonlFile")
    try:
        file = up.open(filename)
        #print("Start to read {}" .format(filename) )
        gamma_graph = file["gamma"]
        return gamma_graph
    except IOError:
        print ("Error: can not open {}" .format(filename)  )


