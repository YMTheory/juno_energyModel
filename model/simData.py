#coding=utf-8
#!/usr/bin/env python

""" simData module """
""" read simulation samples data """

import uproot as up

def electron_data(filename):
    file = up.open(filename)
    elec_graph = file["elec"]
    return elec_graph
