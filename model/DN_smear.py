#coding=utf-8
#!/usr/bin/env python

""" Photon smearing from DN """
""" use a Poisson fluctuation to describe currently """

import global_param as glo

def DN_sigma2():
    DCRlpmt = glo.get_value("dcr_lpmt")
    length = glo.get_value("window_length")
    Nlpmt = glo.get_value("Nlpmt")
    val = Nlpmt * length * DCRlpmt   ## lpmt DN smear
    return val


