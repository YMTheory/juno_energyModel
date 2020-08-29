# -*- coding: utf-8 -*-

import sys

def _init():  # Initialzation
    global _global_dict
    _global_dict = {}

def set_value(key,value):
    """ Define a global variable """
    _global_dict[key] = value

def get_value(key, defValue=None):
    #""" get a global variable,  default value if does not exist """
    try:
        return _global_dict[key]
    except KeyError:
        print("Do not have this key : {}, try to update parameters!" .format(key))
        sys.exit(1)
        #return defValue

def set_dict_value(key, value, val) :
    """ Set the value of a key as a dictionary """
    _global_dict.setdefault(key, {})[value] = val

def get_dict_value(key1, key2, defValue=None):
    try:
        return _global_dict[key1][key2]
    except KeyError:
        print("Do not have this key : {} {}, try to update parameters!" .format(key1, key2))
        sys.exit(1)
        #return defValue
        
