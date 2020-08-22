# -*- coding: utf-8 -*-

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
        return defValue
