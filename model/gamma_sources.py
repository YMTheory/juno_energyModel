#coding=utf-8
#!/usr/bin/env python

""" prediciont mode module """
""" gamma prediction sources list  """

import global_param as glo
import gamma_prediction as gp
import uproot as up
import uproot_methods

class gamma_sources(object):

    """A list to contain all calibration sources data """

    calib_sources = []  # list for all calibration sources
    sources_name = {
            'Cs137' : 0.662,
            'Mn54'  : 0.834,
            'K40'   : 1.461,
            'nH'    :  2.222,
            'Tl208' : 2.614,
            'nC12'  : 4.945,
            'O16'   : 6.130,
            'nFe56' : 7.637
            }


    pred_nonl_data = {}
    sim_nonl_data = {}
    pred_resol_data = {}
    sim_resol_data = {}


    @staticmethod
    def all_sources():
        print(self.sources_name)
        return self.sources_name


    def initialize(self):
        self.calib_sources.clear()
        self.sim_nonl_data.clear()
        self.sim_resol_data.clear()
        self.pred_nonl_data.clear()
        self.pred_resol_data.clear()


    def add_source(self, name):
        try:
            source = gp.gamma_prediction(name, self.sources_name[name])
            calib_sources.append(source)
        except KeyError:
            print("No such calibration sources in list: {}" .format(name))



    def add_all_sources(self):
        for key in self.sources_name:
            source = gp.gamma_prediction(key, self.sources_name[key])
            self.calib_sources.append(source)


    def predict(self, option):
        if option == 1:
            for isource in self.calib_sources:
                evis, sigma = isource.predict_option1()
                key = isource.get_name() 
                self.pred_nonl_data[key]  = evis/isource.get_energy()
                self.pred_resol_data[key] = sigma/evis
        elif option == 2:
            for isource in self.calib_sources:
                nonl = isource.predict_option2()
                key = isource.get_name()
                self.pred_nonl_data[key] = nonl
        else:
            print(" Unknown prediction option !!! ")


    def read_all_sources_nonl(self):
        filename = glo.get_value("simCenterNonlFile")
        try:
            file = up.open(filename)
            gamma_graph = file["gamma"]
            sim_num = len(gamma_graph._fX)
            if sim_num != len(self.sources_name):
                print("SimData and PredData disagree on number")

            """Here we require that calibration sources are consistent in TGraphError data and dict in this file!!! Cs, Mn, K ..."""
            for name, nonl, nonlerr in zip(self.sources_name, gamma_graph._fY, gamma_graph._fEY):
                self.sim_nonl_data[name] = [nonl, nonlerr]
        except IOError:
            print ("Error: can not open {}" .format(filename)  )
    

    def chi2_nonl(self, option):
        chi2 = 0
        self.predict(option)
        # pre-check :
        for key in self.pred_nonl_data:
            if key not in self.sim_nonl_data:
                return chi2
        for key in self.pred_nonl_data:
            chi2 += (self.pred_nonl_data[key] - self.sim_nonl_data[key][0]) **2 / self.sim_nonl_data[key][1]
        return chi2

        
