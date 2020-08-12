# Documents for useful scripts 

nPEMap_corr1.cc  --> use to read samples with certain energy (uniform or center), divide into 40 shells, correct 3DnPE-Map and generate root tree in a file. processing in batches, need to specify start file index, and output file (arg1, arg2). Need to merge (hadd) all the outputs.

batch_process.cc --> read output files from nPEMap_corr1.cc and do gaussian fitting for evis in each shell. Output: evis v.s. radius, resol v.s. radius, evis v.s. radius(within 15m).

uniform_term.cc --> Similarly, do fitting for data and also do MC sampling to get whole CD resolution(w/ FV cut) and compare with data.

residual_uniform_rsl.cc --> compare resol v.s. evis for different cases.