from Parameter_Sweep import *
import subprocess


test_object = Generate_Atmosphere_Parameter_Sweep('MultiNTest', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VeffTestDepos/RunNumber6/PhotochemInputs/', 
                                restart_run= 'VeffTestDepos', 
                                starting_point='Euclidean',
                                hitran_year='2020')

test_object.match_data_multinest()
