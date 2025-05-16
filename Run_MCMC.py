from Parameter_Sweep import *
import subprocess


test_object = Generate_Atmosphere_Parameter_Sweep('MultiNTest', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VeffTestDepos/RunNumber166/PhotochemInputs/', 
                                restart_run= False, 
                                starting_point=None,
                                hitran_year='2020')

test_object.match_data_multinest()
