from Parameter_Sweep import *
import subprocess


test_object = Generate_Atmosphere_Parameter_Sweep('MCMCTest', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/PhotochemInputs/', 
                                restart_run= True, 
                                starting_point=None,
                                hitran_year='2020')

test_object.match_surf_pressure_MCMC()
