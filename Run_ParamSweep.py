from Parameter_Sweep import *
import subprocess


#SmallerRange/RunNumber53/PhotochemInputs/', 
test_object = Generate_Atmosphere_Parameter_Sweep('SensTestH2O', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensitTest/RunNumber61/PhotochemInputs', 
                                restart_run= True, 
                                starting_point=None,
                                hitran_year='2020')

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
