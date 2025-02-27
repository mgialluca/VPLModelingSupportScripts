from Parameter_Sweep import *
import subprocess


#SmallerRange/RunNumber53/PhotochemInputs/', 
test_object = Generate_Atmosphere_Parameter_Sweep('TestVDEP2', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/InputDir_ForTesting/', 
                                restart_run= True, 
                                starting_point=None,
                                hitran_year='2020')

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
