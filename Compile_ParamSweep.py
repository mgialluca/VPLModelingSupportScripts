from Parameter_Sweep import *
import subprocess


test_object = Generate_Atmosphere_Parameter_Sweep('Sens5TestH2O', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/SmallerRange/RunNumber53/PhotochemInputs/', 
                                restart_run= False, 
                                hitran_year='2020')

test_object.compile_info_failed_run(Num_of_Models=192)

#test_object.compile_smart_spectra(Num_of_Models=150)

#test_object.compile_restart_input_options()
#test_object.run_grid_sweep()
#test_object.compile_run_output()
