from Parameter_Sweep import *
import subprocess


#SmallerRange/RunNumber53/PhotochemInputs/', 
test_object = Generate_Atmosphere_Parameter_Sweep('Vdep4e-1T3', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/SmallerRange/RunNumber53/PhotochemInputs/', 
                                restart_run= 'Vdep4e-1T2', 
                                starting_point='Exact',
                                hitran_year='2020')

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
