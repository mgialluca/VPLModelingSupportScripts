from Parameter_Sweep import *
import subprocess


test_object = Generate_Atmosphere_Parameter_Sweep('SweepTry6', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/StartStateSweep1/', 
                                restart_run= 'SweepTry5', 
                                starting_point='Exact',
                                hitran_year='2020')

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()