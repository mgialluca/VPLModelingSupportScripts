from Parameter_Sweep import *

test_object = Generate_Atmosphere_Parameter_Sweep('SweepTry3', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/StartStateSweep1/', 
                                restart_run= 'SweepTry2', 
                                hitran_year='2020')


test_object.run_grid_sweep()
test_object.compile_run_output()