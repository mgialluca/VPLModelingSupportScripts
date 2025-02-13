from Parameter_Sweep import *

test_object = Generate_Atmosphere_Parameter_Sweep('SweepTry1', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/PsurfTest/', 
                                restart_run= True, 
                                hitran_year='2020')


test_object.run_grid_sweep()
test_object.compile_run_output()