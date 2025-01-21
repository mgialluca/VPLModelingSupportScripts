from Parameter_Sweep import *

test_object = Generate_Atmosphere_Parameter_Sweep('TestSweep', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/T1c/T1cOutgas_Testing/', 
                                  find_molecules_of_interest=False, hitran_year='2020')


test_object.run_grid_sweep()