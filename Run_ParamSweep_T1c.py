from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Cinit', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1cSt/PhotochemInputs/', 
                                restart_run= False, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=True,
                                spectra=True,
                                planet='T1c')

# Outgassing rate for T1c
test_object.outgass_species_MinMax_gridsweep['H2O'] = [44552887.2545331, 9.47899801e+11]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
