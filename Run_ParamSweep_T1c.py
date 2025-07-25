from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Cinit', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1cSt/PhotochemInputs/', 
                                restart_run= True, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=True,
                                spectra=True,
                                planet='T1c')

# Outgassing rate for T1c
test_object.outgass_species_MinMax_gridsweep['H2O'] = [44552887.2545331, 9.47899801e11]

'''
test_object.outgass_sample_resolution_gridsweep = [4]

test_object.escape_sample_type_gridsweep = ['UserDef', 'UserDef', 'UserDef', 'UserDef']
test_object.escape_sample_resolution_gridsweep = []

test_object.escape_samples_gridsweep['O'] = [0.01]
test_object.escape_samples_gridsweep['O2'] = [0.001, 0.005, 0.05, 0.15, 0.2]#[1e26, 5e26, 1e27]
test_object.escape_samples_gridsweep['O3'] = [0.02] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02]

#test_object.compile_info_failed_run()
'''

test_object.run_grid_sweep()
test_object.compile_run_output()
