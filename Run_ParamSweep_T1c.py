from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Czoom1', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1cSt/PhotochemInputs/', 
                                restart_run= True, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=True,
                                spectra=True,
                                planet='T1c')

# Outgassing rate for T1c
test_object.outgass_species_MinMax_gridsweep['H2O'] = [34208000000.0, 9.47899801e11]#[44552887.2545331, 9.47899801e11]


test_object.outgass_sample_resolution_gridsweep = [16]

test_object.escape_samples_gridsweep['O'] = [0.01, 0.1]
test_object.escape_samples_gridsweep['O2'] = [0.001, 0.01]
test_object.escape_samples_gridsweep['O3'] = [0.02] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02]


test_object.run_grid_sweep()
test_object.compile_run_output()
