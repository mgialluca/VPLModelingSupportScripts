from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Dzoom1', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1dSt/PhotochemInputs/', 
                                restart_run= True, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=True,
                                spectra=True,
                                planet='T1d')

# Outgassing rate for T1d
test_object.outgass_species_MinMax_gridsweep['H2O'] = [732290000.0, 2.97888812e11]#[36307296.51150426, 2.97888812e11]

test_object.outgass_sample_resolution_gridsweep = [16]

test_object.escape_samples_gridsweep['O'] = [0.0001, 0.001]
test_object.escape_samples_gridsweep['O2'] = [0.0001, 0.001]
test_object.escape_samples_gridsweep['O3'] = [0.02] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02]

test_object.run_grid_sweep()
test_object.compile_run_output()
