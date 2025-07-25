from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('BT2p2', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1bSt/PhotochemInputs/', 
                                restart_run= False, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=True,
                                spectra=True,
                                planet='T1b')

# Outgassing rate for T1b
test_object.outgass_species_MinMax_gridsweep['H2O'] = [36036000000.0, 1011500000000.0]
test_object.outgass_sample_resolution_gridsweep = [4]

test_object.escape_sample_type_gridsweep = ['Linear', 'UserDef', 'UserDef', 'UserDef']
test_object.escape_sample_resolution_gridsweep = [5, 0, 0, 0]

test_object.escape_samples_gridsweep['O'] = []
test_object.escape_species_MinMax_gridsweep['O'] = [0.175, 0.835]

test_object.escape_samples_gridsweep['O2'] = [0.1]#[1e26, 5e26, 1e27]
test_object.escape_samples_gridsweep['O3'] = [0.02] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
