from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Ezoom2', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1eSt/PhotochemInputs/', 
                                restart_run= True, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=False,
                                spectra=True,
                                planet='T1e')

# Outgassing rate for T1h
test_object.outgass_species_MinMax_gridsweep['H2O'] = [22228000000.0, 5.1613525e11]#[41227874.20930379, 5.1613525e11]

test_object.outgass_sample_resolution_gridsweep = [16]

test_object.escape_samples_gridsweep['O'] = [0.00001, 0.00005]
test_object.escape_samples_gridsweep['O2'] = [0.00001, 0.00005]
test_object.escape_samples_gridsweep['O3'] = [0.02] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
