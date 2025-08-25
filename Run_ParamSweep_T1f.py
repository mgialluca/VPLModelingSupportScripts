from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Finit', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1fSt/PhotochemInputs/', 
                                restart_run= True, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=False,
                                spectra=True,
                                planet='T1f')

# Outgassing rate for T1h
test_object.outgass_species_MinMax_gridsweep['H2O'] = [48906843.61425449, 8.1495329e11]
#test_object.adjust_N2 = 0.05 # 0.05 bars of N2 fixed

test_object.outgass_sample_resolution_gridsweep = [4]

test_object.escape_samples_gridsweep['O'] = [0.01, 0.1]
test_object.escape_samples_gridsweep['O2'] = [0.01, 0.05]
test_object.escape_samples_gridsweep['O3'] = [0.02, 0.4] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02, 0.4]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
