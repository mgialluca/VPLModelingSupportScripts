from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Finit2', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1fInpR/PhotochemInputs/', 
                                restart_run= True, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=False,
                                spectra=True,
                                planet='T1f')

# Outgassing rate for T1h
#test_object.adjust_N2 = 0.05 # 0.05 bars of N2 fixed

#test_object.outgass_species_gridsweep = ['H2O', 'CO2']
test_object.outgass_species_MinMax_gridsweep['H2O'] = [48906843.61425449, 8.1495329e11]
#test_object.outgass_species_MinMax_gridsweep['CO2'] = [97544.61583408, 4.7368185e+10]
#test_object.outgass_species_molarmass['CO2'] = [44.01]*(u.g/u.mol)
#test_object.outgass_sample_type_gridsweep = ['Log', 'Log']

test_object.outgass_sample_resolution_gridsweep = [4]#, 4]

test_object.escape_samples_gridsweep['O'] = [0.01, 0.1]
test_object.escape_samples_gridsweep['O2'] = [0.01, 0.05]
test_object.escape_samples_gridsweep['O3'] = [0.02, 0.4] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02, 0.4]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
