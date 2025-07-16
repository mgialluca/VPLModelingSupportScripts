from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Dinit', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1dSt/PhotochemInputs/', 
                                restart_run= False, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=True,
                                spectra=True,
                                planet='T1d')

# Outgassing rate for T1d
test_object.outgass_species_MinMax_gridsweep['H2O'] = [36307296.51150426, 2.97888812e11]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
