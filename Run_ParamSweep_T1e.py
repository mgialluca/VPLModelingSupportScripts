from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Einit', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1eSt/PhotochemInputs/', 
                                restart_run= False, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=False,
                                spectra=True,
                                planet='T1e')

# Outgassing rate for T1h
test_object.outgass_species_MinMax_gridsweep['H2O'] = [41227874.20930379, 5.1613525e11]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
