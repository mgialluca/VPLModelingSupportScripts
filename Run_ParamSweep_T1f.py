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
test_object.outgass_species_MinMax_gridsweep['H2O'] = [44380169.91551189, 7.39523609e11]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
