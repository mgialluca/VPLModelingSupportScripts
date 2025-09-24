from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Bco2wL', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1bco2/', 
                                restart_run= True, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=True,
                                spectra=True,
                                planet='T1b')

# Outgassing rate for T1b
test_object.outgass_species_gridsweep = ['H2O', 'CO2']
test_object.outgass_species_MinMax_gridsweep['H2O'] = [44193675.62502126, 9.77334769e11]#[34819000000.0, 9.77334769e11]
test_object.outgass_species_MinMax_gridsweep['CO2'] = [85186.07150893, 5.83683757e10]
test_object.outgass_species_molarmass['CO2'] = [44.01]*(u.g/u.mol)
test_object.outgass_sample_type_gridsweep = ['Log', 'Log']

test_object.outgass_sample_resolution_gridsweep = [4, 4]

test_object.escape_species_gridsweep = ['O', 'O2', 'O3', 'H2O2', 'CO', 'CO2'] # Species to vary escape rates of
test_object.escape_species_losstype = ['Veff', 'Veff', 'Vdep', 'Vdep', 'Vdep', 'Veff']
test_object.escape_species_molarmass['CO'] = [28.01]*(u.g/u.mol) 
test_object.escape_species_molarmass['CO2'] = [44.01]*(u.g/u.mol)
test_object.escape_sample_type_gridsweep = ['UserDef', 'UserDef', 'UserDef', 'UserDef', 'UserDef', 'UserDef'] 

test_object.escape_samples_gridsweep['O'] = [0.01] #[0.01, 0.1] #[0.01, 1]
test_object.escape_samples_gridsweep['O2'] = [0.01] #[0.01, 0.1]
test_object.escape_samples_gridsweep['O3'] = [0.02] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02]
test_object.escape_samples_gridsweep['CO'] = [0.0, 0.01]
test_object.escape_samples_gridsweep['CO2'] = [0.01, 0.1]

test_object.run_grid_sweep()
test_object.compile_run_output()

