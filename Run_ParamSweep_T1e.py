from Parameter_Sweep import *
import subprocess


#'/gscratch/vsm/gialluca/VPLModelingTools_Dev/SensTestH2O/RunNumber175/
test_object = Generate_Atmosphere_Parameter_Sweep('Eco2wL2', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1eco2/', 
                                restart_run= True, 
                                starting_point='Exact',
                                hitran_year='2020',
                                climate2col=False,
                                spectra=True,
                                planet='T1e')

# Outgassing rate for T1e
test_object.outgass_species_gridsweep = ['H2O', 'CO2']
test_object.outgass_species_MinMax_gridsweep['H2O'] = [58617668.80239138, 7.33839562e11] #[31604000000.0, 7.33839562e11]
test_object.outgass_species_MinMax_gridsweep['CO2'] = [118488.3029984, 4.22198983e+10]
test_object.outgass_species_molarmass['CO2'] = [44.01]*(u.g/u.mol)
test_object.outgass_sample_type_gridsweep = ['Log', 'Log']

test_object.outgass_sample_resolution_gridsweep = [4, 4]

test_object.escape_species_gridsweep = ['O', 'O2', 'O3', 'H2O2', 'CO', 'CO2'] # Species to vary escape rates of
test_object.escape_species_losstype = ['Veff', 'Veff', 'Vdep', 'Vdep', 'Vdep', 'Veff']
test_object.escape_species_molarmass['CO'] = [28.01]*(u.g/u.mol) 
test_object.escape_species_molarmass['CO2'] = [44.01]*(u.g/u.mol)
test_object.escape_sample_type_gridsweep = ['UserDef', 'UserDef', 'UserDef', 'UserDef', 'UserDef', 'UserDef'] 
test_object.escape_species_units_gridsweep = [u.cm/u.s, u.cm/u.s, u.cm/u.s, u.cm/u.s, u.cm/u.s, u.cm/u.s]

test_object.escape_samples_gridsweep['O'] = [0.01]
test_object.escape_samples_gridsweep['O2'] = [0.01]
test_object.escape_samples_gridsweep['O3'] = [0.02] 
test_object.escape_samples_gridsweep['H2O2'] = [0.02]
test_object.escape_samples_gridsweep['CO'] = [0.001, 0.0001]
test_object.escape_samples_gridsweep['CO2'] = [0.01, 0.1]

#test_object.compile_info_failed_run()
test_object.run_grid_sweep()
test_object.compile_run_output()
