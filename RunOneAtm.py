from Pipeline import *
import copy
import subprocess
import os

master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/SZADependenceTest/'

def set_pipeline_vars(casename, pipelineobj, master_out=master):

    # Paths are the main thing to set, because they will be massive amounts of running/files, want to keep each sweep colocated in one master dir
    atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'
    # Create casename dir
    #if not os.path.exists(master_out+casename+'/'):
    #    os.mkdir(master_out+casename+'/')
    
    shutil.copytree(atmos_Dir,  master_out + casename + '/atmos/')

    pipelineobj.photochemDir = master_out+casename+'/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
    pipelineobj.atmosDir = master_out+casename+'/atmos/' # path to atmos/ dir
    pipelineobj.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
    pipelineobj.OutPath = master_out+casename+'/' # path for the raw model run outputs (NOT for created data products like dictionaries)
    pipelineobj.DataOutPath = master_out+casename+'/' # path for created data products like dictionaries
    pipelineobj.AtmProfPath = master_out+casename+'/' # path to put atmospheric profile files (.pt files really)
    pipelineobj.BackupPhotochemRuns = False # Make backups of individual photochem runs
    pipelineobj.photochemBackupDir = master_out+casename+'/PhotochemBackup/' # path to save output from each photochem run
    pipelineobj.LBLABC_AbsFilesDir = master_out+casename+'/ABSFiles/' # path to put the created lbl .abs files in 
    pipelineobj.lblabc_RunScriptDir = master_out+casename+'/' # path to put lbl runscripts in
    pipelineobj.vplclimate_RunScriptDir = master_out+casename+'/' # path to put vpl climate runscripts in
    pipelineobj.photochem_InputsDir = master_out+casename+'/PhotochemInputs/' # The path to create new photochem inputs in
    pipelineobj.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found
    pipelineobj.SMART_RunScriptDir = master_out+casename+'/'

    # Adjust the atmospheric pressure
    pipelineobj.adjust_atmospheric_pressure = True
    pipelineobj.include_2column_climate = True
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = True
    pipelineobj.dayside_starting_PT = None
    pipelineobj.nightside_starting_PT = None
    pipelineobj.NewPressure_Psurf_tolerance = 0.035

    pipelineobj.c_NumberSolarZeniths = 1

    pipelineobj.adjust_atmospheric_pressure = True
    pipelineobj.suppress_IOerrors = True
    pipelineobj.MCMC_pressure_only = False
    pipelineobj.MultiNest_DataFit = False
    pipelineobj.rerun_smart_for_2col = True
    
    if pipelineobj.MCMC_pressure_only == True:
        pipelineobj.include_2column_climate = False
        pipelineobj.run_spectra = False
    else:
        pipelineobj.include_2column_climate = True
        pipelineobj.run_spectra = True

    if pipelineobj.MultiNest_DataFit == True:
        pipelineobj.multinest_climate_copycase = 'ClimTestMulti/Run99'
        copycase = pipelineobj.multinest_climate_copycase.split('/')[1]
        pipelineobj.dayside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/PT_profile_dayside_'+copycase+'.pt'
        pipelineobj.nightside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/PT_profile_nightside_'+copycase+'.pt'
        pipelineobj.run_spectra = True

        climoutcopy = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/vpl_2col_climate_output_'+copycase+'.run'
        fi = open(climoutcopy, 'r')
        lines = fi.readlines()
        fi.close()

        # Want to get the last output trop heating rate and avg flux, should be last two lines
        # so loop in reversed order, break loop after to conserve efficiency

        nightside_found = False
        for i in reversed(range(len(lines))):
            hold = lines[i].split()
            if len(hold) > 2:
                if hold[0] == 'surface:':
                    if nightside_found == False:
                        pipelineobj.surface_temp_nightside = float(hold[8])
                        nightside_found = True
                    else:
                        pipelineobj.surface_temp_dayside = float(hold[8])
                        # After retrieving surface temp for nightside, will have all values, break
                        break

    pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate_supernode'

    # Molecules for the type of atmosphere we're interested in 

    pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
    gas_names = ['O2', 'H2O', 'O3']
    pipelineobj.molecule_dict['Gas_names'] = gas_names
    for m in range(len(gas_names)):
        pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
        pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2

case = 'TestRun1'

pipelineobj = VPLModelingPipeline(case, 
                                  master+case+'/PhotochemInputs/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')
    
set_pipeline_vars(case, pipelineobj)

converged = pipelineobj.run_automatic()