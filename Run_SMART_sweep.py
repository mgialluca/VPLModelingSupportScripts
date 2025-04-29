
import copy
from Pipeline import *
import subprocess
import shutil
import os
from multiprocessing import Pool


master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/SO2-H2OT1c/'

def set_pipeline_vars(casename, pipelineobj, master_out=master):

    # Paths are the main thing to set, because they will be massive amounts of running/files, want to keep each sweep colocated in one master dir
    atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'
    # Create casename dir
    if not os.path.exists(master_out+casename+'/'):
        os.mkdir(master_out+casename+'/')
    
    #shutil.copytree(atmos_Dir,  master_out + casename + '/atmos/')

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
    pipelineobj.photochem_InputsDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/T100mbar/' # The path to create new photochem inputs in
    pipelineobj.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found
    pipelineobj.SMART_RunScriptDir = master_out+casename+'/'

    # Adjust the atmospheric pressure
    pipelineobj.adjust_atmospheric_pressure = False
    pipelineobj.global_convergence = True # BYPASS PHOTOCHEM
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = False
    pipelineobj.dayside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cComparison/PT_profile_dayside_'+casename+'.pt'
    pipelineobj.nightside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cComparison/PT_profile_nightside_'+casename+'.pt'

    pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate_supernode'

    # Molecules for the type of atmosphere we're interested in 

    pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
    gas_names = ['O2', 'H2O', 'O3', 'SO2']# 'CO2', 'CO']
    pipelineobj.molecule_dict['Gas_names'] = gas_names
    for m in range(len(gas_names)):
        pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
        pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2

def rerun_with_so2(runname):

    os.mkdir(master+runname)

    copyfrom = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cComparison/'

    # Copy files over needed for LBLABC and SMART run
    fis_to_copy = ['PT_profile_'+runname+'.pt', 'MixingRs_'+runname+'.dat']
    
    for fi in fis_to_copy:
        shutil.copyfile(copyfrom+runname+'/'+fi, master+runname+'/'+fi)

    pipelineobj = VPLModelingPipeline(runname, 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/T100mbar/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')
    
    set_pipeline_vars(runname, pipelineobj)
    
    pipelineobj.run_automatic()

    pipelineobj.make_lblabc_runscripts(whichcol='dayside')

    # Now Run LBLABC for all the gases of interest
    for gas in pipelineobj.molecule_dict['Gas_names']:
        pipelineobj.run_lblabc_1instance(pipelineobj.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+pipelineobj.casename+'.script', gas)
    pipelineobj.num_lblabc_runs += 1

    pipelineobj.make_smart_runscript(whichcol='dayside')

    pipelineobj.run_smart_1instance(pipelineobj.SMART_RunScriptDir+'RunSMART_dayside_'+pipelineobj.casename+'.run', whichcol='dayside')


names = ['Run'+str(i) for i in range(80)]

with Pool() as p:
    models = p.map(rerun_with_so2, names)