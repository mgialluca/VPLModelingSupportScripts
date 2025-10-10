from Pipeline import *
import copy
import subprocess
import os
from multiprocessing import Pool

master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'

def set_pipeline_vars(casename, pipelineobj, gas_names, master_out=master):

    # Paths are the main thing to set, because they will be massive amounts of running/files, want to keep each sweep colocated in one master dir
    atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'
    # Create casename dir
    if not os.path.exists(master_out+casename+'/'):
        os.mkdir(master_out+casename+'/')
    
    if not os.path.exists(master_out+casename+'/atmos/'):
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
    pipelineobj.adjust_atmospheric_pressure = False
    pipelineobj.include_2column_climate = True
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = True
    pipelineobj.dayside_starting_PT = None
    pipelineobj.nightside_starting_PT = None
    pipelineobj.NewPressure_Psurf_tolerance = 0.035

    pipelineobj.c_NumberSolarZeniths = 4

    pipelineobj.MCMC_pressure_only = False
    pipelineobj.MultiNest_DataFit = False
    pipelineobj.rerun_smart_for_2col = True

    pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate/vpl_climate'

    # Molecules for the type of atmosphere we're interested in 

    pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
    #gas_names = ['O2', 'H2O', 'O3', 'H2O2']#, 'CO2', 'CO', 'CH4', 'N2O']
    pipelineobj.molecule_dict['Gas_names'] = gas_names
    for m in range(len(gas_names)):
        pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
        pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2


'''
# To run one atmosphere:
case = 'Baseline'

pipelineobj = VPLModelingPipeline(case, 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/TEMPLATES/ModernEarth/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020', planet='Earth')
    
set_pipeline_vars(case, pipelineobj)

converged = pipelineobj.run_automatic()
'''


def run_starting_points(case):

    
    master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'

    if case == 'T1cH2O10ppm':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cCycle5/H2OOutg_10ppm/'
        planet = 'T1c'
        gas_names = ['O2', 'H2O', 'O3', 'H2', 'SO2']

    elif case == 'T1cH2O100ppm':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cCycle5/H2OOutg_100ppm/'
        planet = 'T1c'
        gas_names = ['O2', 'H2O', 'O3', 'H2', 'SO2']

    elif case == 'T1cCO210ppm':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cCycle5/CO2Outg_10ppm/'
        planet = 'T1c'
        gas_names = ['O2', 'H2O', 'O3', 'CO2', 'CO', 'SO2']

    elif case == 'T1cCO2100ppm':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cCycle5/CO2Outg_100ppm/'
        planet = 'T1c'
        gas_names = ['O2', 'H2O', 'O3', 'CO2', 'CO', 'SO2']

    elif case == 'T1cCO10ppm':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cCycle5/COdom_10ppm/'
        planet = 'T1c'
        gas_names = ['O2', 'H2O', 'O3', 'CO2', 'CO', 'SO2']

    elif case == 'T1cCO100ppm':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cCycle5/COdom_100ppm/'
        planet = 'T1c'
        gas_names = ['O2', 'H2O', 'O3', 'CO2', 'CO', 'SO2']

    pipelineobj = VPLModelingPipeline(case, 
                                    initin, 
                                    True, find_molecules_of_interest=False, hitran_year='2020', planet=planet)
    
    set_pipeline_vars(case, pipelineobj, gas_names)

    '''
    for sdshol, dshol, fishol in os.walk(master+case+'/'):
        fishol = fishol
        break

    for fhold in fishol:
        subprocess.run('rm '+master+case+'/'+fhold, shell=True)

    subprocess.run('rm -rf '+pipelineobj.LBLABC_AbsFilesDir+'*.abs', shell=True)
    '''

    converged = pipelineobj.run_automatic()

    return pipelineobj




inputs = ['T1cH2O10ppm', 'T1cH2O100ppm', 'T1cCO210ppm', 'T1cCO2100ppm', 'T1cCO10ppm', 'T1cCO100ppm']#['b01', 'b1', 'b10']

with Pool() as p:
    models = p.map(run_starting_points, inputs)


