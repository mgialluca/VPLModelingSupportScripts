from Pipeline import *
import copy
import subprocess
import os
from multiprocessing import Pool

master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'

def set_pipeline_vars(casename, pipelineobj, master_out=master):

    # Paths are the main thing to set, because they will be massive amounts of running/files, want to keep each sweep colocated in one master dir
    atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'
    # Create casename dir
    if not os.path.exists(master_out+casename+'/'):
        os.mkdir(master_out+casename+'/')
    
    if not os.path.exists(master_out+casename+'/atmos/'):
        shutil.copytree(atmos_Dir,  master_out + casename + '/atmos/')

    # Remove PhotGrid.f and replace with GJ12 updated one; put spectrum in DATA/FLUX
    os.remove(master_out+casename+'/atmos/PHOTOCHEM/SUBROUTINES/Photgrid.f')
    shutil.copy('/gscratch/vsm/alinc/VPL_RUNS/gj12b/Photgrid.f', master_out+casename+'/atmos/PHOTOCHEM/SUBROUTINES/Photgrid.f')
    shutil.copy('/gscratch/vsm/alinc/VPL_RUNS/gj12b/steam/1bar/atmos/PHOTOCHEM/DATA/FLUX/gj12.txt', master_out+casename+'/atmos/PHOTOCHEM/DATA/FLUX/gj12.txt')

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
    pipelineobj.include_2column_climate = False
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = True
    pipelineobj.dayside_starting_PT = None
    pipelineobj.nightside_starting_PT = None
    pipelineobj.NewPressure_Psurf_tolerance = 0.035

    pipelineobj.c_NumberSolarZeniths = 4

    pipelineobj.MCMC_pressure_only = False
    pipelineobj.MultiNest_DataFit = False
    pipelineobj.rerun_smart_for_2col = False

    pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate_supernode'

    # Molecules for the type of atmosphere we're interested in 

    pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
    gas_names = ['O2', 'H2O', 'O3', 'H2O2']#, 'CO2', 'CO', 'CH4', 'N2O']
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

    if case == 'GJ12b01':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/GJ12b_Starts/O2_01bar/'
        planet = 'GJ12b'

    elif case == 'GJ12b1':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/GJ12b_Starts/O2_1bar/'
        planet = 'GJ12b'

    elif case == 'GJ12b10':
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/GJ12b_Starts/O2_10bar/'
        planet = 'GJ12b'

    pipelineobj = VPLModelingPipeline(case, 
                                    initin, 
                                    True, find_molecules_of_interest=False, hitran_year='2020', planet=planet)
    
    set_pipeline_vars(case, pipelineobj)

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




inputs = ['GJ12b01', 'GJ12b1', 'GJ12b10']#['b01', 'b1', 'b10']

with Pool() as p:
    models = p.map(run_starting_points, inputs)


