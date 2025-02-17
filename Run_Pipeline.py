
from Pipeline import *
from multiprocessing import Pool

casename = 'climatebug'
master_out = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/testclimatebug/'
pipelineobj = VPLModelingPipeline('climatebug', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/SweepTry4/RunNumber71/PhotochemInputs/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

if not os.path.exists(master_out):
    os.mkdir(master_out)

if not os.path.exists(master_out+casename+'/'):
    os.mkdir(master_out+casename+'/')

shutil.copytree(pipelineobj.atmosDir, master_out + casename + '/atmos/')

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
pipelineobj.suppress_IOerrors = True
pipelineobj.run_spectra = True

# Testing if climate executable needs to be copied
shutil.copy(pipelineobj.vplclimate_executable, pipelineobj.OutPath+'vplclimate')
pipelineobj.vplclimate_executable = pipelineobj.OutPath+'vplclimate'

# Molecules for the type of atmosphere we're interested in 

pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
gas_names = ['O2', 'H2O', 'O3']
pipelineobj.molecule_dict['Gas_names'] = gas_names
for m in range(len(gas_names)):
    pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
    pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2


casename2 = 'climatebug2'
pipelineobj2 = VPLModelingPipeline('climatebug2', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/SweepTry4/RunNumber76/PhotochemInputs/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

#####


if not os.path.exists(master_out):
    os.mkdir(master_out)

if not os.path.exists(master_out+casename2+'/'):
    os.mkdir(master_out+casename2+'/')

shutil.copytree(pipelineobj2.atmosDir, master_out + casename2 + '/atmos/')

pipelineobj2.photochemDir = master_out+casename2+'/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
pipelineobj2.atmosDir = master_out+casename2+'/atmos/' # path to atmos/ dir
pipelineobj2.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
pipelineobj2.OutPath = master_out+casename2+'/' # path for the raw model run outputs (NOT for created data products like dictionaries)
pipelineobj2.DataOutPath = master_out+casename2+'/' # path for created data products like dictionaries
pipelineobj2.AtmProfPath = master_out+casename2+'/' # path to put atmospheric profile files (.pt files really)
pipelineobj2.BackupPhotochemRuns = False # Make backups of individual photochem runs
pipelineobj2.photochemBackupDir = master_out+casename2+'/PhotochemBackup/' # path to save output from each photochem run
pipelineobj2.LBLABC_AbsFilesDir = master_out+casename2+'/ABSFiles/' # path to put the created lbl .abs files in 
pipelineobj2.lblabc_RunScriptDir = master_out+casename2+'/' # path to put lbl runscripts in
pipelineobj2.vplclimate_RunScriptDir = master_out+casename2+'/' # path to put vpl climate runscripts in
pipelineobj2.photochem_InputsDir = master_out+casename2+'/PhotochemInputs/' # The path to create new photochem inputs in
pipelineobj2.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found
pipelineobj2.SMART_RunScriptDir = master_out+casename2+'/'

# Adjust the atmospheric pressure
pipelineobj2.adjust_atmospheric_pressure = True
pipelineobj2.suppress_IOerrors = True
pipelineobj2.run_spectra = True

# Testing if climate executable needs to be copied
shutil.copy(pipelineobj2.vplclimate_executable, pipelineobj2.OutPath+'vplclimate')
pipelineobj2.vplclimate_executable = pipelineobj2.OutPath+'vplclimate'

# Molecules for the type of atmosphere we're interested in 

pipelineobj2.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
gas_names = ['O2', 'H2O', 'O3']
pipelineobj2.molecule_dict['Gas_names'] = gas_names
for m in range(len(gas_names)):
    pipelineobj2.molecule_dict[gas_names[m]] = pipelineobj2.hitran_lookup.loc[gas_names[m]]['HitranNumber']
    pipelineobj2.molecule_dict[gas_names[m]+'_RmixCol'] = m+2

def run_it(obj):
    obj.run_automatic()

with Pool() as p:
    models = p.map(run_it, [pipelineobj, pipelineobj2])