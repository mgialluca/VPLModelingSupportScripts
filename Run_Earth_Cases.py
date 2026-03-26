from Pipeline import *
import shutil
from multiprocessing import Pool

pipelineobj = VPLModelingPipeline('MEarth', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/TEMPLATES/ModernEarth/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020', planet='Earth')

master_out = '/gscratch/vsm/gialluca/PostDocPropose/Earth_Isotherm_Test/'
atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'

if not os.path.exists(master_out+pipelineobj.casename+'/'):
        os.mkdir(master_out+pipelineobj.casename+'/')

shutil.copytree(atmos_Dir, master_out + pipelineobj.casename + '/atmos/')

pipelineobj.photochemDir = master_out+pipelineobj.casename+'/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
pipelineobj.atmosDir = master_out+pipelineobj.casename+'/atmos/' # path to atmos/ dir
pipelineobj.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
pipelineobj.OutPath = master_out+pipelineobj.casename+'/' # path for the raw model run outputs (NOT for created data products like dictionaries)
pipelineobj.DataOutPath = master_out+pipelineobj.casename+'/' # path for created data products like dictionaries
pipelineobj.AtmProfPath = master_out+pipelineobj.casename+'/' # path to put atmospheric profile files (.pt files really)
pipelineobj.BackupPhotochemRuns = False # Make backups of individual photochem runs
pipelineobj.photochemBackupDir = master_out+pipelineobj.casename+'/PhotochemBackup/' # path to save output from each photochem run
pipelineobj.LBLABC_AbsFilesDir = master_out+pipelineobj.casename+'/ABSFiles/' # path to put the created lbl .abs files in 
pipelineobj.lblabc_RunScriptDir = master_out+pipelineobj.casename+'/' # path to put lbl runscripts in
pipelineobj.vplclimate_RunScriptDir = master_out+pipelineobj.casename+'/' # path to put vpl climate runscripts in
pipelineobj.photochem_InputsDir = master_out+pipelineobj.casename+'/PhotochemInputs/' # The path to create new photochem inputs in
pipelineobj.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found
pipelineobj.SMART_RunScriptDir = master_out+pipelineobj.casename+'/'

pipelineobj.run_spectra = True
pipelineobj.adjust_atmospheric_pressure = False
pipelineobj.rerun_smart_for_2col = False

pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
gas_names = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2']
pipelineobj.molecule_dict['Gas_names'] = gas_names
for m in range(len(gas_names)):
    pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
    pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2


#### Stripped down automatic pipeline run:
###############################################################

if not os.path.exists(pipelineobj.LBLABC_AbsFilesDir):
    os.mkdir(pipelineobj.LBLABC_AbsFilesDir)
#else:
#    subprocess.run('rm -rf '+pipelineobj.LBLABC_AbsFilesDir+'*', shell=True)
# Prepare directory for storing lblabc run script files
if not os.path.exists(pipelineobj.lblabc_RunScriptDir):
    os.mkdir(pipelineobj.lblabc_RunScriptDir)
# Prepare directory for storing vpl climate run script files
if not os.path.exists(pipelineobj.vplclimate_RunScriptDir):
    os.mkdir(pipelineobj.vplclimate_RunScriptDir)
# Prepare directory for storing new photochem inputs
#if not os.path.exists(pipelineobj.photochem_InputsDir):
#    os.mkdir(pipelineobj.photochem_InputsDir)
pipelineobj.setup_intial_photochem_dir()
# Make sure model and data output dirs exist
if not os.path.exists(pipelineobj.OutPath):
    os.mkdir(pipelineobj.OutPath)
if not os.path.exists(pipelineobj.DataOutPath):
    os.mkdir(pipelineobj.DataOutPath)
if not os.path.exists(pipelineobj.SMART_RunScriptDir):
    os.mkdir(pipelineobj.SMART_RunScriptDir)
# Prepare the Hyak environment

pipelineobj.run_photochem_1instance(CleanMake=True, InputCopy=pipelineobj.photochem_InputsDir, trynum=1)

# Get the radius and gravity:

planet = open(pipelineobj.photochemDir+'INPUTFILES/PLANET.dat', 'r')
planet_lines = planet.readlines()
planet.close()
grav = None
rad = None
for i in planet_lines:
    if len(i.split('= G')) > 1:
        grav = float(i.split()[0])*(u.cm*u.s**-2).to(u.m*u.s**-2) #*1e-2 # Get the gravity from the first line of PLANET.dat and convert to m/s**2 (should be originally cm/s**2)
    elif len(i.split('= R0')) > 1:
        rad = float(i.split()[0])*u.cm.to(u.km) #*1e-5 # Get radius from PLANET.dat and convert from cm to km
    elif grav != None and rad != None:
        break
    
# Set object values
pipelineobj.planetary_gravity = grav
pipelineobj.planetary_radius = rad

# Get MMW:
fi = open(pipelineobj.OutPath+'photochem_run_output_'+pipelineobj.casename+'.run', 'r')
lines = fi.readlines()
fi.close()
for i in lines:
    if len(i.split('Molecular weight of atmosphere')) > 1:
        pipelineobj.MMW = float(i.split()[len(i.split())-1])
        break


# Make LBLABC run files
pipelineobj.degrade_PT()
pipelineobj.prep_rmix_file(pipelineobj.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist')

pipelineobj.make_lblabc_runscripts()

lblabc_input = []
for gas in pipelineobj.molecule_dict['Gas_names']:
    lblabc_input.append([pipelineobj.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+pipelineobj.casename+'.script', gas, 'Avg'])

with Pool() as p:
    lblruns = p.map(pipelineobj.run_lblabc_1instance_Parallel, lblabc_input)

pipelineobj.make_smart_runscript()

pipelineobj.run_smart_1instance(pipelineobj.SMART_RunScriptDir+'RunSMART_'+pipelineobj.casename+'.run')