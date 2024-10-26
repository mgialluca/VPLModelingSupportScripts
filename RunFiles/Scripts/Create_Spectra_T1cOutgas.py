import numpy as np
from astropy.io import ascii
from astropy.table import Table
import astropy.units as u
import os, sys, subprocess
from multiprocessing import Pool


use_path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cOutgassing_Cycle4/Debug/'
lblabc_template_path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/T1c0-01barO2-1ppmCO2_hitran2020/'
lblabc_newscripts_path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/T1c0-1barO2_EarthH2OOutgass_hitran2020/Debug/'

mmw_t1c = 31.995902507889991

# Taken from EZ_Photochem
#
### Change the case for the lblabc run
##
## Input:
# template - script to copy, 
# casename - name of the next case you're running (e.g., T1c1barO2-100ppmCO2)
# molecule - molecule of the lblabc script in question (e.g., ch4)
# MMW - atmosphere mean molecular weight (default is keep same from template)
# hitran_gas_code - hitran gas index for the molecule in question (default is same from template) 
# rmix_col - mixing ratio column of molecule in question (default is keep same from template)
# outpath - path where script should be written too
##
## *** Note, outfile (or name of new script) in this case is automatically generated from the case name (less human error)
## *** Note also, this is easiest to run when just using the lblabc script for the molecule as a template...
##     ... but if you want to make it a new molecule this is easily done with the gas code and rmix col options
def lblabc_script_change_case(template, casename, molecule, MMW='same', hitran_gas_code='same', rmix_col='same', 
                              outpath=lblabc_newscripts_path, 
                              absoutpath=use_path):
    f = open(template, 'r')
    l = f.readlines()
    f.close()
    out = open(outpath+'runlblabc_'+molecule+'_'+casename+'_hitran2020.script', 'w')

    out.write(l[0])
    out.write(absoutpath+'PT_profile_'+casename+'.pt\n')

    for i in range(len(l)):
        if i == 0 or i == 1:
            pass
        else:
            if hitran_gas_code != 'same' and len(l[i].split('gas index')) > 1:
                out.write(str(hitran_gas_code)+'                                       hitran gas index\n')
            elif len(l[i].split('/AtmProfiles/')) > 1:
                out.write(absoutpath+'MixingRs_'+casename+'.dat\n')
            elif rmix_col != 'same' and len(l[i].split('columns')) > 1 and len(l[i].split('rmix')) > 1:
                out.write('1,'+str(rmix_col)+'                                     columns of p and rmix\n')
            elif MMW != 'same' and len(l[i].split('mol. wgt.')) > 1:
                out.write(str(MMW)+'                      mol. wgt. of atmosphere (kg/kmole)\n')
            elif len(l[i].split('/LinebyLine_absFiles/')) > 1:
                out.write(absoutpath+casename+'_'+molecule+'_hitran20_50_1e5cm-1.abs\n')
            else:
                out.write(l[i])
    out.close()


# taken from pipeline_raw functions
#
### Take the PT profile output from photochem and degrade it to a specified number of layers ...
### ... to create PT profile for LBLABC and SMART
##
## Inputs:
# casename - name of case on your grid you're doing
# nlevel_new - number of layers you want in your degraded atmosphere
# Prof - PTZ_mixingratios_out WITH PATH output by photochem
# PressUnits - can do Bar or Pa, Megan stick with bar for the forseeable forever
# AtmProfPath - path to output new PT profile to
##
def degrade_PT(casename, nlevel_new=65, ptzFile='/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/T1c/T1cOutgas_Testing/PTZ_mixingratios_out_fordebug.dist', 
               PressUnits='Bar', AtmProfPath=use_path):
    atm = ascii.read(ptzFile, delimiter=' ')
    alt = atm['ALT']
    pres = atm['PRESS']
    temp = atm['TEMP']
    #alt = atm['Alt[km]']
    #new_pres = atm['P[Pa]']*1e-5
    #new_temp = atm['T[K]']

    new_grid = np.linspace(alt[0], alt[len(alt)-1], nlevel_new)
    new_temp = np.interp(new_grid, alt, temp)
    new_pres = np.interp(new_grid, alt, pres)

    if PressUnits == 'Pa':
        new_pres = new_pres*u.bar.to(u.Pa)

    dat = Table([new_pres[::-1], new_temp[::-1]], names=('Press', 'Temp'))
    dat = Table([new_pres, new_temp], names=('Press', 'Temp'))
    ascii.write(dat, AtmProfPath+'PT_profile_'+casename+'.pt', overwrite=True)

### Make pressure increase in the column (reverse all columns) for smart
##
## Inputs:
# casename - name of the atmosphere case on your grid 
# Prof - PTZ_mixingratios_out WITH PATH that was output by photochem
# outputpath - path to write new mixing ratio file to
##
def prep_p_rmix_files_smart(casename, Prof='/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/T1c/T1cOutgas_Testing/PTZ_mixingratios_out_fordebug.dist', 
                            outputpath=use_path, 
                            gases=['O2', 'H2O', 'O3', 'H2']): # 'OH', 'H2', 'HO2', 'H2O2', 'CO', 'O3']):
    atm = ascii.read(Prof, delimiter=' ')
    datfortab = [atm['PRESS'][::-1]]
    #datfortab = [atm['P[Pa]']*1e-5]
    namesfortab = ['Press']
    for i in gases:
        datfortab.append(atm[i])#[::-1])
        namesfortab.append(i)

    dat = Table(datfortab, names=namesfortab)
    ascii.write(dat, outputpath+'MixingRs_'+casename+'.dat', overwrite=True)

# Taken from EZ photochem
#
### Run LBLABC for a given runscript (just useful to get the output w/e)
##
## Inputs:
# runscript - name of LBLABC runscript WITH PATH
# casename - name of case you're running (to name output file)
# outpath - path to put run output in
##
def run_lblabc_1instance(inputs):
    casename, molecule = inputs
    outpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cOutgassing_Cycle4/Debug/'
    scriptpath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/T1c0-1barO2_EarthH2OOutgass_hitran2020/Debug/'
    runscript = scriptpath+'runlblabc_'+molecule+'_'+casename+'_hitran2020.script'

    f = open(outpath+'lblabc_run_output_'+molecule+'_'+casename+'.run', 'w')
    workdir = os.getcwd()
    os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/')
    subprocess.run('/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/lblabc < '+runscript, shell=True, stdout=f)
    os.chdir(workdir)

### Run SMART for a given runscript (just useful to get the output or be fully in python w/e)
##
## Inputs:
# runscript - name of SMART runscript WITH PATH
# casename - name of case you're running (to name output file)
# outpath - path to put run output in
##
def run_smart_1instance(runscript, casename, outpath=use_path):
    workdir = os.getcwd()
    os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/smart/')
    subprocess.run('/gscratch/vsm/gialluca/VPLModelingTools_Dev/smart/smart_spectra < '+runscript+' > '+outpath+'smart_run_output_'+casename+'.run', shell=True)
    os.chdir(workdir)

# Make PT prof and mixing ratios for lblabc, climate & smart preferences

degrade_PT('EarthLikeOutgasT1c_debug')
prep_p_rmix_files_smart('EarthLikeOutgasT1c_debug')

# make lblabc runscripts for relevant constituents

constituents = ['h2o', 'o2', 'o3', 'h2']#, 'oh', 'h2', 'ho2', 'h2o2']

lblabc_script_change_case(lblabc_template_path+'runlblabc_h2o_T1c0-01barO2-1ppmCO2_hitran2020.script', 'EarthLikeOutgasT1c_debug', 
                          'h2o', MMW=mmw_t1c, rmix_col=3)
lblabc_script_change_case(lblabc_template_path+'runlblabc_o2_T1c0-01barO2-1ppmCO2_hitran2020.script', 'EarthLikeOutgasT1c_debug', 
                          'o2', MMW=mmw_t1c, rmix_col=2)
#lblabc_script_change_case(lblabc_template_path+'runlblabc_co_T1c0-01barO2-1ppmCO2_hitran2020.script', 'EarthLikeOutgasT1c_debug', 
#                          'co', MMW=mmw_t1c, rmix_col=8)
lblabc_script_change_case(lblabc_template_path+'runlblabc_o3_T1c0-01barO2-1ppmCO2_hitran2020.script', 'EarthLikeOutgasT1c_debug', 
                          'o3', MMW=mmw_t1c, rmix_col=4)

#lblabc_script_change_case(lblabc_template_path+'runlblabc_o3_T1c0-01barO2-1ppmCO2_hitran2020.script', 'EarthLikeOutgasT1c_debug', 
#                          'oh', MMW=mmw_t1c, hitran_gas_code=13, rmix_col=4)
lblabc_script_change_case(lblabc_template_path+'runlblabc_o3_T1c0-01barO2-1ppmCO2_hitran2020.script', 'EarthLikeOutgasT1c_debug', 
                          'h2', MMW=mmw_t1c, hitran_gas_code=45, rmix_col=5)
#lblabc_script_change_case(lblabc_template_path+'runlblabc_o3_T1c0-01barO2-1ppmCO2_hitran2020.script', 'EarthLikeOutgasT1c_debug', 
#                          'ho2', MMW=mmw_t1c, hitran_gas_code=33, rmix_col=6)
#lblabc_script_change_case(lblabc_template_path+'runlblabc_o3_T1c0-01barO2-1ppmCO2_hitran2020.script', 'EarthLikeOutgasT1c_debug', 
#                          'h2o2', MMW=mmw_t1c, hitran_gas_code=25, rmix_col=7)

# Run lblabc for relevant constituents

lblabc_inputs = [['EarthLikeOutgasT1c_debug', molecule] for molecule in constituents]
with Pool() as p:
    runlblabc = p.map(run_lblabc_1instance, lblabc_inputs)

run_smart_1instance('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/SMART/runsmart_EarthLikeH2OOutgas_T1c0-1barO2_debug.run', 'EarthLikeOutgasT1C_debug_withh2')