#import matplotlib.pyplot as plt
import numpy as np
import os,sys,shutil,subprocess
from astropy.io import ascii
import sys
import astropy.units as u
import astropy.constants as const
from astropy.table import Table

################################
##
## Big Ol' python wrapper to do various things with photochem
## Author: Megan Gialluca
##
################################

###################### STOLEN FUNCTIONS ADAPTED FOR USE W/O THOSE SCRIPTS/DEPENDENCIES

### Taken from: vplc_photochem.py
## Reads an in.dist/out.dist and returns useful values, mixin are LL species in order of out.params
def read_indist(atmosfile,nz,nq,npar,ncol=10):
    '''Parses the photochem in.dist/out.dist files.
    
    Returns
    -------
    object array
        T,edd,den,O3,denco2,aerosol,wfall,rpar,mixes
    '''

    
    T = np.zeros(nz,dtype=float)
    edd = np.zeros(nz,dtype=float)
    den = np.zeros(nz,dtype=float)
    o3 = np.zeros(nz,dtype=float)
    denco2 = np.zeros(nz,dtype=float)
    aersol = np.empty(npar,dtype=object)
    wfall = np.empty(npar,dtype=object)
    rpar = np.empty(npar,dtype=object)

    for i in range(npar):
        aersol[i] = np.zeros(nz,dtype=float)
        wfall[i] = np.zeros(nz,dtype=float)
        rpar[i] = np.zeros(nz,dtype=float)

    mixin = np.empty(nq,dtype=object)
    for i in range(nq):
        mixin[i] = np.zeros(nz,dtype=float)

    f = open(atmosfile,'r')

    ncols = nq%ncol

    # Gases
    for n in range(nq//ncol):
        for k in range(nz):
            ss = f.readline().split()
            for i in range(ncol):
                try:
                    mixin[i+n*ncol][k] = ss[i]
                except:
                    mixin[i+n*ncol][k] = 0.
                if(mixin[i+n*ncol][k] < 1.e-99): mixin[i+n*ncol][k] = 0.
    n = nq//ncol
    if(ncols >= 1):
        for k in range(nz):
            ss = f.readline().split()
            for i in range(ncols):
                try:
                    mixin[i+n*ncol][k] = ss[i]
                except:
                    mixin[i+n*ncol][k] = 0.
                if(mixin[i+n*ncol][k] < 1.e-99): mixin[i+n*ncol][k] = 0.
        
    # Other info (sec 2)
    for k in range(nz):
        ss = f.readline().split()
        T[k] = ss[0]
        edd[k] = ss[1]
        den[k] = ss[2]
        o3[k] = ss[3]
        denco2[k] = ss[4]
    
    # Sec 3 aerosols
    for k in range(nz):
        ss = f.readline().split()
        for p in range(npar):
            aersol[p][k] = ss[p*3]
            wfall[p][k] = ss[p*3+1]
            rpar[p][k] = ss[p*3+2]
        
    f.close()

    return T,edd,den,o3,denco2,aersol,wfall,rpar,mixin


#######################################################################################


### Easy wrapper to just run 1 photochem run for Megan's Klone dev environment
##
## Inputs:
# casename - name of case you're runnin (to name make and run output files)
# CleanMake - do a clean make of photochem or no
# InputCopy - the path to the input photochem files, if false assumes they are already in the PHOTOCHEM/INPUTS directory
# OutPath - path to write model make and run outputs to
##
def run_photochem_1instance(casename, CleanMake=True, InputCopy=False, OutPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):

    # If you have new input files to use, give 'InputCopy' the dir path
    if InputCopy != False:
        if CleanMake != True:
            choice = input('Youve selected new inputs but do not want a clean make, this is typically incorrect. Do you still want to continue? [y/n] \n')
            if choice == 'n' or choice == 'N':
                sys.exit('Try to Clean Make with new inputs')
            elif choice == 'y' or choice == 'Y':
                pass
            else:
                choice = input('I SAID continue? - """ y """ or """ n """"\n')
                if choice == 'n' or choice == 'N':
                    sys.exit('Try to Clean Make with new inputs')
                elif choice == 'y' or choice == 'Y':
                    pass
                else:
                    sys.exit('User doesnt follow instructions, terminating')

        # Remove old input files if they exist
        subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/input_photchem.dat', shell=True)
        subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/parameters.inc', shell=True)
        #subprocess.run('rm ./atmos/PHOTOCHEM/INPUTFILES/params.dat', shell=True)
        subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/PLANET.dat', shell=True)
        subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/reactions.rx', shell=True)
        subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/species.dat', shell=True)
        subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/in.dist', shell=True)

        # Copy new input files to the right places
        subprocess.run('cp '+InputCopy+'input_photchem.dat /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'parameters.inc /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/', shell=True)
        #subprocess.run('cp '+InputCopy+'params.dat ./atmos/PHOTOCHEM/INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'PLANET.dat /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'reactions.rx /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'species.dat /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'in.dist /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/', shell=True)

    # Clear the outputs
    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/OUTPUT/*', shell=True)
    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/PTZ_mixingratios_in.dist', shell=True)

    # Clean make, if requested
    if CleanMake:
        fmake = open(OutPath+'photochem_make_output_'+casename+'.txt', 'w')
        os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/')
        subprocess.run('make -f ./PhotoMake clean', shell=True)
        subprocess.run('make -f ./PhotoMake', shell=True, stdout=fmake)
        os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/')

    # Run photochem
    f = open(OutPath+'photochem_run_output_'+casename+'.run', 'w')
    os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/')
    subprocess.run('./Photo.run', shell=True, stdout=f)
    os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/')

### Run LBLABC for a given runscript (just useful to get the output w/e)
##
## Inputs:
# runscript - name of LBLABC runscript WITH PATH
# casename - name of case you're running (to name output file)
# outpath - path to put run output in
##
def run_lblabc_1instance(runscript, casename, outpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):
    f = open(outpath+'lblabc_run_output_'+casename+'.run', 'w')
    subprocess.run('/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/lblabc < '+runscript, shell=True, stdout=f)

### Run SMART for a given runscript (just useful to get the output or be fully in python w/e)
##
## Inputs:
# runscript - name of SMART runscript WITH PATH
# casename - name of case you're running (to name output file)
# outpath - path to put run output in
##
def run_smart_1instance(runscript, casename, outpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):
    subprocess.run('/gscratch/vsm/gialluca/VPLModelingTools_Dev/smart/smart_spectra < '+runscript+' > '+outpath+'smart_run_output_'+casename+'.run', shell=True)

### Get the NZ, NQ, NQ1, NSP2, NR, KJ, NP from out.params
def basic_params(Prms='./atmos/PHOTOCHEM/OUTPUT/out.params'):
    prms = open(Prms, 'r')
    fstl_unformatted = prms.readlines()[0]
    fstl = []
    for i in fstl_unformatted.split(' '):
        if i != '':
            fstl.append(i)
    fstl[len(fstl)-1] = fstl[len(fstl)-1].split('\n')[0]
    return fstl

### Take the PT profile output from photochem and degrade it to a specified number of layers ...
### ... to create PT profile for LBLABC and SMART
##
## Inputs:
# nlayer - number of layers you want in your degraded atmosphere
# casename - name of case on your grid you're doing
# Prof - PTZ_mixingratios_out WITH PATH output by photochem
# outputunits - can do Bar or Pa, Megan stick with bar for the forseeable forever
# outputpath - path to output new PT profile to
##
def degrade_PT(nlayer, casename, Prof='/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist', 
               outputunits='Bar', outputpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/'):
    atm = ascii.read(Prof, delimiter=' ')
    alt = atm['ALT']
    pres = atm['PRESS']
    temp = atm['TEMP']

    new_grid = np.linspace(alt[0], alt[len(alt)-1], nlayer)
    new_temp = np.interp(new_grid, alt, temp)
    new_pres = np.interp(new_grid, alt, pres)

    if outputunits == 'Pa':
        new_pres = new_pres*u.bar.to(u.Pa)

    dat = Table([new_pres[::-1], new_temp[::-1]], names=('Press', 'Temp'))
    ascii.write(dat, outputpath+'PT_profile_'+casename+'.pt', overwrite=True)

### Make pressure increase in the column (reverse all columns) for smart
##
## Inputs:
# casename - name of the atmosphere case on your grid 
# Prof - PTZ_mixingratios_out WITH PATH that was output by photochem
# outputpath - path to write new mixing ratio file to
##
def prep_p_rmix_files_smart(casename, Prof='/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist', 
                            outputpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/', 
                            gases=['O2', 'H2O', 'O3', 'CO2', 'CO', 'SO2', 'H2SO4', 'N2O', 'NO2', 'HNO3']):
    atm = ascii.read(Prof)
    datfortab = [atm['PRESS'][::-1]]
    namesfortab = ['Press']
    for i in gases:
        datfortab.append(atm[i][::-1])
        namesfortab.append(i)

    dat = Table(datfortab, names=namesfortab)
    ascii.write(dat, outputpath+'MixingRs_'+casename+'.dat', overwrite=True)


############## BELOW FUNCTIONS FOR EDITING/CREATING LBLABC SCRIPTS ###########

### Switch the hitran used by an lblabc runscript
# Actualy fuck this, I'll just use HITRAN2020 4ever
#def lblabc_script_change_hitran(template, outfile, hitran=2020, outpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/'):

### Where you left off:
## The below finds the hitran relevant lines of a lblabc run script
## Your goal is to rewrite an identical runscript that changes hitran to 2016 or 2020

#    Careful on spaces versus tabs (\t)
#    for i in l:
#        if len(i.split(' ')) > 1 and len(i.split('HITRAN')) > 1:
#            print('Guessing this is the integer hitran line:')
#            print(i)
#            print('\n')
#        elif len(i.split(' ')) == 1 and len(i.split('HITRAN')) > 1:
#            print('Guessing this is the file path to .par:')
#            print(i)
#            print('\n')
#        elif len(i.split('fundam')) > 1:
#            print('Guessing this is the fundamental dat file path:')
#            print(i)
#            print('\n')
#        elif len(i.split('hitran')) > 1 and len(i.split(' ')) == 1:
#            print('Guessing this is the output file name:')
#            print(i)
#            print('\n')

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
                              outpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/'):
    f = open(template, 'r')
    l = f.readlines()
    f.close()
    out = open(outpath+'runlblabc_'+molecule+'_'+casename+'_hitran2020.script', 'w')

    out.write(l[0])
    out.write('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/PT_profile_'+casename+'.pt\n')

    for i in range(len(l)):
        if i == 0 or i == 1:
            pass
        else:
            if hitran_gas_code != 'same' and len(l[i].split('gas index')) > 1:
                out.write(str(hitran_gas_code)+'                                       hitran gas index\n')
            elif len(l[i].split('/AtmProfiles/')) > 1:
                out.write('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/MixingRs_'+casename+'.dat\n')
            elif rmix_col != 'same' and len(l[i].split('columns')) > 1 and len(l[i].split('rmix')) > 1:
                out.write('1,'+str(rmix_col)+'                                     columns of p and rmix\n')
            elif MMW != 'same' and len(l[i].split('mol. wgt.')) > 1:
                out.write(str(MMW)+'                      mol. wgt. of atmosphere (kg/kmole)\n')
            elif len(l[i].split('/LinebyLine_absFiles/')) > 1:
                out.write('/gscratch/vsm/gialluca/VPLModelingTools_Dev/LinebyLine_absFiles/'+casename+'_'+molecule+'_hitran20_50_1e5cm-1.abs\n')
            else:
                out.write(l[i])
    out.close()


############## BELOW FUNCTIONS FOR EDITING/CREATING SMART SCRIPTS ###########

### Add a stellar SED
##
## Input:
# noSEDfi - script u want to copy without an stellar SED
# SED - solar SED file
# outfile - name of output script WITHOUT path
# units - flux units of the stellar SED (default [W/m**2/cm])
# skiplines - # of lines smart to skip in SED file
# cols_wn_flux - columns of wavenumber & flux in SED file (1 indexed)
# outpath - path where output script should be written to
##
def smart_script_addSED(noSEDfi, SED, outfile, units=1, skiplines=1, cols_wn_flux=[1,2], 
                        outpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/SMART/'):
    f = open(noSEDfi, 'r')
    l = f.readlines()
    f.close()
    out = open(outpath+outfile, 'w')

    for i in range(len(l)):
        if len(l[i].split('Source')) > 1 or len(l[i].split('source')) > 1:
            out.write('3\t\t\tSource - Solar\n')
            out.write('3\t\t\tFile Format - List Directed\n')
            out.write(SED+'\n')
            out.write(str(skiplines)+'\t\t\tLines to Skip\n')
            out.write(str(units)+'\t\t\tUnits solar flux\n')
            out.write('1\t\t\tSolar spectral units\n')
            out.write('1.0\t\t\tMicron Conversion Factor\n')
            out.write(str(cols_wn_flux[0])+','+str(cols_wn_flux[1])+' \t\tColumns of wn and Flux\n')
            out.write('1\t\t\tNumber of Solar Zenith Angles\n')
            out.write('60,0\t\t\tZenith and Azimuth Angles\n')
            out.write('0.01\t\t\tConvergence Criterion\n')
        else:
            out.write(l[i])

    out.close()

### Change the case for the SMART run, dumb function version
### This will just take a template and replace all the case names with options to change the:
#### surface temp, albedo file (but NOT the lines to skip, columns, grid type, or scaling of that file), ...
#### ... and MMW of atmosphere
### IF YOU WANT TO CHANGE THE ABSORBERS, their rmix columns or anything like that, use the verbose function
## 
## Input:
# template - script to copy and replace case name in
# casename - name of next grid case you're running
# Tsurf - surface temperature [K] (default same as template)
# AlbedoFi - Albedo File WITH FULL PATH (default same as template)
# MMW - atmosphere mean molecular weight (default same as template)
# LBL_identifier - identifier of LBLABC abs files that should be identical for every absorber ...
#     ... after that absorbers gas name (e.g., _hitran20_50_1e5cm-1.abs) > if ur megan this should be the same most of the time
# outpath - path where output script should be written to
##
## *** Note, outfile (or name of new script) in this case is automatically generated from the case name (less human error)
## *** Note also, this quick function really only works for the file management structure Megan impliments in her Klone environment
def smart_script_change_case_quick(template, casename, Tsurf='same', AlbdedoFi='same', MMW='same', 
                                   LBL_identifier='_hitran20_50_1e5cm-1.abs', outpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/SMART/'):
    f = open(template, 'r')
    l = f.readlines()
    f.close()
    out = open(outpath+'runsmart_'+casename+'.run', 'w')

    for i in range(len(l)):
        if len(l[i].split('PT_profile')) > 1:
            a = l[i].split('PT_profile')
            out.write(a[0]+'PT_profile_'+casename+'.pt\n')
        elif Tsurf != 'same' and len(l[i].split('Surface Temp')) > 1:
            out.write(str(Tsurf)+' \t\tSurface Temperature\n')
        elif len(l[i].split('LinebyLine_absFiles')) > 1:
            a = l[i].split(LBL_identifier)
            b = a[0].split('_')
            gas = b[len(b)-1]
            out.write('/gscratch/vsm/gialluca/VPLModelingTools_Dev/LinebyLine_absFiles/'+casename+'_'+gas+LBL_identifier+'\n')
        elif len(l[i].split('MixingRs')) > 1:
            out.write('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/MixingRs_'+casename+'.dat\n')
        elif AlbdedoFi != 'same' and len(l[i].split('/albedo/')) > 1:
            out.write(AlbdedoFi+'\n')
        elif MMW != 'same' and len(l[i].split('Mean Molecular Weight')) > 1:
            out.write(str(MMW)+'\t\tMean Molecular Weight\n')
        elif len(l[i].split('/Spectra/')) > 1:
            out.write('/gscratch/vsm/gialluca/VPLModelingTools_Dev/Spectra/'+casename+'\n')
        else:
            out.write(l[i])
    out.close()


# Function to edit to just call on hyak lol
def go_go_hyak_do_it(case, mmw):
    flag = True
    if flag == True:
        lblrunscriptpath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/'
        smartrunscriptpath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/SMART/'
        gases = ['o2', 'h2o', 'o3', 'co2', 'co', 'so2', 'n2o', 'no2', 'hno3']
        gases_capp = ['O2', 'H2O', 'O3', 'CO2', 'CO', 'SO2', 'N2O', 'NO2', 'HNO3']
        rmixcols = [2, 3, 4, 5, 6, 7, 9, 10, 11]

        for i in range(len(gases)):
            lblabc_script_change_case(lblrunscriptpath+'T1cFidO21bar500ppmCO2_hitran2020/runlblabc_'+gases_capp[i]+'_T1cFidO21bar_hitran2020.script', case, gases[i], MMW=mmw, rmix_col=rmixcols[i])

        smart_script_change_case_quick(smartrunscriptpath+'runsmart_T1c_TRUEFIDUCIAL.run', case, MMW=mmw)

### Write a new smart run script for a new case with a template file for moderate guidance
### Template file should be everything u want but the [required input] absorber descriptions, PT profile, ...
### ... [optional input] pressure scaling, pressure column, temp column in pt prof, surf temp, albedo file, ...
### ... albedo file description, MMW of atmosphere, & Min/Max wavenumber
### (if it doesnt have an SED and u want one, just do this and then use the addSED function above)
### Note the dictionary option also allows u to specify every input in a dictionary with keys given by the input names
##
## Input:
## [REQUIRED INPUTS]
# template - script to copy and modify, defines defaults
# casename - name of the next grid case you're running
# PTProf - P-T profile file WITH PATH
# MixingRsFile - File with mixing ratios WITH PATH
# LBLabsFiles_path - path to where your LBL abs files live
# LBLabsFiles - array of LBL abs file names
# xsec_path - path to where your cross section files live
# xsecFiles - array of xsec files (should correspond to same gases in LBL array)
# hitran_codes - array of hitran codes defining order of absorbers
# rmix_cols - columns of rmix for each absorber in the MixingRsFile
# abs_types - number of absorption types for each absorber (1 - xsec only, 2 - both xsec and LBL, 3 - LBL only)
## [OPTIONAL INPUTS]
# P_col - column of Pressure in the PT and MixingRs files (default is template)
# T_col_inPTProf - column of temperature in the PT profile
# Pscaling - pressure scaling factor
# Tsurf - surface temperature
# AlbedoFile - albedo file WITH PATH
# SkipLines_Albedo - # of lines to skip in albedo file
# WLAlbedoCols - columns of wavelength and albedo
# wlGridType_albedo - wavelength grid type in albedo file
# wlScaling_albedo - wavelength scaling in albedo file
# MMW - mean molecular weight of atmosphere
# outpath - path to write output script to
##
#def smart_script_change_case_verbose(dict_defined=False, template=None, casename=None, PTProf=None, MixingRsFile=None, P_col='same', 
##                                     T_col_inPTProf='same', Pscaling='same', LBLabsFiles_path=None, LBLabsFiles=None, xsec_path=None, 
#                                     xsecFiles=None, hitran_codes=None, rmix_cols=None, abs_types=None, Tsurf='same', AlbedoFile='same', 
#                                     SkipLines_Albedo='same', WlAlbedoCols='same', wlGridType_albedo='same', wlScaling_albedo='same', 
#                                     MMW='same', outpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/SMART/')
#
# NOTE NOTE I started doing this then decided it wasn't worth the effort yet because of the grid I'm trying to run
# will probably come back and write this in the future





# Example of wrapper: Run the ModernEarth template
#run_photochem_1instance(CleanMake=True, InputCopy='atmos/PHOTOCHEM/INPUTFILES/TEMPLATES/ModernEarth/')
