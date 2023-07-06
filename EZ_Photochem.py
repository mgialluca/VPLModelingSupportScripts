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

### Taken from vplc_photo.py
## Plots gases of interest from the in.dist and the output after photochem runs
def plot_dists_b4andafter(nz, nq, npar, dz, atmosfile1='./atmos/PHOTOCHEM/OUTPUT/out.dist', atmosfile0='./atmos/PHOTOCHEM/in.dist', atmosfile='./atmos/PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist'):
    #atmosfile1 = './PHOTOCHEM/OUTPUT/out.dist'
    #atmosfile0 = './PHOTOCHEM/in.dist'

    T,edd,den,o3,denco2,aersol,wfall,rpar,mix0 = read_indist(atmosfile0,nz,nq,npar)
    T,edd,den,o3,denco2,aersol,wfall,rpar,mix1 = read_indist(atmosfile1,nz,nq,npar)

    #atmosfile = './PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist'
    f = open(atmosfile)
    ss = f.readline().split()
    gasnames = ss[3:]

    plot_gases = ['O','O2','H2O','H2','CO2','CO','SO4AER', \
                  'SO2','OH','N2O','NO','NO2','O3','NO3','N2O5']

    ng = len(plot_gases)
    ncol = 5
    nrow = 3
    
    fig,axes = plt.subplots(nrow,ncol,figsize=[20,16])

    alt = np.arange(dz,dz*nz+dz,dz)

    print('shape',mix0.shape)
    k = 0

    for i in range(nrow):
        for j in range(ncol):
            
            ax = axes[i,j]
            print('plotting',plot_gases[k])
            ig = gasnames.index(plot_gases[k])
            ax.plot(mix0[ig],alt,color='C3',ls='-')
            ax.plot(mix1[ig],alt,color='k',ls='--')
            ax.set_title(plot_gases[k])
            ax.set_xscale('log')
            minval = np.max(mix0[ig])*1e-6
            minval = np.max([minval,np.min(mix0[ig])])
            ax.set_xlim(left=minval)
            k = k + 1
            if(k >= ng):
                break
    plt.show()

### Taken from vplc_photochem.py
## Take the PT profile output from photochem and degrade it to a specified number of layers
#def degrade_for_smart():

#######################################################################################


### Easy wrapper to just run 1 photochem run
def run_photochem_1instance(CleanMake=True, InputCopy=False, OutPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):

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
        fmake = open(OutPath+'photochem_make_output.txt', 'w')
        os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/')
        subprocess.run('make -f ./PhotoMake clean', shell=True)
        subprocess.run('make -f ./PhotoMake', shell=True, stdout=fmake)
        os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/')

    # Run photochem
    f = open(OutPath+'photochem_run_output.run', 'w')
    os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/')
    subprocess.run('./Photo.run', shell=True, stdout=f)
    os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/')

### Plot the P-T profile and alt vs pressure
def pt_profile(Prof='./atmos/PHOTOCHEM/OUTPUT/profile.pt'):
    atm = ascii.read(Prof, delimiter=' ')
    fig = plt.figure(figsize=(8, 10))
    pt = fig.add_subplot(111)
    palt = pt.twiny()

    palt.plot(atm['Alt']*1e-5, atm['Press'], color='xkcd:blue', label='Alt [km]')
    palt.set_xlabel('Altitude [km]', size=20)
    pt.plot(atm['Temp'], atm['Press'], color='xkcd:orange', label='Temp [K]')
    pt.set_xlabel('Temperature [K]', size=20)
    pt.set_ylabel('Pressure [Bar]', size=20)

    subs = [pt, palt]
    for ax in subs:
        ax.tick_params(length=6, width=2, labelsize=18)
        ax.tick_params(which='minor', width=1, length=4)
        ax.spines['top'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)
        #ax.legend(fontsize=18)
        ax.set_yscale('log')

    pt.legend(fontsize=18, loc='lower right')
    palt.legend(fontsize=18, loc='upper right')
    pt.set_ylim([1e-5, 10])
    pt.set_xlim([220, 760])
    pt.invert_yaxis()
    plt.tight_layout()
    plt.show()

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

### Take the PT profile output from photochem and degrade it to a specified number of layers
def degrade_PT(nlayer, Prof='/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/OUTPUT/profile.pt', outputunits='Bar', outputpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/', outputname='PT_profile_T1c1bar.pt'):
    atm = ascii.read(Prof, delimiter=' ')
    alt = atm['Alt']
    pres = atm['Press']
    temp = atm['Temp']

    new_grid = np.linspace(alt[0], alt[len(alt)-1], nlayer)
    new_temp = np.interp(new_grid, alt, temp)
    new_pres = np.interp(new_grid, alt, pres)

    if outputunits == 'Pa':
        new_pres = new_pres*u.bar.to(u.Pa)

    dat = Table([new_pres[::-1], new_temp[::-1]], names=('Press', 'Temp'))
    ascii.write(dat, outputpath+outputname, overwrite=True)

### Make pressure increase in the column (reverse all columns) for smart
def prep_p_rmix_files_smart(Prof='/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/OUTPUT/profile.pt', outputpath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/', outputname='MixingRs_T1c1bar.dat'):
    atm = ascii.read(Prof, delimiter=' ')
    datfortab = [atm['Press'][::-1]]
    namesfortab = ['Press']
    for i in list(atm.columns):
        if i not in ['Alt', 'Temp', 'Den', 'Press']:
            datfortab.append(atm[i][::-1])
            namesfortab.append(i)

    dat = Table(datfortab, names=namesfortab)
    ascii.write(dat, outputpath+outputname, overwrite=True)


### Make a SMART spectra from the photochem run
#def create_trnst_smart():

### Plot a smart run
def plot_photochems_smart(type='toarad', file='./T1c_1barO2_onlycrosssec_toa.rad', rp=1.097, rs=0.11):
    if type == 'trnst':
        trnst_spec = ascii.read(file)

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111)

        trnst_dep = trnst_spec['col4'] #- (((rp*u.Rearth)/(rs*u.Rsun))**2).decompose()

        ax.plot(trnst_spec['col1'], trnst_dep)
        ax.set_xlabel('Wavelength [um]', size=20)
        ax.set_ylabel('Fractional Transit Depth', size=20)
        
        #if rp is not None and rs is not None:
        #    frac = ((rp*u.Rearth)/(rs*u.Rsun))**2
        #    frac = frac.decompose()
        #    ax.axhline(frac, label='(Rp/Rs)**2')
        #    ax.legend(fontsize=16)
        ax.set_xlim([0,10])

        plt.tight_layout()
        plt.show()

    elif type=='toarad':
        toa = ascii.read(file)

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111)

        ax.plot(toa['col1'], toa['col4'])
        ax.set_xlabel('Wavelength [um]', size=20)
        ax.set_ylabel('Planetary Flux [W/m**2/um]', size=20)
        #ax.set_xlim([0,10])

        plt.tight_layout()
        plt.show()



# Example of wrapper: Run the ModernEarth template
#run_photochem_1instance(CleanMake=True, InputCopy='atmos/PHOTOCHEM/INPUTFILES/TEMPLATES/ModernEarth/')

#run_photochem_1instance(CleanMake=True, InputCopy='./Andrew_T1-c/', OutPath='./T1c_Runs')
#pt_profile()
#prms = basic_params()
#plot_dists_b4andafter(int(prms[0]), int(prms[1]), 0, 0.5)

#degrade_PT(69)