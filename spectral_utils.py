import numpy as np

### SMART-related functions

def read_hrt(filename,nlev,nsza,intsza=False,retsza=False):
    '''
    Reads and parses the VPL Climate heating rate files.
    '''
    
    p_q = np.zeros(nlev)
    alt_q = np.zeros(nlev)
    t_q = np.zeros(nlev)
    p = np.zeros(nlev)
    alt = np.zeros(nlev)
    t = np.zeros(nlev)
    
    q_t = np.zeros(nlev)
    ft_dn = np.zeros(nlev)
    ft_up = np.zeros(nlev)
    ft_net = np.zeros(nlev)

    q_s = np.zeros([nsza,nlev])
    fs_dir = np.zeros([nsza,nlev])
    fs_dn = np.zeros([nsza,nlev])
    fs_up = np.zeros([nsza,nlev])
    fs_net = np.zeros([nsza,nlev])
    
    q_t_net = np.zeros(nlev)
    q_s_net = np.zeros(nlev)
    fs_avg = np.zeros(nlev)

    ft_net_avg = np.zeros(nlev)
    fs_net_avg = np.zeros(nlev)

    
    nskip = 5
    
    f = open(filename,'r')
    
    #### THERMAL ####
    
    for i in range(nskip): f.readline()
        
    for i in range(nlev):
        ss = f.readline().split()
        p_q[i] = ss[0]
        alt_q[i] = ss[1]
        t_q[i] = ss[2]
        q_t[i] = ss[3]
        ft_dn[i] = ss[4]
        ft_up[i] = ss[5]
        ft_net[i] = ss[6]
        p[i] = ss[7]
        alt[i] = ss[8]
        t[i] = ss[9]
        
    
    #### SOLAR ####
    for n in range(nsza):
        for i in range(nskip): f.readline()
            
        for i in range(nlev):
            ss = f.readline().split()
            q_s[n,i] = ss[3]
            fs_dir[n,i] = ss[4]
            fs_dn[n,i] = ss[5]
            fs_up[n,i] = ss[6]
            fs_net[n,i] = ss[7]
            
    #### NET ####
    
    for i in range(nskip): f.readline()
        
    for i in range(nlev):
        ss = f.readline().split()
        q_t_net[i] = ss[3]
        q_s_net[i] = ss[4]
        fs_net_avg[i] = ss[6]
        ft_net_avg[i] = ss[7]

    f.close()

    #### INTEGRATE SZAS ####

    if(intsza and nsza > 1):

        print('Integrating/averaging SZA data')

        umu,wgts = np.polynomial.legendre.leggauss(nsza)
        umu = umu*0.5 + 0.5
        wgts = wgts*0.5
        print('gaussian points for n =',nsza,':',umu)
        print('gaussian weights for n =',nsza,':',wgts)

        fs_dir_avg = np.zeros(nlev)
        fs_dn_avg = np.zeros(nlev)
        fs_up_avg = np.zeros(nlev)

        for n in range(nsza):
            fs_dir_avg = fs_dir_avg + fs_dir[n,:]*wgts[n]
            fs_dn_avg = fs_dn_avg + fs_dn[n,:]*wgts[n]
            fs_up_avg = fs_up_avg + fs_up[n,:]*wgts[n]

        
    elif(intsza and nsza == 1):
        fs_dir_avg = fs_dir[0]
        fs_dn_avg = fs_dn[0]
        fs_up_avg = fs_up[0]
    elif(intsza):
        fs_dir_avg = fs_dir
        fs_dn_avg = fs_dn
        fs_up_avg = fs_up
        
    if(intsza):
        retval = p_q,alt_q,t_q,p,alt,t,ft_dn,ft_up,ft_net,fs_dir_avg,fs_dn_avg,fs_up_avg,fs_net_avg,q_t_net,q_s_net
    else:
        retval = p_q,alt_q,t_q,p,alt,t,ft_dn,ft_up,ft_net,fs_dir,fs_dn,fs_up,fs_net,q_t_net,q_s_net
        
    return retval
    
# load VPL climate spectral outputs (binary)

import numpy as np
import sys

def get_spec(prefix,res,wnt_lim,wns_lim,nlev,nsza=1,plotsza=False,get_solar=False,get_thermal=False,write_solar=False,write_thermal=False,trim_zeros=True):
    '''
    Gets binary VPL climate model spectrum for solar and thermal and combines them.
    
    tflx : thermal emitted flux in thermal wavelengths (nt)
    eflx : Solar outgoing flux in Solar wavelengths (ns)
    sflx : Solar TOA irradiation in Solar wavelengths (ns)
    
    s_tflx : thermal emitted flux in Solar wavelengths (ns)
    t_eflx : Solar outgoing flux in thermal wavelengths (nt)
    t_sflx : Solar TOA irradiation in thermal wavelengths (nt)
    
    res : thermal, solar resolution in wn
    '''
    
    if(nsza >= 1):
        # get gaussian quadrature weights
        umu,wgts = np.polynomial.legendre.leggauss(nsza)
        umu = umu*0.5 + 0.5
        wgts = wgts*0.5
        print('gaussian points for n =',nsza,':',umu)
        print('gaussian weights for n =',nsza,':',wgts)
        ns = int((wns_lim[1]-wns_lim[0])//res[1])

    nt = int((wnt_lim[1]-wnt_lim[0])//res[0])
    
    
    #### FLUX ####
    # Load flux files - climate

    print('plot_spec opening files ',prefix)

    sfix = '.sflx'
    tfix = '.tflx'


    wlt = np.zeros(nt)
    wnt = np.zeros(nt)
    albt = np.zeros(nt)
    tflx = np.zeros(nt)
    
    if(nsza >= 1):
        wlsn = np.zeros([nsza,ns])
        wnsn = np.zeros([nsza,ns])
        eflxn = np.zeros([nsza,ns])
        sflxn = np.zeros([nsza,ns])
        wls = np.zeros(ns)
        wns = np.zeros(ns)
        eflx = np.zeros(ns)
        sflx = np.zeros(ns)

    # Load binary thermal flux data
    tdata = np.fromfile(prefix+tfix,dtype=np.float32)

    # Write ASCII thermal flux data
    if(write_thermal):
        f = open(prefix+'_tflx.txt','w')
        for i in range(nt):
            try:
                f.write('{:12.4e} {:12.4e} {:12.4e} {:12.4e}\n'.format(tdata[1+6*i],tdata[2+6*i],tdata[3+6*i],tdata[4+6*i]))
            except:
                print('could not print next wn,wl',tdata[2+6*(i-1)],tdata[1+6*(i-1)])
                nt = i-1
                break
        f.close()

    if(nsza >= 1):
        sdata = np.empty(nsza,dtype=object)
        for i in range(nsza):
            print('opening sflx'+'{:02d}'.format(i+1))
            sdata[i] = np.fromfile(prefix+sfix+'{:02d}'.format(i+1),dtype=np.float32)

    # parse data
    for i in range(nt):
        try:
            wlt[i] = tdata[1+6*i]
            wnt[i] = tdata[2+6*i]
            albt[i] = tdata[3+6*i]
            tflx[i] = tdata[4+6*i]
        except:
            print('Broke parsing thermal at wn =',wnt[i-1])
            nt = i-1
            break

    wlt = np.trim_zeros(wlt,trim='b')
    #wnt = np.trim_zeros(wnt,trim='b')
    #albt = np.trim_zeros(albt,trim='b')
    #tflx = np.trim_zeros(tflx,trim='b')
    nt = len(wlt)
    wnt = wnt[:nt]
    albt = albt[:nt]
    tflx = tflx[:nt]

    print(len(wlt),len(wnt),len(albt),len(tflx))
    
    if(nsza >= 1):
        for i in range(ns):
            try:
                wls[i] = sdata[0][1+6*i]
                for j in range(nsza):
                    wlsn[j,i] = sdata[j][1+6*i]
            except:
                break
            wns[i] = sdata[0][2+6*i]
            for j in range(nsza):
                wnsn[j,i] = sdata[j][2+6*i]

            for j in range(nsza):
                sflxn[j,i] = sdata[j][3+6*i]
                eflxn[j,i] = sdata[j][4+6*i]
                
            eflx[i] = np.sum(eflxn[:,i]*wgts)
            sflx[i] = np.sum(sflxn[:,i]*wgts)

            #if(trim_zeros):
        wls = np.trim_zeros(wls,trim='b')
        
        #wns = np.trim_zeros(wns,trim='b')
        #eflx = np.trim_zeros(eflx,trim='b')
        #sflx = np.trim_zeros(sflx,trim='b')
        ns = len(wls)
        wns = wns[:ns]
        eflx = eflx[:ns]
        sflx = sflx[:ns]
                    
    # Write ASCII thermal flux data
    if(write_solar and nsza>=1):
        for n in range(nsza):
            f = open(prefix+'_sflx'+str(n)+'.txt','w')
            for i in range(ns):
                try:
                    f.write('{:12.4e} {:12.4e} {:12.4e} {:12.4e}\n'.format(sdata[n][1+6*i],sdata[n][2+6*i],sdata[n][3+6*i],sdata[n][4+6*i]))
                except:
                    print('could not print next wn,wl',sdata[n][2+6*(i-1)],sdata[n][1+6*(i-1)])
                    f.close()
                    break
            f.close()
        f = open(prefix+'_sflx.txt','w')
        for i in range(ns):
            try:
                f.write('{:12.4e} {:12.4e} {:12.4e} {:12.4e}\n'.format(wls[i],wns[i],sflx[i],eflx[i]))
            except:
                print('could not print next wn,wl',wls[i-1],wns[i-1])
                f.close()
                break
        f.close()
            


    # Cross-load fluxes

    if(nsza >=1):
        s_tflx = np.zeros(ns)
        for i in range(ns):
            for j in range(nt):
                if(wlt[j] == wls[i]): 
                    s_tflx[i] = tflx[j]
                    break

        try:
            t_eflx = np.interp(wnt,wns,eflx)
        except:
            print('error in interp',wnt,wns,eflx,len(wnt),len(wns),len(eflx))
            sys.exit()
        try:            
            t_sflx = np.interp(wnt,wns,sflx)
        except:
            print('error in interp',wnt,wns,sflx,len(wnt),len(wns),len(sflx))
            sys.exit()


        s_refl = (eflx + s_tflx)/sflx
        flux = tflx + t_eflx
        fluxs = eflx + s_tflx

    
    if(get_solar):
        return wns,wls,sflx,fluxs
    elif(get_thermal):
        return wnt,wlt,tflx
    else:
        return wnt,wlt,t_sflx,flux
    
    
def open_MIRI_filters(path="/home/mgialluca/Nextcloud/VPL_Modeling/Atmos_Dev/Andrew_HelperFunctions/MIRI.xlsx"):
    import pandas as pd
    df_filt = pd.read_excel(path)
    vals = np.array(df_filt.values[1:,:], dtype=float)
    wl_filt = vals[:,0]
    MIRI_filters = vals[:,1:]
    MIRI_filters[MIRI_filters < 0.001] = 0.0
    filt_names = np.array(df_filt.columns[1:], dtype=str)
    return MIRI_filters, wl_filt, filt_names

def plot_miri_filters(wl, filters, names, ax=None, ylim=[0.0,1.0],same=True):
    if ax is None:
        fig, ax1 = plt.subplots(figsize=(16,7))
        rotation=90
        labelpad=None
    elif(same == False):
        ax1 = ax.twinx()
        rotation=270
        labelpad=25
    else:
        ax1 = ax
        rotation=270
        labelpad=25
    for i in range(len(names)):
        # Plot filter throughputs
        ax1.fill_between(wl, filters[:,i], label=names[i], alpha=0.3)
        ax1.set_ylabel("Filter Throughput", rotation=rotation, labelpad=labelpad)
    leg=ax1.legend(loc=2, fontsize=16, ncol=1)
    leg.get_frame().set_alpha(0.0)
    ax1.set_ylim(ylim)
    plt.show()

def rebin_spec(specHR,lamHR,dlam=None,bin_edges=None):
    
    """
    Re-bin spectum to lower resolution by trapezoidal integration.

    Parameters
    ----------
    specHR : array-like
        Spectrum to be degraded
    lamHR : array-like
        High-res wavelength grid
    lamLR : array-like
        Low-res wavelength grid
    dlam : array-like, optional
        Low-res wavelength width grid
    bin_edges : array-like, optional
        Low-res wavelength grid edges, should be len(lamLR)+1

    Returns
    -------
    specLR : :py:obj:`numpy.ndarray`
        Low-res spectrum
    """

    if dlam is None and bin_edges is None:
        ValueError("Please supply dlam or bin_edges in rebin_spec()")

    elif(dlam is None):
        opt = 'dlam'
    else:
        opt = 'edges'
        
    nbins = len(bin_edges)-1
    
    specLR = np.zeros(nbins)
    
    for i in range(nbins):
        mask = (lamHR >= bin_edges[i]) & (lamHR <= bin_edges[i+1])
        if(np.sum(mask) == 1):
            specLR[i] = specHR[mask] #/(bin_edges[i]-bin_edges[i+1])
        elif(np.sum(mask) < 1):
            specLR[i] = np.interp(0.5*(bin_edges[i]+bin_edges[i+1]),lamHR,specHR)
        elif(np.max(lamHR[mask]) == np.min(lamHR[mask])):
            specLR[i] = np.average(specHR[mask])
        else:            
            #specLR[i] = np.trapz(specHR[mask],lamHR[mask])/(bin_edges[i+1]-bin_edges[i])
            specLR[i] = np.trapz(specHR[mask],lamHR[mask])/(np.max(lamHR[mask])-np.min(lamHR[mask]))
            
        
    return specLR

################################################################################
#### Radiance Transmittance
################################################################################
def read_radtrn(path, Nstr = 8, nlout = 1):
    """
    Reads SMART `*.trn` (transmission) output files
    Note
    ----
    These files only get produced by SMART when specifically asked for
    e.g. ``smartin.out_format = 2``
    Parameters
    ----------
    path : str
        Pathname for *.trn file
    Nstr : int
        Number of RT streams
    nlout : int
        Number of output levels
    Returns
    -------
    lam, wno, trn_ray_0, P_ray_0, P_ray_str, \
                     trn_gas_0, P_gas_0, P_gas_str, \
                     trn_aero_0, P_aero_0, P_aero_str, \
                     trn_all_0
    """
    # Convert each line to vector, compose array of vectors
    arrays = np.array([np.array(list(map(float, line.split()))) for line in open(path)])
    # Flatten into a vector
    arr = np.hstack(arrays)
    # Derive number of columns
    dim = 2 + 3 * (2 + int(Nstr / 2))
    # Compose rectangular grid using known number of columns
    rect = arr.reshape((dim, -1), order='F')
    # Parse file:
    # Wavelength and wavenumbers always first
    lam = rect[0,:]
    wno = rect[1,:]
    # Rayleigh scattering first
    trn_ray_0 = rect[2, :]
    P_ray_0 =   rect[3, :] * 1e5
    P_ray_str = rect[4:4+int(Nstr/2), :] * 1e5
    # Gas absorption second
    trn_gas_0 = rect[4+int(Nstr/2), :]
    P_gas_0 =   rect[4+int(Nstr/2)+1, :] * 1e5
    P_gas_str = rect[4+int(Nstr/2)+2:4+2*int(Nstr/2)+2, :] * 1e5
    # Aerosols third
    trn_aero_0 = rect[4+2*int(Nstr/2)+2, :]
    P_aero_0 =   rect[4+2*int(Nstr/2)+3, :] * 1e5
    P_aero_str = rect[4+2*int(Nstr/2)+4:4+3*int(Nstr/2)+4, :] * 1e5
    # Calculate total
    trn_all_0 = trn_ray_0 * trn_gas_0 * trn_aero_0
    return lam, wno, trn_ray_0, P_ray_0, P_ray_str, \
                     trn_gas_0, P_gas_0, P_gas_str, \
                     trn_aero_0, P_aero_0, P_aero_str, \
                     trn_all_0
class RadTrn(object):
    """
    SMART RadTrn object to contain all SMART `*.trn` (transmission) output
    files.
    Note
    ----
    These files only get produced by SMART when specifically asked for
    e.g. ``smartin.out_format = 2``.
    Parameters
    ----------
    path : str
        Pathname for *.trn file
    Nstr : int
        Number of RT streams
    Attributes
    ----------
    lam : array-like
        Wavelength grid [microns]
    wno : array-like
        Wavenumber grid [cm$^{-1}$]
    trn_ray_0 : array-like
        Transmittance due to Rayleigh scattering
    P_ray_0 : array-like
        Pressure [Pa] at tau=1 due to Rayleigh scattering at normal incidence
    P_ray_str : array-like
        Pressure [Pa] at tau=1 due to Rayeigh scattering for upwelling streams
    trn_gas_0 : array-like
        Transmittance due to gas opacity
    P_gas_0 : array-like
        Pressure [Pa] at tau=1 due to gas opacity at normal incidence
    P_gas_str : array-like
        Pressure [Pa] at tau=1 due to gas opacity for upwelling streams
    trn_areo_0 : array-like
        Transmittance due to aerosols
    P_aero_0 : array-like
        Pressure [Pa] at tau=1 due to aerosols at normal incidence
    P_aero_str : array-like
        Pressure [Pa] at tau=1 due to aerosols for upwelling streams
    trn_all_0 : array-like
        Total transmittance
    """
    def __init__(self, path, Nstr = 8):
        self._path = path
        self.Nstr = Nstr
        self._opened = False
        if path is not None: self._open_path(path)
    @property
    def path(self):
        return self._path
    @path.setter
    def path(self, value):
        self._path = value
        try:
            self._open_path(value)
        except:
            print("Error opening path")
    def _open_path(self, path):
        """
        """
        # Read trn file
        self.lam, \
        self.wno, \
        self.trn_ray_0, \
        self.P_ray_0, \
        self.P_ray_str, \
        self.trn_gas_0, \
        self.P_gas_0, \
        self.P_gas_str, \
        self.trn_aero_0, \
        self.P_aero_0, \
        self.P_aero_str, \
        self.trn_all_0 \
            = read_radtrn(path, Nstr = self.Nstr)
        self._opened = True
    def plot_transmission(self, lammin = 0.2, lammax = 2.5, ax0 = None):
        """
        Plot the transmission at :math:`\\tau = 1`.
        Parameters
        ----------
        lammin : float
        lammax : float
        ax0 : matplotlib.Axis
        """
        # Plot it
        if ax0 is None:
            fig, ax = plt.subplots(figsize = (12,8))
            ax.set_ylabel(r"Transmission")
            ax.set_xlabel(r"Wavelength [$\mu$m]")
        else:
            ax = ax0
        mask = (self.lam > lammin) & (self.lam < lammax)
        #ax.plot(lam[mask],  trn_ray_0[mask])
        #ax.plot(lam[mask],  trn_gas_0[mask])
        #ax.plot(lam[mask],  trn_aero_0[mask])
        ax.plot(self.lam[mask],  self.trn_all_0[mask])
        if ax0 is None:
            return fig, ax
        else:
            return
    def plot_tau1_pressure(self, lammin = .2, lammax = 2.5, Pmin = 1e-2, Pmax = 1e5,
                            ax0 = None, show_rads = False):
        """
        Plot the pressure at :math:`\\tau = 1`.
        Parameters
        ----------
        lammin : float
        lammax : float
        Pmin : float
        Pmax : float
        ax0 : matplotlib.Axis
        show_rads : bool
        """
        if ax0 is None:
            fig, ax = plt.subplots(figsize = (12,8))
            ax.set_ylabel(r"$\tau = 1$ Pressure [Pa]")
            ax.set_xlabel(r"Wavelength [$\mu$m]")
            ax.set_yscale("log")
            ax.set_ylim(Pmax, Pmin)
        else:
            ax = ax0
        mask = (self.lam > lammin) & (self.lam < lammax)
        #ax.plot(self.lam[mask],  self.P_all_0[mask], c="k")
        ax.plot(self.lam[mask],  self.P_ray_0[mask])
        ax.plot(self.lam[mask],  self.P_gas_0[mask])
        ax.plot(self.lam[mask],  self.P_aero_0[mask])
        if show_rads:
            for i in range(P_gas_str.shape[0]):
                ax.plot(self.lam[mask],  self.P_ray_str[i,mask] , c = "C%i" %i, ls = "dotted")
                ax.plot(self.lam[mask],  self.P_gas_str[i,mask] , c = "C%i" %i, ls = "dotted")
                ax.plot(self.lam[mask],  self.P_aero_str[i,mask], c = "C%i" %i, ls = "dotted")
                #ax.plot(self.lam[mask],  self.P_all_str[i,mask], c = "C%i" %i, ls = "-")
                pass
        if ax0 is None:
            return fig, ax
        else:
            return


########################################################
################### MICROWAVE UTILS ####################
########################################################


def load_hitran(gfile='/home/linc/Documents/ASTROBIO/fixed_input/HITRAN2016.par',
                nlines=5701285):

    '''
    Loads the HITRAN line list.

    Outputs:
    --------
    molid : integer array
       HITRAN molecular ID number
    iso : integer array
       HITRAN isotopologue number
    nu : float array
       wavenumber of line [cm^-1]
    gint : float array
       intensity [cm^-2 / (molec^-1 cm^-2)]
    gam_air : float array
       air broadening parameter [cm^-1 atm^-1]
    gam_self : float array
       self broadening parameter [cm^-1 atm^-1]
    
    '''
    
    import numpy as np

    molid = np.zeros(nlines,dtype=np.int32)
    iso = np.zeros(nlines,dtype=np.int32)
    nu = np.zeros(nlines)
    gint = np.zeros(nlines)
    gam_air = np.zeros(nlines)
    gam_self = np.zeros(nlines)

    # read file format
    f = open(gfile)

    for i in range(nlines):
    
        s = f.readline()
        #print(len(s),s[:2])
    
        molid[i] = s[:2]
        try:
            iso[i] = s[2]
        except:
            iso[i] = -1
        nu[i]  = s[3:15]
        gint[i] = s[15:25]
        gam_air[i] = s[33:38]
        gam_self[i] = s[38:43]
    
    f.close()

    return molid,iso,nu,gint,gam_air,gam_self

def load_hitran_all(gfile='/home/linc/Documents/ASTROBIO/fixed_input/HITRAN2016.par',
                nlines=5701285):

    '''
    Loads the HITRAN line list.

    Outputs:
    --------
    molid : integer array
       HITRAN molecular ID number
    iso : integer array
       HITRAN isotopologue number
    nu : float array
       wavenumber of line [cm^-1]
    gint : float array
       intensity [cm^-2 / (molec^-1 cm^-2)]
    EA : float array
       Einstein A coefficient [s^-1]
    gam_air : float array
       air broadening parameter [cm^-1 atm^-1]
    gam_self : float array
       self broadening parameter [cm^-1 atm^-1]
    E0 : float array
       Lower state energy (E") [cm^-1]
    n_air : float array
       Temperature exponent for air broadening
    delta_air : float array
       Pressure shift induced by air at 1 atm [cm^-1 atm^-1]
    gu_q : str array
       Global Upper Quanta
    gl_q : str array
       Global Lower Quanta
    lu_q : str array
       Local Upper Quanta
    ll_q : str array
       Local Lower Quanta
    err : int array
       Error indices
    refs : int array
       References
    lmix : str array
       Line mixing flag
    g1 : float array
       Upper state degeneracy (g')
    g0 : float array
       Lower state degeneracy (g")
    
    '''
    
    import numpy as np

    molid = np.zeros(nlines,dtype=np.int32)
    iso = np.zeros(nlines,dtype=np.int32)
    nu = np.zeros(nlines)
    gint = np.zeros(nlines)
    EA = np.zeros(nlines)
    gam_air = np.zeros(nlines)
    gam_self = np.zeros(nlines)
    E0 = np.zeros(nlines)
    n_air= np.zeros(nlines)
    delta_air= np.zeros(nlines)
    #gu_q= np.zeros(nlines,dtype=object)
    #gl_q= np.zeros(nlines,dtype=object)
    #lu_q= np.zeros(nlines,dtype=object)
    #ll_q= np.zeros(nlines,dtype=object)
    gu_q= np.zeros(nlines,dtype='|S15')
    gl_q= np.zeros(nlines,dtype='|S15')
    lu_q= np.zeros(nlines,dtype='|S15')
    ll_q= np.zeros(nlines,dtype='|S15')
    #gu_q= np.chararray([nlines,15])
    #gl_q= np.chararray([nlines,15])
    #lu_q= np.chararray([nlines,15])
    #ll_q= np.chararray([nlines,15])
    err= np.zeros([6,nlines],dtype=np.int32)
    refs=np.zeros([6,nlines],dtype=np.int32)
    lmix= np.zeros(nlines,dtype=str)
    g1= np.zeros(nlines)
    g0= np.zeros(nlines)

    # read file format
    f = open(gfile)

    for i in range(nlines):
    
        s = f.readline()
        #print(len(s),s[:2])
    
        molid[i] = s[:2]
        try:
            iso[i] = s[2]
        except:
            iso[i] = -1
        nu[i]  = s[3:15]
        gint[i] = s[15:25]
        EA[i] = s[25:35]
        gam_air[i] = s[35:40]
        gam_self[i] = s[40:45]
        E0[i] = s[45:55]
        n_air[i] = s[55:59]
        delta_air[i] = s[59:67]
        gu_q[i] = s[67:82]
        gl_q[i] = s[82:97]
        lu_q[i] = s[97:112]
        ll_q[i] = s[112:127]
        idx = 127
        for j in range(6):
            err[j,i] = s[idx+j]
        idx = 133
        for j in range(6):
            refs[j,i] = s[idx+j*2:idx+j*2+2]
        lmix[i] = s[145]
        g1[i] = s[146:153]
        g0[i] = s[153:]
    
    f.close()

    return molid,iso,nu,gint,EA,gam_air,gam_self,E0,n_air,delta_air,gu_q,gl_q,lu_q,ll_q,err,refs,lmix,g1,g0

def load_hitran_broadening(gfile):

    '''
    Loads custom HITRAN data with different foreign broadening parameters.

    Outputs:
    --------
    molid : integer array
       HITRAN molecular ID number
    iso : integer array
       HITRAN isotopologue number
    nu : float array
       wavenumber of line [cm^-1]
    gint : float array
       intensity [cm^-2 / (molec^-1 cm^-2)]
    gam_for : float array
       foreign gas broadening parameter [cm^-1 atm^-1]
    n_for : float array
       Temperature exponent for foreign gas broadening
    delta_for : float array
       Pressure shift induced by air at 1 atm [cm^-1 atm^-1]
    err : int array
       Error indices    
    '''

    molid,iso,nu,gint,gam_for,n_for,delta_for,err1,err2,err3,err4,err5,ref1,ref2,ref3,ref4,ref5 = np.loadtxt(gfile,skiprows=1,unpack=True)

    err = np.array(np.vstack([err3,err4,err5]),dtype=np.int32)
    refs = np.array(np.vstack([ref3,ref4,ref5]),dtype=np.int32)

    return molid,iso,nu,gint,gam_for,n_for,delta_for,err,refs

def write_hitran(gfile,molid,iso,nu,gint,EA,gam_air,gam_self,E0,n_air,delta_air,gu_q,gl_q,lu_q,ll_q,err,refs,lmix,g1,g0,binary=False):

    '''
    Loads the HITRAN line list.

    Inputs:
    --------
    gfile : string
       path to output file
    molid : integer array
       HITRAN molecular ID number
    iso : integer array
       HITRAN isotopologue number
    nu : float array
       wavenumber of line [cm^-1]
    gint : float array
       intensity [cm^-2 / (molec^-1 cm^-2)]
    EA : float array
       Einstein A coefficient [s^-1]
    gam_air : float array
       air broadening parameter [cm^-1 atm^-1]
    gam_self : float array
       self broadening parameter [cm^-1 atm^-1]
    E0 : float array
       Lower state energy (E") [cm^-1]
    n_air : float array
       Temperature exponent for air broadening
    delta_air : float array
       Pressure shift induced by air at 1 atm [cm^-1 atm^-1]
    gu_q : str array
       Global Upper Quanta
    gl_q : str array
       Global Lower Quanta
    lu_q : str array
       Local Upper Quanta
    ll_q : str array
       Local Lower Quanta
    err : int array
       Error indices
    refs : int array
       References
    lmix : str array
       Line mixing flag
    g1 : float array
       Upper state degeneracy (g')
    g0 : float array
       Lower state degeneracy (g")
    
    '''
    
    import numpy as np

    # read file format
    f = open(gfile,'w')

    nlines = len(nu)

    for i in range(nlines):

        f.write('{:2d}{:1d}{:12.6F}{:10.3E}{:10.3E}'.format(molid[i],iso[i],nu[i],gint[i],EA[i]))
        f.write('{:5.4F}'.format(gam_air[i]).lstrip('0'))
        f.write('{:5.3F}{:10.4F}{:4.2F}'.format(gam_self[i],E0[i],n_air[i]))
        if(delta_air[i] < 0.):
            f.write('-'+'{:8.6F}'.format(-delta_air[i]).lstrip('0'))
        else:
            f.write('{:8.6F}'.format(delta_air[i]))
        f.write('{:15}{:15}{:15}{:15}{:1d}{:1d}{:1d}{:1d}{:1d}{:1d}'.format(gu_q[i].decode('UTF-8'),gl_q[i].decode('UTF-8'),lu_q[i].decode('UTF-8'),ll_q[i].decode('UTF-8'),err[0,i],err[1,i],err[2,i],err[3,i],err[4,i],err[5,i]))
        f.write('{:2d}{:2d}{:2d}{:2d}{:2d}{:2d}{:1}{:7.1F}{:7.1F}\n'.format(refs[0,i],refs[1,i],refs[2,i],refs[3,i],refs[4,i],refs[5,i],lmix[i],g1[i],g0[i]))
        
    f.close()


def calc_continuum(x,y,nfreq,ifreq,nfreq1=1,ifreq1=0,freqcut0=None):
    '''
    Calculates the continuum from a sample of the full spectrum and interpolate using spline functions. 
    This uses nfreq and ifreq to determine sampling.

    Inputs
    ------
    x : float array
       wavelength / frequency / energy etc vector
    y : float array
       flux / intensity etc vector
    nfreq : integer
       sampling frequency of points
    ifreq : integer
       initial offset of point sampling
    '''

    from scipy.interpolate import interp1d

    continuum_pts = x[ifreq::nfreq]
    splfunc = interp1d(continuum_pts,y[ifreq::nfreq],kind='cubic',bounds_error=False)
    
    if(freqcut0):
        
        c = 2.99792458e10
        
        freqcut1 = freqcut0 * 1e9 / c
        continuum_pts1 = x[ifreq1::nfreq1]
        splfunc1 = interp1d(continuum_pts1,y[ifreq1::nfreq1],kind='cubic',bounds_error=False)
        
        # Now interpolate to continuum
        continuum0 = splfunc(x)
        continuum1 = splfunc1(x)

        continuum = np.hstack([continuum0[x < freqcut1],continuum1[x >= freqcut1]])
    else:
        
        continuum = splfunc(x)

    return continuum

def calc_continuum_mask(x,y,nfreq,ifreq,masks=[],kind='cubic'):
    '''
    Calculates the continuum from a sample of the full spectrum and interpolate using spline functions. 
    This uses nfreq and ifreq to determine sampling.

    Inputs
    ------
    x : float array
       wavelength / frequency / energy etc vector
    y : float array
       flux / intensity etc vector
    nfreq : integer
       sampling frequency of points
    ifreq : integer
       initial offset of point sampling
    '''
    
    
    #                     if(hasattr(self.lblcutoff,"__len__")):


    from scipy.interpolate import interp1d

    mask = np.ones(len(x),dtype=np.bool)
    
    # Now determine data mask

    print('Masking data ranges:')
    for i in range(len(masks)):
        print(masks[i])
        mask[(x <= masks[i][1]) & (x >= masks[i][0])] = False
    

    continuum_pts = x[mask][ifreq::nfreq]
    splfunc = interp1d(continuum_pts,y[mask][ifreq::nfreq],kind=kind,bounds_error=False)
    
    continuum = splfunc(x)

    return continuum



def plot_flux_stems(x0,y,
                    nu,molid,gint,gases,gasnames,
                    yc=None,
                    nfreq=1,ifreq=0,freq0=0.,nfreq1=1,ifreq1=0,freqcut1=None,masks=[],
                    figsize=[20,10],height_ratios=[5,2],ylabel='',yclabel='',ycolor='k',yccolor='grey',
                    xlim=None,ylim=None,xunits='GHz',yglim=None,rvxlim=None,
                    plotcont=False,adjcont=False,plotrv=False,plotstep=False,plotbt=False,fig=None,axes=None,
                    ls='-',lw=1.5,plotstems=True,zorder=10):
    '''
    Generalized plotting w/ HITRAN lines. 
    Options to plot as flux or continuum-adjusted. Option to plot dots for continuum.
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from brightness_temp import brightness_temperature
    
    # Speed of light in cm/s
    c = 2.99792458e10

    if(plotrv and xunits=='GHz'):
        x = c / 100. *(freq0-x0)/freq0 * 1e-3 # km/s
    elif(plotrv):
        print('Error, RV not set up with non-GHz units!')
        return 0.,0.
    else:
        x = np.array(x0)


    if(fig is None):
        fig,axes = plt.subplots(2,figsize=figsize,sharex=True,gridspec_kw={'height_ratios': height_ratios})

    ax = axes[0]

    if(freqcut1 is not None):
        xc0 = x[ifreq::nfreq]
        xc1 = x[ifreq1::nfreq1]
        yc0 = yc[ifreq::nfreq][xc0 < freqcut1]
        yc1 = yc[ifreq1::nfreq1][xc1 >= freqcut1]

        ycplot = np.hstack([yc0,yc1])
        xcplot = np.hstack([xc0[xc0 < freqcut1],xc1[xc1 >= freqcut1]])
        #ycplot = np.hstack([yc[x<freqcut1][ifreq::nfreq],yc[x>=freqcut1][ifreq1::nfreq1]])
        #xcplot = np.hstack([x[x<freqcut1][ifreq::nfreq],x[x>=freqcut1][ifreq1::nfreq1]])
    elif(len(masks) > 0):
        
        masklines = np.ones(len(x),dtype=np.bool)
        
        print('Masking line ranges:')
        for i in range(len(masks)):
            print(masks[i])
            masklines[(x <= masks[i][1]) & (x >= masks[i][0])] = False
            
        xcplot = x[masklines][ifreq::nfreq]
        ycplot = yc[masklines][ifreq::nfreq]
        
    elif(yc is not None):
        ycplot = yc[ifreq::nfreq]
        xcplot = x[ifreq::nfreq]


    if(plotbt):
        adjcont = False
        plotstep = False
        plotcont = False

        if(xunits=='GHz'):
            bt = brightness_temperature(c/x*1e-5,y)
        elif(xunits=='wn'):
            bt = brightness_temperature(1e4/x,y)
        else:
            bt = brightness_temperature(1e4/x,y)

        ax.plot(x,bt,color=ycolor,label=ylabel,ls=ls,lw=lw,zorder=zorder)

    else:
        if(adjcont):
            if(plotstep):
                ax.step(x,y/yc-1,color=ycolor,label=ylabel,where='mid',ls=ls,lw=lw,zorder=zorder)
            else:
                ax.plot(x,y/yc-1,color=ycolor,label=ylabel,ls=ls,lw=lw,zorder=zorder)
        else:
            if(plotstep):
                ax.step(x,y,color=ycolor,label=ylabel,where='mid',ls=ls,lw=lw,zorder=zorder)
            else:
                ax.plot(x,y,color=ycolor,label=ylabel,ls=ls,lw=lw,zorder=zorder)
        
    if(plotcont):

        if(adjcont):
            ax.plot(x,yc/yc-1.,color='grey',label=yclabel,zorder=zorder-1,ls=ls,lw=lw)  
            ax.plot(xcplot,ycplot/ycplot-1.,color=yccolor,label=yclabel,zorder=zorder-1,lw=0.,marker='o',ls=ls)
        else:
            ax.plot(x,yc,color='grey',label=yclabel,zorder=zorder-1,ls=ls,lw=lw)
            ax.plot(xcplot,ycplot,color=yccolor,label=yclabel,zorder=zorder-1,lw=0.,marker='o',ls=ls)

    if(adjcont):
        ax.axhline(0.,color='grey',lw=1.)

    if(plotbt):
        ax.set_ylabel(r'Brightness Temperature [K]')
    elif(adjcont):
        ax.set_ylabel(r'Line/Continuum--1')
    else:
        ax.set_ylabel(r'Emergent Specific Flux [W/m$^2$/$\upmu$m]')
        
    ax.minorticks_on()

    ax.legend()

    if(xlim is None):
        xlim = np.zeros(2)
        xlim[0] = np.min(x)
        xlim[1] = np.max(x)

    if(ylim is not None):
        ax.set_ylim(ylim)


    print(xlim)
    
    if(rvxlim is None):
        ax.set_xlim(xlim)
    else:
        ax.set_xlim(rvxlim)
        
        
    if(adjcont is False or plotbt is True):
        ax.set_yscale('log')

    if(plotrv):
        ax.set_title('{:7.3f} GHz center'.format(freq0),fontsize=14)
        

    ### LINE STEMS ###

    if(plotstems):
        ax = axes[1]


        # nu_run is in units of cm^-1 (like HITRAN)

        if(plotrv): # must get limits 
            if(xunits == 'GHz'):
                nu_run = np.zeros(2)
                nu_run[0] = np.min(x0)
                nu_run[1] = np.max(x0)

    
            else:
                print('error nu_run!')
                return fig,axes
        else:
            if(xunits == 'GHz'):
                nu_run = xlim / c * 1e9
            elif(xunits == 'wn'):
                nu_run = np.array(xlim)
            else:
                print('Error, xunits!')
                return fig,axes

        print(nu_run)

        wl_run = c / 100. / (nu_run* 1e9) * 1e6
        wn_run = 1e4/wl_run

        if(plotstems):
            print(nu_run)
            print(freq0)
            for i in range(len(gases)):
                print('gas ',i,gases[i],gasnames[i])
                try:
                    if(plotrv):
                        rv1 = c /100. *(freq0-nu[molid == gases[i]][(nu[molid == gases[i]] < wn_run[1]) & (nu[molid == gases[i]] > wn_run[0])]*c*1e-9)/freq0 * 1e-3
                        print(rv1)
                        ax.stem(rv1,
                                gint[molid == gases[i]][(nu[molid == gases[i]] < wn_run[1]) & (nu[molid == gases[i]] > wn_run[0])],
                                use_line_collection=True,markerfmt='C'+str(i) + '.',label=r'\ce{'+gasnames[i]+r'}',basefmt='C'+str(i),linefmt='C'+str(i))
                    else:
                        ax.stem(nu[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])]*c*1e-9,
                                gint[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])],
                                use_line_collection=True,markerfmt='C'+str(i) + '.',label=r'\ce{'+gasnames[i]+r'}',basefmt='C'+str(i),linefmt='C'+str(i))
                    print(gases[i],nu[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])]*c*1e-9)
                except:
                    print('gas',gases[i],'has no lines')
                    continue
    
        ax.set_yscale('log')
    
        if(yglim):
            ax.set_ylim(yglim)


        ax.legend()

        ax.minorticks_on()

        if(rvxlim):
            ax.set_xlim(rvxlim)
        else:
            ax.set_xlim(xlim)
        

    if(plotrv):
        ax.set_xlabel('velocity [km/s]')
    elif(xunits=='GHz'):
        ax.set_xlabel('Frequency [GHz]')
    elif(xunits=='wn'):
        ax.set_xlabel(r'Wavenumber [cm$^{-1}$]')
    else:
        ax.set_xlabel(r'Wavelength [$\upmu$m]')

    plt.subplots_adjust(hspace=0.)

    return fig,axes

def plot_flux_stems_old(x0,y,
                    nu,molid,gint,gases,gasnames,
                    yc=None,
                    nfreq=1,ifreq=0,freq0=0.,nfreq1=1,ifreq1=0,freqcut1=None,
                    figsize=[20,10],height_ratios=[5,2],ylabel='',yclabel='',ycolor='k',yccolor='grey',
                    xlim=None,ylim=None,xunits='GHz',yglim=None,rvxlim=None,
                    plotcont=False,adjcont=False,plotrv=False,plotstep=False,plotbt=False,fig=None,axes=None):
    '''
    Generalized plotting w/ HITRAN lines. 
    Options to plot as flux or continuum-adjusted. Option to plot dots for continuum.
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from brightness_temp import brightness_temperature
    
    # Speed of light in cm/s
    c = 2.99792458e10

    if(plotrv and xunits=='GHz'):
        x = c / 100. *(freq0-x0)/freq0 * 1e-3 # km/s
    elif(plotrv):
        print('Error, RV not set up with non-GHz units!')
        return 0.,0.
    else:
        x = np.array(x0)


    if(fig is None):
        fig,axes = plt.subplots(2,figsize=figsize,sharex=True,gridspec_kw={'height_ratios': height_ratios})

    ax = axes[0]

    if(freqcut1 is not None):
        xc0 = x[ifreq::nfreq]
        xc1 = x[ifreq1::nfreq1]
        yc0 = yc[ifreq::nfreq][xc0 < freqcut1]
        yc1 = yc[ifreq1::nfreq1][xc1 >= freqcut1]

        ycplot = np.hstack([yc0,yc1])
        xcplot = np.hstack([xc0[xc0 < freqcut1],xc1[xc1 >= freqcut1]])
        #ycplot = np.hstack([yc[x<freqcut1][ifreq::nfreq],yc[x>=freqcut1][ifreq1::nfreq1]])
        #xcplot = np.hstack([x[x<freqcut1][ifreq::nfreq],x[x>=freqcut1][ifreq1::nfreq1]])
    elif(yc is not None):
        ycplot = yc[ifreq::nfreq]
        xcplot = x[ifreq::nfreq]


    if(plotbt):
        adjcont = False
        plotstep = False
        plotcont = False

        if(xunits=='GHz'):
            bt = brightness_temperature(c/x*1e-5,y)
        elif(xunits=='wn'):
            bt = brightness_temperature(1e4/x,y)
        else:
            bt = brightness_temperature(1e4/x,y)

        ax.plot(x,bt,color=ycolor,label=ylabel,zorder=10)

    else:
        if(adjcont):
            if(plotstep):
                ax.step(x,y/yc-1,color=ycolor,label=ylabel,zorder=10,where='mid')
            else:
                ax.plot(x,y/yc-1,color=ycolor,label=ylabel,zorder=10)
        else:
            if(plotstep):
                ax.step(x,y,color=ycolor,label=ylabel,zorder=10,where='mid')
            else:
                ax.plot(x,y,color=ycolor,label=ylabel,zorder=10)
        
    if(plotcont):

        if(adjcont):
            ax.plot(x,yc/yc-1.,color='grey',label=yclabel,zorder=8)    
            ax.plot(xcplot,ycplot/ycplot-1.,color=yccolor,label=yclabel,zorder=8,lw=0.,marker='o')
        else:
            ax.plot(x,yc,color='grey',label=yclabel,zorder=8)
            ax.plot(xcplot,ycplot,color=yccolor,label=yclabel,zorder=8,lw=0.,marker='o')

    if(adjcont):
        ax.axhline(0.,color='grey',lw=1.)

    if(plotbt):
        ax.set_ylabel(r'Brightness Temperature [K]')
    elif(adjcont):
        ax.set_ylabel(r'Line/Continuum--1')
    else:
        ax.set_ylabel(r'Emergent Specific Flux [W/m$^2$/$\upmu$m]')
        
    ax.minorticks_on()

    ax.legend()

    if(xlim is None):
        xlim = np.zeros(2)
        xlim[0] = np.min(x)
        xlim[1] = np.max(x)

    if(ylim is not None):
        ax.set_ylim(ylim)


    print(xlim)
    
    if(rvxlim is None):
        ax.set_xlim(xlim)
    else:
        ax.set_xlim(rvxlim)
        
        
    if(adjcont is False):
        ax.set_yscale('log')

    if(plotrv):
        ax.set_title('{:7.3f} GHz center'.format(freq0),fontsize=14)
        

    ### LINE STEMS ###
    
    ax = axes[1]


    # nu_run is in units of cm^-1 (like HITRAN)

    if(plotrv): # must get limits 
        if(xunits == 'GHz'):
            nu_run = np.zeros(2)
            nu_run[0] = np.min(x0)
            nu_run[1] = np.max(x0)

    
        else:
            print('error nu_run!')
            return fig,axes
    else:
        if(xunits == 'GHz'):
            nu_run = xlim / c * 1e9
        elif(xunits == 'wn'):
            nu_run = np.array(xlim)
        else:
            print('Error, xunits!')
            return fig,axes

    print(nu_run)

    wl_run = c / 100. / (nu_run* 1e9) * 1e6
    wn_run = 1e4/wl_run

    for i in range(len(gases)):
        print('gas ',i,gases[i],gasnames[i])
        try:
            if(plotrv):
                rv1 = c /100. *(freq0-nu[molid == gases[i]][(nu[molid == gases[i]] < wn_run[1]) & (nu[molid == gases[i]] > wn_run[0])]*c*1e-9)/freq0 * 1e-3
                ax.stem(rv1,
                        gint[molid == gases[i]][(nu[molid == gases[i]] < wn_run[1]) & (nu[molid == gases[i]] > wn_run[0])],
                        use_line_collection=True,markerfmt='C'+str(i) + '.',label=r'\ce{'+gasnames[i]+r'}',basefmt='C'+str(i),linefmt='C'+str(i))
            else:
                ax.stem(nu[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])]*c*1e-9,
                        gint[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])],
                        use_line_collection=True,markerfmt='C'+str(i) + '.',label=r'\ce{'+gasnames[i]+r'}',basefmt='C'+str(i),linefmt='C'+str(i))
            print(gases[i],nu[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])]*c*1e-9)
        except:
            print('gas',gases[i],'has no lines')
            continue
    
    ax.set_yscale('log')
    
    if(yglim):
        ax.set_ylim(yglim)

    ax.legend()

    ax.minorticks_on()

    if(rvxlim):
        ax.set_xlim(rvxlim)
    else:
        ax.set_xlim(xlim)
        

    if(plotrv):
        ax.set_xlabel('velocity [km/s]')
    elif(xunits=='GHz'):
        ax.set_xlabel('Frequency [GHz]')
    elif(xunits=='wn'):
        ax.set_xlabel(r'Wavenumber [cm$^{-1}$]')
    else:
        ax.set_xlabel(r'Wavelength [$\upmu$m]')

    plt.subplots_adjust(hspace=0.)

    return fig,axes


def plot_flux_stems2(x0,y,
                    nu,molid,gint,gases,gasnames,
                     yc=None,freq0=0.,
                     figsize=[20,10],height_ratios=[5,2],width_ratios=[2,5],ylabel='',yclabel='',ycolor='k',yccolor='grey',
                     xlim0=None,ylim=None,xunits='GHz',yglim=None,rvxlim=None,rvylim=None,
                     plotstep=False):
    '''
    Specialized plotting for paper, w/ HITRAN lines. 

    Inputs
    ------
    x0 : Frequency (GHz)

    freq1 : location to change to different baseline sampling

    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator,MultipleLocator,LogLocator
    
    # Speed of light in cm/s
    c = 2.99792458e10

    # calculate RV
    x = c / 100. *(freq0-x0)/freq0 * 1e-3 # km/s
    
    fig,axes = plt.subplots(2,2,figsize=figsize,gridspec_kw={'height_ratios': height_ratios, 'width_ratios': width_ratios})

    ax = axes[0,0]

    if(plotstep):
        ax.step(x,y/yc-1,color=ycolor,label=ylabel,zorder=10,where='mid')
    else:
        ax.plot(x,y/yc-1,color=ycolor,label=ylabel,zorder=10)
        
    ax.axhline(0.,color='grey',lw=1.)

    ax.set_ylabel(r'Line/Continuum--1')
    plt.setp(ax.get_xticklabels(), visible=False)

    ax.legend()

    if(xlim0 is None):
        xlim = np.zeros(2)
        xlim[0] = np.min(x)
        xlim[1] = np.max(x)
    else:
        xlim = np.array(xlim0)

    if(ylim is not None):
        ax.set_ylim(ylim)

   
    if(rvxlim is None):
        ax.set_xlim(xlim)
    else:
        ax.set_xlim(rvxlim)

    if(rvylim is None):
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(rvylim)

    ax.set_title('{:7.3f} GHz center'.format(freq0),fontsize=14)
        

    ### LINE STEMS ###
    
    ax = axes[1,0]


    # nu_run is in units of cm^-1 (like HITRAN)
    
    if(xunits == 'GHz'):
        nu_run = np.zeros(2)
        nu_run[0] = np.min(x0)
        nu_run[1] = np.max(x0)

    
    else:
        print('error nu_run!')
        return fig,axes

    print(nu_run)

    wl_run = c / 100. / (nu_run* 1e9) * 1e6
    wn_run = 1e4/wl_run

    if(rvxlim):
        xlim_stems = rvxlim
    else:
        xlim_stems = xlim

    for i in range(len(gases)):
        print('gas ',i,gases[i],gasnames[i])
        maskgas = (molid == gases[i])
        maskstems = [(nu[molid == gases[i]] < wn_run[1]) & (nu[molid == gases[i]] > wn_run[0])]
        rv1 = c /100. *(freq0-nu[maskgas][maskstems]*c*1e-9)/freq0 * 1e-3
        maskstems1 = (rv1 <= xlim_stems[1]) & (rv1 >= xlim_stems[0])
        print('rv:',rv1[maskstems1])
        print('gint:',gint[maskgas][maskstems][maskstems1])
        print('nu:',nu[maskgas][maskstems][maskstems1]*c*1e-9)
        try:
            ax.stem(rv1[maskstems1],gint[maskgas][maskstems][maskstems1],
                    use_line_collection=True,markerfmt='C'+str(i) + '.',label=r'\ce{'+gasnames[i]+r'}',basefmt='C'+str(i),linefmt='C'+str(i))
            print(gases[i],nu[maskstems][maskstems1]*c*1e-9)
        except:
            print('gas code',gases[i],'has no lines')
            continue
    
    ax.set_yscale('log')

    if(yglim):
        ax.set_ylim(yglim)

    ax.legend()

    ax.set_xlabel('Velocity [km/s]')

    if(rvxlim):
        ax.set_xlim(rvxlim)
        ax.get_xaxis().set_major_locator(MultipleLocator(20.))
        ax.get_xaxis().set_minor_locator(AutoMinorLocator(2))
    else:
        ax.set_xlim(xlim)
        
    #######################
    ##### LARGER PLOT #####
    #######################

    x = np.array(x0)
    
    ax = axes[0,1]

    ax.plot(x,y/yc-1,color=ycolor,label=ylabel,zorder=10)
        
    ax.axhline(0.,color='grey',lw=1.)

    ax.set_ylabel(r'Line/Continuum--1')
    ax.minorticks_on()

    plt.setp(ax.get_xticklabels(), visible=False)

    ax.legend()

    if(xlim0 is None):
        xlim = np.zeros(2)
        xlim[0] = np.min(x)
        xlim[1] = np.max(x)
    else:
        xlim = np.array(xlim0)

    if(ylim is not None):
        ax.set_ylim(ylim)

    ax.set_xlim(xlim)


    ### LINE STEMS ###
    
    ax = axes[1,1]


    # nu_run is in units of cm^-1 (like HITRAN)

    if(xunits == 'GHz'):
        nu_run = xlim / c * 1e9
    elif(xunits == 'wn'):
        nu_run = np.array(xlim)
    else:
        print('Error, xunits!')
        return fig,axes

    print(nu_run)

    wl_run = c / 100. / (nu_run* 1e9) * 1e6
    wn_run = 1e4/wl_run

    for i in range(len(gases)):
        print('gas ',i,gases[i],gasnames[i])
        try:
            ax.stem(nu[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])]*c*1e-9,
                gint[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])],
                use_line_collection=True,markerfmt='C'+str(i) + '.',label=r'\ce{'+gasnames[i]+r'}',basefmt='C'+str(i),linefmt='C'+str(i))
            print(gases[i],nu[molid == gases[i]][(nu[molid == gases[i]] < nu_run[1]) & (nu[molid == gases[i]] > nu_run[0])]*c*1e-9)
        except:
            print('gas',gases[i],'has no lines')
            continue
    
    ax.set_yscale('log')
    
    if(yglim):
        ax.set_ylim(yglim)

    ax.legend(ncol=4,fontsize=12)

    ax.minorticks_on()

    ax.set_xlim(xlim)

    ax.set_xlabel('Frequency [GHz]')
    
    plt.subplots_adjust(hspace=0.,wspace=0.25)

    return fig,axes
