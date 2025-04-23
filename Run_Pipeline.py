
from Pipeline import *

'''
pipelineobj = VPLModelingPipeline('RunNumber189', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/O2Sens4TestH2O/RunNumber189/PhotochemInputs/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

pipelineobj.OutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/O2Sens4TestH2O/RunNumber189/'
trn = '14_Subtry3'
pipelineobj.photochem_InputsDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/O2Sens4TestH2O/RunNumber189/PhotochemInputs/'
oldptz = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/O2Sens4TestH2O/RunNumber189/atmos/PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist'
subprocess.run('cp '+pipelineobj.photochem_InputsDir+'in.dist '+pipelineobj.photochem_InputsDir+'save_from_before_thickatmclimate_in.dist', shell=True)
dat = pipelineobj.get_final_climate_output_temp_profile(trynum=trn)
pipelineobj.update_indist_T_EDD(oldptz, dat)
'''

#/gscratch/vsm/gialluca/VPLModelingTools_Dev/VeffTestDepos/RunNumber6/PhotochemInputs

pipelineobj = VPLModelingPipeline('Template100mbar', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/InputDir_100mbar/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

pipelineobj.run_spectra = True
pipelineobj.adjust_atmospheric_pressure = False

pipelineobj.run_automatic()


### JUST RUN SMART SPECTRA:
'''

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
pipelineobj.MMW = 31.63843524999342

pipelineobj.set_climate_settings()
pipelineobj.set_smart_settings()


#pipelineobj.run_automatic()

pipelineobj.surface_temp = 3.6214E+02
pipelineobj.surface_temp_dayside = 4.1843E+02
pipelineobj.surface_temp_nightside = 2.5046E+02

# First, get the final profile from the last climate run to use in the restart
climate_profile = pipelineobj.get_final_2column_climate_output_temp_profile(trynum=1)

# Update temperature in the PT profile and update the surface temperature
pipelineobj.replace_PT_tempcol(climate_profile['Dayside']['T[K]'], whichcolumn='dayside')
pipelineobj.replace_PT_tempcol(climate_profile['Nightside']['T[K]'], whichcolumn='nightside')

pipelineobj.make_lblabc_runscripts()

# Now Run LBLABC for all the gases of interest
for gas in pipelineobj.molecule_dict['Gas_names']:
    pipelineobj.run_lblabc_1instance(pipelineobj.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+pipelineobj.casename+'.script', gas)
    if pipelineobj.verbose == True:
        print('LBLABC run for '+gas+' complete, LBLABC iteration '+str(pipelineobj.num_lblabc_runs+1))
        #ftestingoutput.write('LBLABC run for '+gas+' complete, LBLABC iteration '+str(self.num_lblabc_runs+1)+'\n')
pipelineobj.num_lblabc_runs += 1

pipelineobj.make_smart_runscript()

pipelineobj.run_smart_1instance(pipelineobj.SMART_RunScriptDir+'RunSMART_'+pipelineobj.casename+'.run')

pipelineobj.make_lblabc_runscripts(whichcol='dayside')

# Now Run LBLABC for all the gases of interest
for gas in pipelineobj.molecule_dict['Gas_names']:
    pipelineobj.run_lblabc_1instance(pipelineobj.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+pipelineobj.casename+'.script', gas)
    if pipelineobj.verbose == True:
        print('LBLABC run for '+gas+' complete, LBLABC iteration '+str(pipelineobj.num_lblabc_runs+1))
        #ftestingoutput.write('LBLABC run for '+gas+' complete, LBLABC iteration '+str(self.num_lblabc_runs+1)+'\n')
pipelineobj.num_lblabc_runs += 1

pipelineobj.make_smart_runscript(whichcol='dayside')

pipelineobj.run_smart_1instance(pipelineobj.SMART_RunScriptDir+'RunSMART_dayside_'+pipelineobj.casename+'.run', whichcol='dayside')

pipelineobj.make_lblabc_runscripts(whichcol='nightside')

# Now Run LBLABC for all the gases of interest
for gas in pipelineobj.molecule_dict['Gas_names']:
    pipelineobj.run_lblabc_1instance(pipelineobj.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+pipelineobj.casename+'.script', gas)
    if pipelineobj.verbose == True:
        print('LBLABC run for '+gas+' complete, LBLABC iteration '+str(pipelineobj.num_lblabc_runs+1))
        #ftestingoutput.write('LBLABC run for '+gas+' complete, LBLABC iteration '+str(self.num_lblabc_runs+1)+'\n')
pipelineobj.num_lblabc_runs += 1

pipelineobj.make_smart_runscript(whichcol='nightside')

pipelineobj.run_smart_1instance(pipelineobj.SMART_RunScriptDir+'RunSMART_nightside_'+pipelineobj.casename+'.run', whichcol='nightside')

'''