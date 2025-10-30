from astropy.io import ascii 
import json

# Compile all the tested input options for outgassing, escape, etc to later make a histogram 

def compile_for_hist_inputs(planet, sweepdirs, atmtype):

    h2ooutg = []
    co2outg = []
    otoa = []
    o2toa = []
    o3dep = []
    h2o2dep = []
    codep = []
    co2toa = []

    h2ooutg_stable = []
    co2outg_stable = []
    otoa_stable = []
    o2toa_stable = []
    o3dep_stable = []
    h2o2dep_stable = []
    codep_stable = []
    co2toa_stable = []

    for sd, sweep_dir in enumerate(sweepdirs):
        partab = ascii.read(sweep_dir+'ParameterSweep_RunStats_failedrun.dat', format='fixed_width', delimiter=' ')
            
        for i, status in enumerate(partab['FinalState']):
            h2ooutg.append(float(partab['H2O_OutgassRate'][i]))
            otoa.append(float(partab['O_EscapeRate'][i]))
            o2toa.append(float(partab['O2_EscapeRate'][i]))
            o3dep.append(float(partab['O3_EscapeRate'][i]))
            h2o2dep.append(float(partab['H2O2_EscapeRate'][i]))

            if atmtype[i] == 'CO2':
                co2outg.append(float(partab['CO2_OutgassRate'][i]))
                co2toa.append(float(partab['CO2_EscapeRate'][i]))
                codep.append(float(partab['CO_EscapeRate'][i]))

            if status == 'Converged':
                h2ooutg_stable.append(float(partab['H2O_OutgassRate'][i]))
                otoa_stable.append(float(partab['O_EscapeRate'][i]))
                o2toa_stable.append(float(partab['O2_EscapeRate'][i]))
                o3dep_stable.append(float(partab['O3_EscapeRate'][i]))
                h2o2dep_stable.append(float(partab['H2O2_EscapeRate'][i]))

                if atmtype[i] == 'CO2':
                    co2outg_stable.append(float(partab['CO2_OutgassRate'][i]))
                    co2toa_stable.append(float(partab['CO2_EscapeRate'][i]))
                    codep_stable.append(float(partab['CO_EscapeRate'][i]))
            
        
    d = {}
    d['h2ooutg'] = h2ooutg
    d['otoa'] = otoa
    d['o2toa'] = o2toa
    d['o3dep'] = o3dep
    d['h2o2dep'] = h2o2dep
    d['co2outg'] = co2outg
    d['co2toa'] = co2toa
    d['codep'] = codep

    d['h2ooutg_stable'] = h2ooutg_stable
    d['otoa_stable'] = otoa_stable
    d['o2toa_stable'] = o2toa_stable
    d['o3dep_stable'] = o3dep_stable
    d['h2o2dep_stable'] = h2o2dep_stable
    d['co2outg_stable'] = co2outg_stable
    d['co2toa_stable'] = co2toa_stable
    d['codep_stable'] = codep_stable

    f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/TestedInputsHistogramData_'+planet+'.json', 'w')
    dh = json.dumps(d)
    json.dump(dh, f)
    f.close()
