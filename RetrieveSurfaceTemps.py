import json

def get_surface_temp_one_sim(smart_runscript):

    f = open(smart_runscript, 'r')
    lines = f.readlines()

    for l in lines:
        if len(l.split('Surface temperature')) > 1:
            temp = float(l.split()[0])

    return temp

def get_surfaceTemps(planet):

    if planet == 'T1b' or planet == 'TRAPPIST-1b':
        finaldb = '../VPLModelingSupportScripts/T1b_FinalSpectra_Database.json'
    elif planet == 'T1c' or planet == 'TRAPPIST-1c':
        finaldb = '../VPLModelingSupportScripts/T1c_FinalSpectra_Database.json'
    elif planet == 'T1d' or planet == 'TRAPPIST-1d':
        finaldb = '../VPLModelingSupportScripts/T1d_FinalSpectra_Database.json'
    elif planet == 'T1e' or planet == 'TRAPPIST-1e':
        finaldb = '../VPLModelingSupportScripts/T1e_FinalSpectra_Database.json'
    elif planet == 'T1f' or planet == 'TRAPPIST-1f':
        finaldb = '../VPLModelingSupportScripts/T1f_FinalSpectra_Database.json'
    elif planet == 'T1g' or planet == 'TRAPPIST-1g':
        finaldb = '../VPLModelingSupportScripts/T1g_FinalSpectra_Database.json'
    elif planet == 'T1h' or planet == 'TRAPPIST-1h':
        finaldb = '../VPLModelingSupportScripts/T1h_FinalSpectra_Database.json'

    # Read in the final DB
    f = open(finaldb, 'r')
    fj = json.load(f)
    fj = json.loads(fj)
    fh2o = fj['H2O-O2']
    fco2 = fj['CO2']
    if planet not in ['T1f', 'T1g', 'T1h']:
        fso2h2o = fj['SO2-H2O']
    fso2co2 = fj['SO2-CO2']
    f.close()

    d = {}
    d['H2O-O2'] = {}
    d['CO2'] = {}
    d['SO2-H2O'] = {}
    d['SO2-CO2'] = {}

    if planet not in ['T1f', 'T1g', 'T1h']:
        for atm in fh2o['AtmIDs']:

            path = fh2o['Atm'+str(atm)]['OriginalPath']
            splpath = path.split('/')
            runnum = splpath[len(splpath)-2]
            smartscript = 'RunSMART_'+runnum+'.run'
            temp = get_surface_temp_one_sim(path+smartscript)

            d['H2O-O2']['Atm'+str(atm)] = temp

            if planet in ['T1b', 'T1c', 'T1d']:
                daysidesmartscript = 'RunSMART_dayside_'+runnum+'.run'
                temp = get_surface_temp_one_sim(path+daysidesmartscript)
                d['H2O-O2']['Atm'+str(atm)+'_Day'] = temp

                nightsidesmartscript = 'RunSMART_nightside_'+runnum+'.run'
                temp = get_surface_temp_one_sim(path+nightsidesmartscript)
                d['H2O-O2']['Atm'+str(atm)+'_Night'] = temp

    for atm in fco2['AtmIDs']:

        path = fco2['Atm'+str(atm)]['OriginalPath']
        splpath = path.split('/')
        runnum = splpath[len(splpath)-2]
        smartscript = 'RunSMART_'+runnum+'.run'
        temp = get_surface_temp_one_sim(path+smartscript)

        d['CO2']['Atm'+str(atm)] = temp

        if planet in ['T1b', 'T1c', 'T1d']:
                daysidesmartscript = 'RunSMART_dayside_'+runnum+'.run'
                temp = get_surface_temp_one_sim(path+daysidesmartscript)
                d['CO2']['Atm'+str(atm)+'_Day'] = temp

                nightsidesmartscript = 'RunSMART_nightside_'+runnum+'.run'
                temp = get_surface_temp_one_sim(path+nightsidesmartscript)
                d['CO2']['Atm'+str(atm)+'_Night'] = temp

    if planet not in ['T1f', 'T1g', 'T1h']:
        for atm in fso2h2o['AtmIDs']:

            path = fso2h2o['Atm'+str(atm)]['OriginalPath']
            splpath = path.split('/')
            runnum = splpath[len(splpath)-2]
            smartscript = 'RunSMART_'+runnum+'.run'
            temp = get_surface_temp_one_sim(path+smartscript)

            d['SO2-H2O']['Atm'+str(atm)] = temp

            if planet in ['T1b', 'T1c', 'T1d']:
                daysidesmartscript = 'RunSMART_dayside_'+runnum+'.run'
                temp = get_surface_temp_one_sim(path+daysidesmartscript)
                d['SO2-H2O']['Atm'+str(atm)+'_Day'] = temp

                nightsidesmartscript = 'RunSMART_nightside_'+runnum+'.run'
                temp = get_surface_temp_one_sim(path+nightsidesmartscript)
                d['SO2-H2O']['Atm'+str(atm)+'_Night'] = temp

    for atm in fso2co2['AtmIDs']:

        path = fso2co2['Atm'+str(atm)]['OriginalPath']
        splpath = path.split('/')
        runnum = splpath[len(splpath)-2]
        smartscript = 'RunSMART_'+runnum+'.run'
        temp = get_surface_temp_one_sim(path+smartscript)

        d['SO2-CO2']['Atm'+str(atm)] = temp

        if planet in ['T1b', 'T1c', 'T1d']:
                daysidesmartscript = 'RunSMART_dayside_'+runnum+'.run'
                temp = get_surface_temp_one_sim(path+daysidesmartscript)
                d['SO2-CO2']['Atm'+str(atm)+'_Day'] = temp

                nightsidesmartscript = 'RunSMART_nightside_'+runnum+'.run'
                temp = get_surface_temp_one_sim(path+nightsidesmartscript)
                d['SO2-CO2']['Atm'+str(atm)+'_Night'] = temp

    
    f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/'+planet+'_SurfaceTemps.json', 'w')
    dh = json.dumps(d)
    json.dump(dh, f)
    f.close()