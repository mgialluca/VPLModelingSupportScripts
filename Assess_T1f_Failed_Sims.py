import os

f_sweep_dirs = ['../Fco2/', '../Fco2lo/', '../Fco2wL/', '../Fco2wL2/']
g_sweep_dirs = ['../Gco2/', '../Gco2lo/', '../Gco2wL/', '../Gco2wL2/']
h_sweep_dirs = ['../Hco2/', '../Hco2lo/', '../Hco2wL/']

planets = ['T1f', 'T1g', 'T1h']

for planet, swpdrs in zip(planets, [f_sweep_dirs, g_sweep_dirs, h_sweep_dirs]):
    passed = []
    timedout = []
    failed_in_finding_surfP = []
    failed_in_photochem_prior_to_surfP = []
    failed_in_photochem_inner_surfP_loop = []
    failed_in_climate = []
    failed_in_twocol_climate = []
    unaccounted_for = []
    failed_with_atleast_1_clim_run = []
    
    for swpdr in swpdrs:
        for ds, sds, fis in os.walk(swpdr):
            break

        for sd in sds:
            path = swpdr+sd+'/'
            for ds2, sds2, fis2 in os.walk(path):
                break

            if 'FINAL_out.out' in fis2:
                passed.append(path)

            else:
            #elif 'FINAL_out_FAILED.out' in fis2:

                # investigate here
                f_investigate = open(path+sd+'_SavingInfoOut.txt', 'r')
                lines = f_investigate.readlines()
                reason_found = False
                for l in reversed(range(len(lines))):
                    if (len(lines[l].split('Max iterations reached')) > 1 and len(lines[l].split('couldnt find new pressure, ending run')) > 1) or (len(lines[l].split('Max Iterations reached')) > 1 and len(lines[l].split('couldnt find new pressure, ending run')) > 1):
                        if len(lines[l-3].split('New Pressure:')) > 1:
                            failed_in_finding_surfP.append(path)
                            if 'vpl_climate_output_'+sd+'.run' in fis2:
                                failed_with_atleast_1_clim_run.append(path)
                            break
                        else:
                            if 'vpl_climate_output_'+sd+'.run' in fis2:
                                failed_with_atleast_1_clim_run.append(path)
                            failed_in_photochem_prior_to_surfP.append(path)
                            break

                    elif len(lines[l].split('Max iterations reached')) > 1 and len(lines[l].split('using new pressure without finding inner convergence')) > 1:
                        if 'vpl_climate_output_'+sd+'.run' in fis2:
                            failed_with_atleast_1_clim_run.append(path)
                        failed_in_photochem_inner_surfP_loop.append(path)
                        break

                    elif len(lines[l].split('First Climate run')) > 1:
                        timedout.append(path)
                        failed_in_climate.append(path)
                        failed_with_atleast_1_clim_run.append(path)
                        break
                    
                    elif len(lines[l].split('Climate convergence NOT')) > 1:
                        failed_in_climate.append(path)
                        failed_with_atleast_1_clim_run.append(path)
                        break

                    if l == len(lines)-25:
                        if len(lines[len(lines)-1].split('runspec: True, include 2 col: False')) <=1:
                            timedout.append(path)
                        else:
                            unaccounted_for.append(path)
                        break

    
    print(planet+' RESULTS:')
    print('Passed Sims: '+str(len(passed)))
    print('Failed in Finding a Correct Surface Pressure: '+str(len(failed_in_finding_surfP)))
    print('Failed in an inner photochem loop after setting new surface pressure: '+str(len(failed_in_photochem_inner_surfP_loop)))
    print('Failed in photochem before setting a new surface pressure: '+str(len(failed_in_photochem_prior_to_surfP)))
    print('Failed in Climate: '+str(len(failed_in_climate)))
    print('Failed with at least 1 climate run completed: '+str(len(failed_with_atleast_1_clim_run)))
    print('Timed out (probably climate though): '+str(len(timedout)))
    print('Unaccounted for: '+str(len(unaccounted_for)))
    print('\n')


