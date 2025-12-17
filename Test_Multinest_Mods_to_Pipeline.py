from Parameter_Sweep import * 
from multiprocessing import Pool


test_object = Generate_Atmosphere_Parameter_Sweep('MultinestModsT5', 
                                'Na', 
                                restart_run= 'Cco2', 
                                starting_point='Euclidean',
                                hitran_year='2020',
                                climate2col=False,
                                spectra=True,
                                planet='T1c')

# last test pretty good - Runs 2, 3, 6, 7 converged; runs 1, 4, 5 did not

# Five Tests
# watflx, co2flx, so2fixmr, oveff, o2veff, o3vdep, h2o2vdep, co2veff, covdep
input_fluxes = [[759600000000.0, 0, 0, 0.01, 0.01, 0.02, 0.02, 0, 0.03, 1], # Equal inputs to Czoom1/RunNumber46
                [759600000000.0, 0.0, 0.0001, 0.01, 0.01, 0.02, 0.02, 0, 0.03, 2], # Equal inputs to AdjSO2/Czoom1/RunNumber19 
                [947900000000.0, 1000.0, 1e-6, 0.01, 0.00678, 0.02, 0.02, 0.001, 0.0000001, 3], # Randomized inputs, low CO2 / SO2
                [900000000000.0, 10.0, 1e-4, 0.0135, 0.0896, 0.02, 0.02, 0.1, 0.002, 4], # Perturbed from Cinit/RunNumber27
                [608710000000.0, 10000, 1e-10, 0.01, 0.044, 0.02, 0.02, 0.088, 0.0001, 5] # Perturbed from Czoom2/RunNumber45
                ] 

test_object.mcmc_pressure_only = True
test_object.multinest_fit_data = True

with Pool() as p:
    models = p.map(test_object.run_one_model, input_fluxes)

