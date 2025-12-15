from Parameter_Sweep import * 
from multiprocessing import Pool


test_object = Generate_Atmosphere_Parameter_Sweep('TestMultinestMods', 
                                'Na', 
                                restart_run= 'Cco2', 
                                starting_point='Euclidean',
                                hitran_year='2020',
                                climate2col=False,
                                spectra=True,
                                planet='T1c')


# Five Tests
# watflx, co2flx, so2fixmr, oveff, o2veff, o3vdep, h2o2vdep, co2veff, covdep
input_fluxes = [[759600000000.0, 0, 0, 0.01, 0.01, 0.02, 0.02, 0, 0.03, 1], # Equal inputs to Czoom1/RunNumber46
                [947900000000.0, 88539.0, 0, 0.01, 0.01, 0.02, 0.02, 0.01, 0, 2], # Equal inputs to Cco2wL/RunNumber3
                [947900000000.0, 88539.0, 0.001, 0.05, 0.05, 0.02, 0.02, 0.0, 0.0, 3], # Equal inputs to AdjSO2/Cco2/RunNumber17
                [759600000000.0, 0.0, 0.0001, 0.01, 0.01, 0.02, 0.02, 0, 0.03, 4], # Equal inputs to AdjSO2/Czoom1/RunNumber19 
                [947900000000.0, 1000.0, 1e-6, 0.01, 0.00678, 0.02, 0.02, 0.001, 0.0000001, 5]] # Randomized inputs, low CO2 / SO2

test_object.mcmc_pressure_only = True
test_object.multinest_fit_data = True

with Pool() as p:
    models = p.map(test_object.run_one_model, input_fluxes)

