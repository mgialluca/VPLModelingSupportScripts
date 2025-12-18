from Parameter_Sweep import *
import subprocess


test_object = Generate_Atmosphere_Parameter_Sweep('T1cMN', 
                                'NA', 
                                restart_run='Cco2', 
                                starting_point='Euclidean',
                                hitran_year='2020',
                                climate2col=False,
                                spectra=True,
                                planet='T1c')

test_object.match_data_multinest()
