
import numpy as np
import astropy.constants as const

### NOTE - Most of these depend on "indistdic" which is the ingested out.dist file turned into a...
###  ... python dictionary found from the ingest_indist(self) function in the pipeline

# Get the true number densities of each species based on their reported mixing ratios
def get_true_number_densities(indistdic):

    d = {}
    for i in indistdic.keys():
        if i not in ['Temp', 'Edd', 'NDens']:
            newndens = [indistdic['NDens'][lvl]*indistdic[i][lvl] for lvl in range(len(indistdic[i]))]
            d[i] = newndens

    return d

# Find the sum of mixing ratios as a function of altitude
# NZ - the number of layers used in photochem (self.nlevel_fine in the pipeline)
# n2mixingrat - fixed mixing ratio of N2
def sum_mixing_ratios(indistdic, NZ, n2mixingrat):

    sumsarr = np.zeros(NZ)

    for lvl in range(len(sumsarr)):
        for i in indistdic.keys():
            if i not in ['Temp', 'Edd', 'NDens']: # Exclude these from the calculation for now ig
                sumsarr[lvl] = sumsarr[lvl] + indistdic[i][lvl]

        # Add the N2
        sumsarr[lvl] = sumsarr[lvl] + n2mixingrat

    return sumsarr

# Get total number density of new atm and the change in number density on that level
def new_total_Ndens(mixings, Ndens):

    new_Ndens = []
    change_in_Ndens = []
    for lvl in range(len(Ndens)):
        new_Ndens_hold = Ndens[lvl]*mixings[lvl]
        change_in_Ndens.append((np.absolute(Ndens[lvl] - new_Ndens_hold)/Ndens[lvl]))
        new_Ndens.append(new_Ndens_hold)

    return new_Ndens, change_in_Ndens

# Get the new mixing ratios for each species
def new_mixing_rats(Ndenses, NdensTot):

    d = {}
    for i in Ndenses.keys():
        d[i] = [Ndenses[i][lvl]/NdensTot[lvl] for lvl in range(len(NdensTot))]

    return d

# Find column mass density of atmosphere
def find_tot_column_mass_dens(Ndenstot, alt, mmw):

    colnumdens = np.trapz(Ndenstot, alt)
    colmas = (1e-3*colnumdens*mmw)/const.N_A.value

    return colmas