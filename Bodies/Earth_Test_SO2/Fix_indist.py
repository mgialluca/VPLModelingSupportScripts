import numpy as np
from astropy.io import ascii
import pandas as pd
from astropy.table import Table, Column

### Get all of the gases in the template in.dist file

#asp = open('../../../Atmos_Dev/atmos/PHOTOCHEM/INPUTFILES/TEMPLATES/Archean/species.dat', 'r')
asp = open('./species.dat', 'r')

done = False
old_gases = []
for l in asp.readlines():
    if l.split()[0][0] == '*':
        if done == True:
            break
    else:
        old_gases.append(l.split()[0])
        done = True

### Get all of the desired gases in the new atmosphere

nsp = open('./species.dat', 'r')

done = False
new_gases = []
for l in nsp.readlines():
    if l.split()[0][0] == '*':
        if done == True:
            break
    else:
        new_gases.append(l.split()[0])
        done = True

#### Figure out which gases are kept and lost

keep_gas = []
for g in old_gases:
    if g in new_gases:
        keep_gas.append(True)
    else:
        keep_gas.append(False)

## We are assuming NZ old = NZ new (in this case 200)
NZ = 200

NQ_old = 41
NQ_new = len(new_gases)

# Find number of blocks of mixing ratios until T/EDD columns
NQblocks_old = np.ceil(NQ_old/10)
NQblocks_new = np.ceil(NQ_new/10)

### Now to create new in.dist

fnew = open('NEWin.dist', 'w')

# Keep track of the lines that are starting and ending the current block in the in.dist file
blockstart = 0
blockend = NZ

flagstart = 0
flagend = 10

olddist = './in.dist'

write_table = Table()
gas_counter = 0

# this should handle all the mixing ratio blocks
for i in range(int(NQblocks_old)):
    curr_nq_block = ascii.read(olddist, data_start=blockstart, data_end=blockend, header_start=None)
    blockstart = blockend
    blockend = blockend+NZ

    curr_flags = keep_gas[flagstart:flagend]
    flagstart = flagend
    flagend = flagend+10        

    for c in range(len(curr_flags)):

        if len(write_table.columns) == 10:
            for line in range(len(write_table)):
                fnew.write('   ')
                for col in write_table.columns:
                    fnew.write("{:.8E}".format(write_table[col][line])+'   ')
                fnew.write('\n')

            write_table = Table()
        
        if curr_flags[c] == True:
            write_table.add_column(curr_nq_block['col'+str(c+1)], name='col'+str(gas_counter))
            if old_gases[gas_counter] == 'S8AER':
                for d in range(len(write_table['col'+str(gas_counter)])):
                    write_table['col'+str(gas_counter)][d] = 0
            gas_counter += 1

if len(write_table.columns) > 0:
    for line in range(len(write_table)):
        fnew.write('   ')
        for col in write_table.columns:
            fnew.write("{:.8E}".format(write_table[col][line])+'   ')
        fnew.write('\n')

# now the T/EDD/DEN/O3/CO2 block
T_edd_block = ascii.read(olddist, data_start=blockstart, data_end=blockend)
blockstart = blockend

T_edd_block['col5'] = 0 # replace CO2 ndens with 0

for line in range(len(T_edd_block)):
    fnew.write('   ')
    for col in T_edd_block.columns:
        fnew.write("{:.8E}".format(T_edd_block[col][line])+'   ')
    fnew.write('\n')

# colums 1 through 4 are aersol number dens - we want first 2 columns
# columns 5 through 8 are particle fall velocities - want col 5 /6 
# columns 9 through 12 are particle radii - want col 9 / 10 
## should wind up with col 1 and 2 (ndens), 3 and 4 (fall velocities), 5 and 6 (part radii)

# now the rest of in.dist
aers_block = ascii.read(olddist, data_start=blockstart)
write_table = Table()
#write_table.add_column(aers_block['col1'])
#write_table.add_column(aers_block['col2'])
#write_table.add_column(aers_block['col5'])
#write_table.add_column(aers_block['col6'])
#write_table.add_column(aers_block['col9'])
#write_table.add_column(aers_block['col10'])

write_table.add_column(aers_block['col1'])
write_table.add_column(aers_block['col2'])
write_table['col2'] = 0
write_table.add_column(aers_block['col3'])
write_table.add_column(aers_block['col4'])
write_table.add_column(aers_block['col5'])
write_table.add_column(aers_block['col6'])

for line in range(len(write_table)):
    fnew.write('   ')
    for col in write_table.columns:
        fnew.write("{:.8E}".format(write_table[col][line])+'   ')
    fnew.write('\n')

fnew.close()

