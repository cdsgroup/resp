import psi4
import resp

# Initialize three different conformations of ethanol
geometry = """C    0.00000000  0.00000000  0.00000000
C    1.48805540 -0.00728176  0.39653260
O    2.04971655  1.37648153  0.25604810
H    3.06429978  1.37151670  0.52641124
H    1.58679428 -0.33618761  1.43102358
H    2.03441010 -0.68906454 -0.25521028
H   -0.40814044 -1.00553466  0.10208540
H   -0.54635470  0.68178278  0.65174288
H   -0.09873888  0.32890585 -1.03449097
"""
mol1 = psi4.geometry(geometry)
mol1.update_geometry()
mol1.set_name('conformer1')

geometry = """C    0.00000000  0.00000000  0.00000000
C    1.48013500 -0.00724300  0.39442200
O    2.00696300  1.29224100  0.26232800
H    2.91547900  1.25572900  0.50972300
H    1.61500700 -0.32678000  1.45587700
H    2.07197500 -0.68695100 -0.26493400
H   -0.32500012  1.02293415 -0.30034094
H   -0.18892141 -0.68463906 -0.85893815
H   -0.64257065 -0.32709111  0.84987482
"""
mol2 = psi4.geometry(geometry)
mol2.update_geometry()
mol2.set_name('conformer2')

geometry = """C    0.00000000  0.00000000  0.00000000
C    1.48805540 -0.00728176  0.39653260
O    2.04971655  1.37648153  0.25604810
H    3.06429978  1.37151670  0.52641124
H    1.58679428 -0.33618761  1.43102358
H    2.03441010 -0.68906454 -0.25521028
H   -0.27127997 -0.97245518 -0.41089913
H   -0.60950912  0.20545158  0.87999334
H   -0.17244493  0.77215758 -0.74975690
"""
mol3 = psi4.geometry(geometry)
mol3.update_geometry()
mol3.set_name('conformer3')



molecules = [mol1, mol2, mol3]

# Specify intermolecular constraints
intermolecular_constraint = {'EQUAL': [[[1, range(1, 10)], [2, range(1, 10)]], 
                                       [[1, range(1, 10)], [3, range(1, 10)]]]}

# Specify options
options1 = {'VDW_SCALE_FACTORS' : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
           'RESTRAINT'          : True,
           'IHFREE'             : False,
           'WEIGHT'             : 1,
           }
options2 = {'WEIGHT': 1}
options3 = {'WEIGHT': 1}
options = [options1, options2, options3]

# Call for first stage fit
charges1 = resp.resp(molecules, options, intermolecular_constraint)

print("Restrained Electrostatic Potential Charges")
print(charges1[0][1])


# Add constraint for atoms fixed in second stage fit
for mol in range(len(molecules)):
    resp.set_stage2_constraint(molecules[mol], charges1[mol][1], options[mol], cutoff=1.2)
    options[mol]['grid'] = '%i_%s_grid.dat' %(mol+1, molecules[mol].name())
    options[mol]['esp'] = '%i_%s_grid_esp.dat' %(mol+1, molecules[mol].name())
    options[0]['resp_a'] = 0.001
    molecules[mol].set_name('conformer' + str(mol+1) + '_stage2')

# Add intermolecular constraints
intermolecular_constraint = resp.stage2_intermolecular_constraint(molecules, cutoff=1.2)

# Call for second stage fit
charges2 = resp.resp(molecules, options, intermolecular_constraint)
print("\nStage Two\n")
print("RESP Charges")
print(charges2[0][1])
