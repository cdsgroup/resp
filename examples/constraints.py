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

# There are examples on how to specify intramolecular
# constraints in the code

###################################################################

# EXAMPLE 1

# Specify intramolcular constraints
# Constrain atoms 1 and 2 to have a total charge 0.5
# Constrain atom 3 to have a charge of -1.0
options = {'constraint_charge': [[0.5, [1, 2]], [-1.0, [3]]]}
charges1 = resp.resp(molecules, options)

###################################################################

# EXAMPLE 2

# Constrain atom 1, 2, and 3 to have equal charges
# Also constrain atoms 4 and 5 to have equal charges
options = {'constraint_group': [[1, 2, 3], [4, 5]]}
charges2 = resp.resp(molecules, options)

###################################################################

# EXAMPLE 3

# The RESP module provides helper script for intramolecular charge constraints
# for the standard two-stage fitting procedure
options = {}

# No costraints
charges3_1 = resp.resp(molecules, options)

# Add constraint for atoms fixed in second stage fit
# Aliphatic carbons and hydrogens are refitted in the second stage
# Additionally, hydrogen atoms connected to the same carbon atom are made to
# have equal charges
resp.set_stage2_constraint(molecules[0], charges3_1[1], options)
charges3_2 = resp.resp(molecules, options)

###################################################################

# Print results
print("\nResults 1: ", charges1[1])
print('\n', '-'*100, '\n')

print("\nResults 2: ", charges2[1])
print('\n', '-'*100, '\n')

print("\nResults 3 stage 1: ", charges3_1[1])
print("\nResults 3 stage 2: ", charges3_2[1])
print('\n', '-'*100, '\n')
