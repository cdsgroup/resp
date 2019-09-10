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

# There are 5 examples on how to specify intra- and inter-molecular
# constraints in the code

###################################################################

# EXAMPLE 1
molecules = [mol1, mol2, mol3]
options = [{}, {}, {}]

# Specify intermolecular constraints
# Make equivalent atoms in the three conformer have identical charges
intermolecular_constraint = {'EQUAL': [[[1, range(1, 10)], [2, range(1, 10)]], 
                                       [[1, range(1, 10)], [3, range(1, 10)]]]}

# Specify intramolcular constraints
# Constrain atoms 1 and 2 to have a total charge 0.5
# Constrain atom 3 to have a charge of -1.0
options[0] = {'constraint_charge': [[0.5, [1, 2]], [-1.0, [3]]]}
charges1 = resp.resp(molecules, options, intermolecular_constraint)

###################################################################

# EXAMPLE 2
molecules = [mol1, mol2, mol3]
options = [{}, {}, {}]

# Specify intermolecular constraints
# Make equivalent atoms in the three conformer have identical charges
intermolecular_constraint = {'EQUAL': [[[1, range(1, 10)], [2, range(1, 10)]], 
                                       [[1, range(1, 10)], [3, range(1, 10)]]]}


# Constrain atom 1, 2, and 3 to have equal charges
# Also constrain atoms 4 and 5 to have equal charges
options[0] = {'constraint_group': [[1, 2, 3], [4, 5]]}
charges2 = resp.resp(molecules, options, intermolecular_constraint)

###################################################################

# EXAMPLE 3
molecules = [mol1, mol2, mol3]
options = [{}, {}, {}]

# Other possible intermolecular constraints
# The sum of the charges on atoms 1 and 2 of molecule 1 and atoms 3 and 4
# of molecule 2 will equal 1.0
# In this case, conformer three will not have any constraints
intermolecular_constraint = {'CHARGE': [[1.0, [[1, [1, 2]], [2, [3, 4]]]]]}
charges3 = resp.resp(molecules, options, intermolecular_constraint)

###################################################################

# EXAMPLE 4

# The RESP module provides helper script for intramolecular charge constraints
# for the standard two-stage fitting procedure
molecules = [mol1]
options = [{}]

# No costraints
charges4_1 = resp.resp(molecules, options)

# Add constraint for atoms fixed in second stage fit
# All atoms are fixed in the second stage except carbon atoms that have hydrogen atoms
# connected to those carbons. The connection is determined by the cutoff distance
# Additionally, hydrogen atoms connected to the same carbon atom are made to
# have equal charges
resp.set_stage2_constraint(molecules[0], charges4_1[0][1], options[0], cutoff=1.2)
charges4_2 = resp.resp(molecules, options)

###################################################################

# EXAMPLE 5
# The RESP module provides helper script for intermolecular charge constraints
# in the second stage fitting
molecules = [mol1, mol2, mol3]
options = [{}, {}, {}]

# Specify intermolecular constraints
# Make equivalent atoms in the three conformer have identical charges
intermolecular_constraint = {'EQUAL': [[[1, range(1, 10)], [2, range(1, 10)]], 
                                       [[1, range(1, 10)], [3, range(1, 10)]]]}

# No intramolecular costraints
charges5_1 = resp.resp(molecules, options, intermolecular_constraint)

# Set intramolecular constraints for second stage atoms
for i in range(len(molecules)):
    resp.set_stage2_constraint(molecules[i], charges5_1[i][1], options[i], cutoff=1.2)

# Set intermolecular constraints for second stage atoms
# If we used the previous intermolecular_constraint, we would get a singular matrix error.
# This is because we introduced more constraints than necessary to specify the charges of
# the molecule. The helper script specifies only the necessary constraints in this case 
intermolecular_constraint = resp.stage2_intermolecular_constraint(molecules, cutoff=1.2)
charges5_2 = resp.resp(molecules, options, intermolecular_constraint)

###################################################################

# Print results
print("\nResults 1: ", charges1[0][1])
print('\n', '-'*100, '\n')

print("\nResults 2: ", charges2[0][1])
print('\n', '-'*100, '\n')

print("\nResults 3 molecule 1: ", charges3[0][1])
print("\nResults 3 molecule 2: ", charges3[1][1])
print("\nResults 3 molecule 3: ", charges3[2][1])
print('\n', '-'*100, '\n')

print("\nResults 4 stage 1: ", charges4_1[0][1])
print("\nResults 4 stage 2: ", charges4_2[0][1])
print('\n', '-'*100, '\n')

print("\nResults 5 stage 1: ", charges5_1[0][1])
print("\nResults 5 stage 2: ", charges5_2[0][1])
print('\n', '-'*100, '\n')
