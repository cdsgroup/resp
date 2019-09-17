import psi4
import resp

mol = psi4.geometry(""" C   1.45051389  -0.06628932   0.00000000
 H   1.75521613  -0.62865986  -0.87500146
 H   1.75521613  -0.62865986   0.87500146
 H   1.92173244   0.90485897   0.00000000
 C  -0.04233122   0.09849378   0.00000000
 O  -0.67064817  -1.07620915   0.00000000
 H  -1.60837259  -0.91016601   0.00000000
 O  -0.62675864   1.13160510   0.00000000""")
mol.update_geometry()

options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
           }

# Call for first stage fit
charges1 = resp.resp([mol], options)
print('Electrostatic Potential Charges')
print(charges1[0])
print('Restrained Electrostatic Potential Charges')
print(charges1[1])

# Change the value of the RESP parameter A
options['RESP_A'] = 0.001

# Add constraint for atoms fixed in second stage fit
constraint_charge = []
for i in range(4, 8):
    constraint_charge.append([charges1[1][i], [i+1]])
options['constraint_charge'] = constraint_charge
options['constraint_group'] = [[2, 3, 4]]
options['grid'] = ['1_%s_grid.dat' %mol.name()]
options['esp'] = ['1_%s_grid_esp.dat' %mol.name()]

# Call for second stage fit
charges2 = resp.resp([mol], options)

# Get RESP charges
print("\nStage Two:\n")
print('RESP Charges')
print(charges2[1])
