import psi4
import resp

# This example was contributed by Karl Kirschner

# ethane-1,2-diol
# Global mimimum (i.e. rel. E = 0.0 kcal/mol) has an
#     intramolecular hydrogen bond.
# Determine charges for conformation without intramolecular
#     hydrogen bond.

# Conf1. HF/6-31G(d)//HF/6-31G(d); rel. E = 2.1 kcal/mol
geometry = """C     0.078467344400     0.008743412200     0.095163192300
H     0.698032751600     0.819309173300    -0.279081297200
H    -0.709956225400    -0.174483577300    -0.629925270600
C     0.922732830900    -1.238543286900     0.240036856600
H     0.303167256400    -2.049109212500     0.614280738400
H     1.711156054700    -1.055316507100     0.965125744700
O     1.447389951300    -1.534620648100    -1.026823868100
H     1.983856853400    -2.312390703000    -0.976649432200
O    -0.446189165900     0.304821153500     1.362024062900
H    -0.982657651400     1.082590195900     1.311849273200
"""
conf1 = psi4.geometry(geometry)
conf1.update_geometry()
conf1.set_name('conformer1')

# Conf2. HF/6-31G(d)//HF/6-31G(d); rel. E = 2.4 kcal/mol
geometry = """C     0.087341104900    -0.008801533700     0.124594260200
H    -0.653760471900     0.010564681500    -0.669677267800
H     0.447760507200    -1.024108696200     0.213349830500
C     1.242702680900     0.907616403200    -0.235629753300
H     2.004683470100     0.843318180600     0.536598411000
H     0.895386596100     1.940251823600    -0.273528261600
O     1.732819730000     0.495058244800    -1.482260250000
H     2.495682683700     1.004676191100    -1.714109526400
O    -0.471405729500     0.337055538700     1.363487493400
H    -0.996210571400     1.119369166500     1.269175064100
"""
conf2 = psi4.geometry(geometry)
conf2.update_geometry()
conf2.set_name('conformer2')

# Conf3. HF/6-31G(d)//HF/6-31G(d); rel. E = 2.5 kcal/mol
geometry = """C     2.366829887600     0.675870944700     0.826419334800
H     2.677198052700    -0.365650047900     0.832645967800
H     3.240760885100     1.274109837100     0.596053933500
C     1.299970661200     0.883928381500    -0.239219141900
H     0.989601941900     1.925449267900    -0.245445092200
H     0.426039927000     0.285688946300    -0.008853953100
O     1.743280010700     0.466329378200    -1.501950606300
H     2.370350061500     1.089512391700    -1.841358439300
O     1.923520556800     1.093470931800     2.089150421800
H     1.296448015400     0.470289968800     2.428557574800
"""
conf3 = psi4.geometry(geometry)
conf3.update_geometry()
conf3.set_name('conformer3')

# Conf4. HF/6-31G(d)//HF/6-31G(d); rel. E = 2.9 kcal/mol
geometry = """C     0.860437394600     2.415873818900    -0.341957689300
H     1.704477000500     3.019997018700    -0.675719391600
H     0.068180619600     2.537578474800    -1.067326537600
C     1.261867500500     0.949056277300    -0.273498319600
H     0.400458948700     0.353443205800    -0.005425084100
H     2.014363546300     0.802469243300     0.501760766200
O     1.713720202000     0.476649292700    -1.512916351600
H     2.557153945700     0.857355849600    -1.713984689700
O     0.363233530000     2.875283966400     0.884919129500
H     1.068107312100     2.938292852400     1.514148167800
"""
conf4 = psi4.geometry(geometry)
conf4.update_geometry()
conf4.set_name('conformer4')

conformers = [conf1, conf2, conf3, conf4]

# Specify options
options = {'VDW_SCALE_FACTORS' : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
           'RESTRAINT'          : True,
           'IHFREE'             : False,
           'WEIGHT'             : [1, 1, 1, 1],
           }
# Specify additional intramolecular constraints
# (i.e. both O (ids 7 and 9) should have equal charges,
#   both hydroxyl H (ids 8 and 10) should equal charges,
#   both C (ids 1 and 4) should have equal charges, and
#   all four aliphatic H (ids 2, 3, 5 and 6) should have equal charges)
options["constraint_group"] = [[7, 9], [8, 10], [1, 4], [2, 3, 5, 6]]

# Call for first stage fit
charges1 = resp.resp(conformers, options)

print("Restrained Electrostatic Potential Charges")
print(charges1[1])

resp.set_stage2_constraint(conformers[0], charges1[1], options)
options['RESP_A'] = 0.001

options['grid'] = []
options['esp'] = []

# Add constraint for atoms fixed in second stage fit
for structure in range(len(conformers)):
    options['grid'].append('%i_%s_grid.dat' % (structure + 1, conformers[structure].name()))
    options['esp'].append('%i_%s_grid_esp.dat' % (structure + 1, conformers[structure].name()))

# Add intermolecular constraints
# Specify additional intramolecular constraints (only need to do for 1 molecule)
#   both C (ids 1 and 4) should have equal charges, and
#   all four aliphatic H (ids 2, 3, 5 and 6) should have equal charges)
options["constraint_group"] = [[1, 4], [2, 3, 5, 6]]

# Call for second stage fit
charges2 = resp.resp(conformers, options)
print("\nStage Two\n")
print("RESP Charges")
print(charges2[1])
