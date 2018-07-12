"""
Driver for the RESP code.
"""
from __future__ import division, absolute_import, print_function

__authors__   =  "Asim Alenaizan"
__credits__   =  ["Asim Alenaizan"]

__copyright__ = "(c) 2014-2018, The Psi4NumPy Developers"
__license__   = "BSD-3-Clause"
__date__      = "2018-04-28"

import subprocess

import numpy as np
import qcdb

from . import espfit
from . import vdw_surface_helper


def resp(name, molecules, intermol_constraints=None, **kwargs):
    """RESP code driver.

    Parameters
    ---------- 
    name : str
        method for computing electrostatic potential
    molecules : list
        list of qcdb.Molecule instances
    options_list : list, optional
        list of dictionaries of user's defined options
    intermol_constraints : dict, optional
        dictionary of options for multi-molecules fitting

    Returns
    -------
    charges : list
        list of charges

    Note
    ----
    output files : mol_results.dat: fitting results
                   mol_grid.dat: grid points in molecule.units
                   mol_grid_esp.dat: QM esp valuese in a.u. 
    """

    if not isinstance(molecules, (list, tuple)):
        molecules = [molecules]

    if intermol_constraints is None:
        intermol_constraints = {}

    # Check options
    intermol_constraints = {k.upper(): v for k, v in intermol_constraints.items()}

    if 'CHARGE' not in intermol_constraints:
        intermol_constraints['CHARGE'] = [] 
    if 'EQUAL' not in intermol_constraints:
        intermol_constraints['EQUAL'] = []

    # Duplicate options if options are not provided for all molecules
    qopt = qcdb.get_active_options().scroll['RESP']
    for opt in qopt:
        if opt in ['RESTRAINT', 'MAX_IT', 'TOLER']: continue 
        if len(qopt[opt].value) != len(molecules):
            qcdb.set_options({'RESP_%s' %opt: qopt[opt].value*len(molecules)})
    qopt = qcdb.get_active_options().scroll['RESP']

    data = {'mol_charge': [], 'coordinates': [], 'symbols': [], 'n_atoms': [], 'invr': [], 'esp_values': []}
    for imol in range(len(molecules)):
        radii = {}
        for i in qopt['RADIUS'].value[imol]:
            radii[i.upper()] = qopt['RADIUS'].value[imol][i]
        qopt['RADIUS'].value[imol] = radii
        data['mol_charge'].append(molecules[imol].molecular_charge())
        data['n_atoms'].append(molecules[imol].natom())
        arrays = molecules[imol].to_arrays()
        coordinates = arrays[0]*qcdb.physconst.psi_bohr2angstroms
        data['coordinates'].append(coordinates)
        symbols = arrays[2]
        data['symbols'].append(symbols)

        if qopt['GRID'].value[imol]:
            # Read grid points
            points = np.loadtxt(qopt['GRID'].value[imol])
            np.savetxt('grid.dat', points, fmt='%15.10f')
            if 'Bohr' == molecules[imol].units():
                points *= qcdb.physconst.psi_bohr2angstroms

        else:
            # Get the points at which we are going to calculate the ESP surface
            points = []
            surface = vdw_surface_helper.vdw_surface_helper()
            for i in range(qopt['N_VDW_LAYERS'].value[imol]):
                scale_factor = qopt['VDW_SCALE_FACTOR'].value[imol] + i * qopt['VDW_INCREMENT'].value[imol]
                surface.vdw_surface(coordinates, symbols, scale_factor,
                                    qopt['VDW_POINT_DENSITY'].value[imol], qopt['RADIUS'].value[imol])
                points.append(surface.shell)
            radii = surface.radii
            points = np.concatenate(points)
            if 'Bohr' == molecules[imol].units():
                points /= qcdb.physconst.psi_bohr2angstroms
                np.savetxt('grid.dat', points, fmt='%15.10f')
                points *= qcdb.physconst.psi_bohr2angstroms
            else:
                np.savetxt('grid.dat', points, fmt='%15.10f')

        # Calculate ESP values at the grid
        if qopt['ESP'].value[imol]:
            # Read electrostatic potential values
            data['esp_values'].append(np.loadtxt(qopt['ESP'].value[imol]))
            np.savetxt('grid_esp.dat', data['esp_values'][imol], fmt='%15.10f')
        else:
            if qopt['QM_PACKAGE'].value[imol] == 'PSI4':
                qcdb.get_active_options().require('PSI4', 'GRIDDAT', open('grid.dat').read(), accession='wert')
                e, jrec = qcdb.properties(name, molecule=molecules[imol], properties=['GRID_ESP'], return_wfn=True)
                with open('grid_esp.dat', 'w') as handle:
                    handle.write(jrec['outfile_grid_esp.dat'])
                data['esp_values'].append(np.loadtxt('grid_esp.dat'))

            elif qopt['QM_PACKAGE'].value[imol] == 'Q-CHEM':
                from .qm_helper import qchem
                data['esp_values'].append(qchem.qchem_esp(molecules[imol], name, qcdb.get_active_options().scroll['QCDB']['BASIS'].value))

            elif qopt['QM_PACKAGE'].value[imol] == 'GAMESS':
                from .qm_helper import gamess
                data['esp_values'].append(gamess.gamess_esp(molecules[imol], name, qcdb.get_active_options().scroll['QCDB']['BASIS'].value)) 

        subprocess.run(['mv', 'grid.dat', '%i_%s_grid.dat' %(imol+1, molecules[imol].name())])
        subprocess.run(['mv', 'grid_esp.dat', '%i_%s_grid_esp.dat' %(imol+1, molecules[imol].name())])

        # Build a matrix of the inverse distances from each ESP point to each nucleus
        invr = np.zeros((len(points), len(coordinates)))
        for i in range(invr.shape[0]):
            for j in range(invr.shape[1]):
                invr[i, j] = 1/np.linalg.norm(points[i]-coordinates[j])
        data['invr'].append(invr*qcdb.physconst.psi_bohr2angstroms) # convert to atomic units
        data['coordinates'][imol] /= qcdb.physconst.psi_bohr2angstroms # convert to angstroms

    # Calculate charges
    qf, labelf, notes = espfit.fit(qopt, data, intermol_constraints)
    index = 0
    charges = []
    
    # Extract the charges
    for imol in range(len(molecules)):
        q = [i[index:index+data['n_atoms'][imol]] for i in qf]
        index += data['n_atoms'][imol]
        charges.append(q)

    for imol in range(len(molecules)):
        # Write the results to disk
        with open(str(imol+1) + '_' + molecules[imol].name() + "_results.out", "w") as f:
            f.write("\n Electrostatic potential parameters\n")
            f.write("\n Grid information (see %i_%s_grid.dat in %s)\n"
                    %(imol+1, molecules[imol].name(), molecules[imol].units()))
            if not qopt['GRID'].value[imol]: 
                f.write("     van der Waals radii (Angstrom):\n")
                for i, j in radii.items():
                    f.write("                                %8s%8.3f\n" %(i, j/scale_factor))
                f.write("     Number of VDW layers:             %d\n" %(qopt["N_VDW_LAYERS"].value[imol]))
                f.write("     VDW scale facotr:                 %.3f\n" %(qopt["VDW_SCALE_FACTOR"].value[imol]))
                f.write("     VDW increment:                    %.3f\n" %(qopt["VDW_INCREMENT"].value[imol]))
                f.write("     VDW point density:                %.3f\n" %(qopt["VDW_POINT_DENSITY"].value[imol]))
                f.write("     Number of grid points:            %d\n" %len(data['esp_values'][imol]))

            f.write("\n Quantum electrostatic potential (see %i_%s_grid_esp.dat)\n" %(imol+1, molecules[imol].name()))
            if not qopt['ESP'].value[imol]:
                f.write("     QM package:                       %s\n" %qopt['QM_PACKAGE'].value[imol])
                f.write("     ESP method:                       %s\n" %name)
                f.write("     ESP basis set:                    %s\n" %qcdb.get_active_options().scroll['QCDB']['BASIS'].value)

            f.write("\n Constraints\n")
            if qopt['CONSTRAINT_CHARGE'].value[imol]:
                f.write("     Charge constraints\n")
                for i in qopt['CONSTRAINT_CHARGE'].value[imol]:
                    f.write("         Total charge of %8.5f on the set" %i[0])
                    for j in i[1]:
                        f.write("%4d" %j)
                    f.write("\n")
            if qopt['CONSTRAINT_GROUP'].value[imol] or qopt['CONSTRAINT_EQUAL'].value[imol]:
                f.write("     Equality constraints\n")
                f.write("         Equal charges on atoms\n")
                for i in qopt['CONSTRAINT_GROUP'].value[imol]:
                    f.write("                              ")
                    for j in i:
                        f.write("%4d" %j)
                    f.write("\n")
                for i in qopt['CONSTRAINT_EQUAL'].value[imol]:
                    for j in range(len(i)):
                        f.write("                              ")
                        f.write("%4d%4d" %(i[0][j], i[1][j]))
                        f.write("\n")
            if intermol_constraints['CHARGE'] or intermol_constraints['EQUAL']:
                f.write('\n     Intermolecular constraints\n')
                if intermol_constraints['CHARGE']:
                    f.write('         Charge constraints\n')
                    for i in intermol_constraints['CHARGE']:
                        f.write('             Total charge of %8.5f on the set:' %i[0])
                        for j in i[1]:
                            f.write('\n                 molecule %4d, atoms' %j[0])
                            for k in j[1]:
                                f.write('%4d' %k)
                        f.write('\n')
                if intermol_constraints['EQUAL']:
                    f.write('         Equality constraints\n')
                    f.write('             Equal charges on\n')
                    for i in intermol_constraints['EQUAL']:
                        f.write('                 ')
                        f.write('molecule %4d, atoms' %i[0][0])
                        for j in i[0][1]:
                            f.write('%4d' %j)
                        f.write('\n                 molecule %4d, atoms' %i[1][0])
                        for j in i[1][1]:
                            f.write('%4d' %j)
                        f.write('\n\n')
            f.write("\n Restraint\n")
            if qopt['RESTRAINT'].value:
                f.write("     Hyperbolic restraint to a charge of zero\n")
                if qopt['IHFREE'].value[imol]:
                    f.write("     Hydrogen atoms are not restrained\n")
                f.write("     resp_a:                           %.4f\n" %(qopt["RESP_A"].value[imol]))
                f.write("     resp_b:                           %.4f\n" %(qopt["RESP_B"].value[imol]))
            f.write("\n Fit\n")
            for i in notes:
                if i:
                    f.write(i+'\n')
            f.write("\n Electrostatic Potential Charges\n")
            f.write("   Center  Symbol")
            for i in labelf:
                f.write("%10s       " %i)
            f.write("\n")
            for i in range(data['n_atoms'][imol]):
                f.write("   %5d    %s     " %(i+1, data['symbols'][imol][i]))
                for j in charges[imol]:
                    f.write("%16.8f" %j[i])
                f.write("\n")
            f.write(" Total Charge:    ")
            for i in charges[imol]:
                f.write("%16.8f" %np.sum(i))
            f.write('\n')

    return charges
