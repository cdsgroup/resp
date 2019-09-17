from __future__ import division, absolute_import, print_function

import numpy as np

import psi4

"""
A helper script to facilitate the use of constraints for two-stage fitting.
"""

def _get_stage2_atoms(molecule):
    """Determines atoms for second stage fit. The atoms
       are identified as sp3 carbons that have one or more hydrogens.

    Parameters
    ----------
    molecule : psi4.Molecule instance

    Returns
    -------
    groups : dict
        a dictionary whose keys are the indecies+1 of carbon
        atoms and whose elements are the indecies+1 of the
        connected hydrogen atoms.

    """

    symbols = []
    for i in range(molecule.natom()):
        symbols.append(molecule.symbol(i))

    groups = {}
    bond_profile = psi4.qcdb.parker._bond_profile(molecule)
    for i in range(molecule.natom()):
        # Find carbon atoms
        if symbols[i] != 'C':
            continue
        # Check that it has 4 bonds
        bonds_for_atom = [j for j in bond_profile if i in j[:2]]
        if len(bonds_for_atom) == 4:
            group = []
            for atoms in bonds_for_atom:
                j = atoms[0] if atoms[0] != i else atoms[1]
                if symbols[j] == 'H':
                    group.append(j + 1)  
            if group:
                groups[i + 1] = group

    return groups


def set_stage2_constraint(molecule, charges, options):
    """Sets default constraints for the second stage fit.

    The default constraints are the following:
    Atoms that are excluded from the second stage fit are constrained
    to their charges from the first stage fit. C-H groups determined 
    by _get_stage2_atoms are refitted and the hydrogen atoms connected
    to the same carbon are constrained to have identical charges.

    Parameters
    ----------
    molecule : psi4.core.Molecule

    charges : :py:class:`numpy.ndarray`
        array containing the charges from the first stage fit
    options : dict
        dictionary of the fitting options. To be modified in place.
    cutoff : float, optional
        cutoff distance in Angstroms, exclusive

    Return
    ------
    None

    """
    second_stage = _get_stage2_atoms(molecule)
    atoms = list(range(1, molecule.natom()+1))
    constraint_group = []
    for i in second_stage.keys():
        atoms.remove(i)
        group = []
        for j in second_stage[i]:
            atoms.remove(j)
            group.append(j)
        constraint_group.append(group)
    constraint_charge = []
    for i in atoms:
        constraint_charge.append([charges[i-1], [i]])
    options['constraint_charge'] = constraint_charge
    options['constraint_group'] = constraint_group
