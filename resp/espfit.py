"""
Fitting procedure for RESP charges.

Reference:
Equations taken from [Bayly:93:10269].
"""
from __future__ import division, absolute_import, print_function

import copy
import warnings

import numpy as np


def esp_solve(A, B):
    """Solves for point charges: A*q = B

    Parameters
    ----------
    A : ndarray
        array of matrix A
    B : ndarray
        array of matrix B

    Return
    ------
    q : ndarray
        array of charges

    """

    q = np.linalg.solve(A, B)
    # Warning for near singular matrix
    # in case np.linalg.solve does not detect singularity
    if np.linalg.cond(A) > 1/np.finfo(A.dtype).eps:
        warnings.warn("Possible fit problem; singular matrix")

    return q


def restraint(q, A_unrestrained, resp_a, resp_b, ihfree, symbols, n_atoms, n_constraints):
    """Adds hyperbolic restraint to matrix A

    Parameters
    ----------
    q : ndarray
        array of charges
    A_unrestrained : ndarray
        array of unrestrained A matrix
    resp_a : list
        list of floats of restraint scale a for each molecule
    resp_b : list
        list of floats of restraint parabola tightness b for each molecule
    ihfree : list
        list of bools on whether hydrogen excluded or included in restraint for each molecule
    symbols : list
        list of arrays of element symbols for each molecule
    n_atoms : list
        list of the number of atoms in each molecule
    n_constraints : list
        list of the number of constraints for each molecule

    Returns
    -------
    a : ndarray
        restrained A array
    """

    # hyperbolic Restraint
    # [Bayly:93:10271] (Eqs. 10, 13)
    A = copy.deepcopy(A_unrestrained)
    n_mol = len(n_atoms)
    index = 0
    # loop over all molecules/conformers
    for mol in range(n_mol):
        for i in range(n_atoms[mol]):
            # if an element is not hydrogen or if hydrogens are to be restrained
            if not ihfree[mol] or symbols[mol][i] != 'H':
                A[index+i, index+i] = (A_unrestrained[index+i, index+i]
                + resp_a[mol]/np.sqrt(q[index+i]**2 + resp_b[mol]**2))
        # move the index to the beginning of the block of A corresponding to molecule mol
        index += n_atoms[mol] + n_constraints[mol]
    return A


def iterate(q, A_unrestrained, B, resp_a, resp_b, ihfree, symbols, toler,
            maxit, n_atoms, n_constraints, atom_indices):
    """Iterates the RESP fitting procedure

    Parameters
    ----------
    q : ndarray
        array of initial charges
    A_unrestrained : ndarray
        array of unrestrained A matrix
    B : ndarray
        array of matrix B
    resp_a : list
        list of floats of restraint scale a for each molecule
    resp_b : list
        list of floats of restraint parabola tightness b for each molecule
    ihfree : list
        list of bools on whether hydrogen excluded or included in restraint for each molecule
    symbols : list
        list of arrays of element symbols for each molecule
    toler : float
        tolerance for charges in the fitting
    maxit : int
        maximum number of iterations
    n_atoms : list
        list of the number of atoms in each molecule
    n_constraints : list
        list of the number of constraints for each molecule
    atom_indices : ndarray
        array of the indices for the atoms in the A and B matrices

    Returns
    -------
    q : ndarray
        array of the fitted charges

    """
    q_last = copy.deepcopy(q)
    niter, dif, note = 0, 2*toler, ''
    while dif > toler and niter < maxit:
        niter += 1
        A = restraint(q, A_unrestrained, resp_a, resp_b, ihfree,
                      symbols, n_atoms, n_constraints)
        q = esp_solve(A, B)
        # Extract vector elements that correspond to charges
        dif = np.sqrt(np.max((q[atom_indices] - q_last[atom_indices])**2))
        q_last = copy.deepcopy(q)

    if dif > toler:
        note += ('\nCharge fitting did not converge; ' + 
               'try increasing the maximum number of iterations to ' +
               '> %i.' %maxit)
    return q[atom_indices], note


def intramolecular_constraints(constraint_charge, constraint_equal, constraint_groups):
    """Extracts intramolecular constraints from user constraint input

    Parameters
    ----------
    constraint_charge : list
        list of lists of charges and atom indices list
        e.g. [[0, [1, 2]], [1, [3, 4]]]
        The sum of charges on 1 and 2 will equal 0
        The sum of charges on 3 and 4 will equal 1
    constraint_equal : list
        list of lists of two lists of indices of atoms to
        have equal charge element by element
        e.g. [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]
        atoms 1 and 3 will have equal charge
        atoms 2 and 4 will have equal charge
        and similarly for 5, 6, 7 and 8
    constraint_group : list
        list of lists of indices of atoms to have equal charge
        e.g. [[1, 2], [3, 4]]
        atoms 1 and 2 will have equal charge
        atoms 3 and 4 will have equal charge

    Returns
    -------
    constrained_charges : list
        list of fixed charges
    constrained_indices : list
        list of lists of indices of atoms in a constraint
        negative number before an index means
        the charge of that atom will be subtracted.

    Notes
    -----
    Atom indices starts with 1 not 0.
    Total charge constraint is added by default for the first molecule.

    """
    constrained_charges = []
    constrained_indices = []
    for i in constraint_charge:
        constrained_charges.append(i[0])
        group = []
        for k in i[1]:
            group.append(k)
        constrained_indices.append(group)

    for i in constraint_equal:
        for j in range(len(i[0])):
            group = []
            constrained_charges.append(0)
            group.append(-i[0][j])
            group.append(i[1][j])
            constrained_indices.append(group)

    for i in constraint_groups:
        for j in range(1, len(i)):
            group = []
            constrained_charges.append(0)
            group.append(-i[j-1])
            group.append(i[j])
            constrained_indices.append(group)
    return constrained_charges, constrained_indices


def intermolecular_constraints(constraint_charge, constraint_equal):
    """Extracts intermolecular constraints from user constraint input

    Parameters
    ----------
    constraint_charge : list
        list of list of lists of charges and atom indices list
        e.g. [[1, [[1, [1, 2]], [2, [3, 4]]]]]
        The sum of charges on atoms 1 and 2 of molecule 1
        and atoms 3 and 4 of molecule 2 will equal 1.
    constraint_equal : list
        list of list of list of indices of atoms to have
        equal charge in two molecules.
        e.g. [[[1, [1, 2]], [2, [3, 4]]]]
        charges on atoms 1 and 2 in molecule 1 will equal
        charges on  atoms 3 and 4 in molecule 2, respectively.

    Returns
    -------
    constrained_charges : list
        list of fixed charges
    constrained_indices : list
        list of lists of indices of atoms in a constraint
        negative number before an index means
        the charge of that atom will be subtracted.
    molecules : list
        list of lists of constrained molecules.

    Note
    ----
    Atom indices starts with 1 not 0

    """
    constrained_charges = []
    constrained_indices = []
    molecules = []
    for i in constraint_charge:
        constrained_charges.append(i[0])
        mol = []
        group_big = []
        for j in i[1]:
            mol.append(j[0])
            group = []
            for k in j[1]:
                group.append(k)
            group_big.append(group)
        constrained_indices.append(group_big)
        molecules.append(mol)

    for i in constraint_equal:
        for j in range(len(i[0][1])):
            molecules.append([i[0][0], i[1][0]])
            group = []
            constrained_charges.append(0)
            group.append([-i[0][1][j]])
            group.append([i[1][1][j]])
            constrained_indices.append(group)
    return constrained_charges, constrained_indices, molecules


def fit(options, data, inter_constraint):
    """Performs ESP and RESP fits.

    Parameters
    ----------
    options : list
        list of dictionaries of fitting options and internal data
    inter_constraint : dict
        dictionary of user-defined intermolecular constraints.

    Returns
    -------
    qf : list
        list of ndarrays of fitted charges
    labelf : list
        list of strings of fitting methods i.e. ESP and RESP
    note : str
        string of notes on the fitting

    """
    rest = options[0]['RESTRAINT']
    n_mols = len(options)
    qf = []
    labelf = []
    invr, coordinates, n_constraint, symbols, n_atoms = [], [], [], [], []
    constrained_charges, constrained_indices = [], []
    ndim = 0
    con_charges_sys, con_indices_sys, con_mol_sys = intermolecular_constraints(
                                                     inter_constraint['CHARGE'],
                                                     inter_constraint['EQUAL'])
    n_sys_constraint = len(con_charges_sys)
    for mol in range(n_mols):
        invr.append(data[mol]['invr'])
        coordinates.append(data[mol]['coordinates'])
        symbols.append(data[mol]['symbols'])
        n_atoms.append(len(symbols[mol]))
        constraint_charge = options[mol]['CONSTRAINT_CHARGE']
        constraint_equal = options[mol]['CONSTRAINT_EQUAL']
        constraint_groups = options[mol]['CONSTRAINT_GROUP']
        # Get user-defined constraints
        charges, indices = intramolecular_constraints(constraint_charge,
                                                       constraint_equal,
                                                       constraint_groups)
        constrained_charges.append(charges)
        constrained_indices.append(indices)
        n_con = len(charges)
        if mol == 0:
            n_con += 1
        n_constraint.append(n_con)
        ndim += n_atoms[mol] + n_constraint[mol]
    n_atoms = np.array(n_atoms)
    n_constraint = np.array(n_constraint)
    symbols = np.array(symbols, dtype='str')
    # Additional constraint to make charges in different molecules equal
    # to the charge in the first molecule
    # Also, Total charges = molecular charge
    ndim += n_sys_constraint
    A = np.zeros((ndim, ndim))
    B = np.zeros(ndim)

    edges_i = 0
    edges_f = 0
    indices = []
    # Bayly:93:10271 (Eqs. 12-14)
    for mol in range(n_mols):
        indices.append(range(edges_i, edges_i+n_atoms[mol]))
        r_inverse, V = invr[mol], data[mol]['esp_values']

        # Lower case a and b are the A matrix and B vector for one molecule
        # and without the addition of constraints

        # Construct a: a_jk = sum_i [(1/r_ij)*(1/r_ik)]
        a = np.einsum("ij, ik -> jk", r_inverse, r_inverse)

        # Construct b: b_j = sum_i (V_i/r_ij)
        b = np.einsum('i, ij->j', V, r_inverse)

        # Weight the moleule 
        a *= options[mol]['WEIGHT']**2
        b *= options[mol]['WEIGHT']**2

        edges_f += n_atoms[mol]

        # set elements in A and B to a and b
        A[edges_i:n_atoms[mol]+edges_i, edges_i:n_atoms[mol]+edges_i] = a
        B[edges_i:n_atoms[mol]+edges_i] = b

        # Sum of point charges = molecular charge
        if mol == 0:
            molecular_charge_constraint, molecular_charge = n_atoms[0], data[0]['mol_charge']
            A[:molecular_charge_constraint, molecular_charge_constraint] = 1
            A[molecular_charge_constraint, :molecular_charge_constraint] = 1
            B[molecular_charge_constraint] = molecular_charge
            edges_f += 1

        # Add constraints to matrices A and B
        for i in range(1, n_constraint[mol]+1):
            if mol == 0 and i == n_constraint[mol]:
                # To account for the total charge constraints in the first molecule
                break
            B[edges_f] = constrained_charges[mol][i-1]
            for k in constrained_indices[mol][i-1]:
                if k > 0:
                    A[edges_f, edges_i+k-1] = 1
                    A[edges_i+k-1, edges_f] = 1
                else:
                    A[edges_f, edges_i-k-1] = -1
                    A[edges_i-k-1, edges_f] = -1
            edges_f += 1
        edges_i = edges_f
    indices = np.array(indices).flatten()

    # Add intermolecular constraints to A and B
    if n_mols > 1:
        for i in range(n_sys_constraint):
            B[edges_f] = con_charges_sys[i]
            for k in range(len(con_indices_sys[i])):
                for l in con_indices_sys[i][k]:
                    index = con_mol_sys[i][k]-1
                    index = int(np.sum(n_atoms[:index]) + np.sum(n_constraint[:index]))
                    if l > 0:
                        A[edges_f, index+l-1] = 1
                        A[index+l-1, edges_f] = 1
                    else:
                        A[edges_f, index-l-1] = -1
                        A[index-l-1, edges_f] = -1
            edges_f += 1

    labelf.append('ESP')
    q = esp_solve(A, B)
    qf.append(q[indices])
    if not rest:
        return qf, labelf, ''
    else:
        ihfree, resp_a, resp_b = [], [], []
        for mol in range(n_mols):
            ihfree.append(options[mol]['IHFREE'])
            resp_a.append(options[mol]['RESP_A'])
            resp_b.append(options[mol]['RESP_B'])
        toler = options[0]['TOLER']
        maxit = options[0]['MAX_IT']
        # Restrained ESP
        labelf.append('RESP')
        q, note = iterate(q, A, B, resp_a, resp_b, ihfree, symbols, toler, maxit,
                    n_atoms, n_constraint, indices)
        qf.append(q)
        return qf, labelf, note
