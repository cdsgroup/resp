from qcdb.moptions.read_options2 import RottenOption

def load_defaults(options):
    """Load RESP code default options and validate user-defined ones
    """

    options.add('RESP', RottenOption(
        keyword='ESP',
        default=[''],
        validator=lambda x: [str(i) for i in x] if isinstance(x, (list, tuple)) else [str(x)],
        glossary='Path to a file containing electrostatic potential values'))

    
    options.add('RESP', RottenOption(
        keyword='GRID',
        default=[''],
        validator=lambda x: [str(i) for i in x] if isinstance(x, (list, tuple)) else [str(x)],
        glossary='Path to a file containing grid points'))


    options.add('RESP', RottenOption(
        keyword='N_VDW_LAYERS',
        default=[4],
        validator=lambda x: [int(i) for i in x] if isinstance(x, (list, tuple)) else [int(x)],
        glossary='Number of van der Waals layers'))


    options.add('RESP', RottenOption(
        keyword='VDW_SCALE_FACTOR',
        default=[1.4],
        validator=lambda x: [float(i) for i in x] if isinstance(x, (list, tuple)) else [float(x)],
        glossary='Scale factor to mutiply van der Waals radii'))


    options.add('RESP', RottenOption(
        keyword='VDW_INCREMENT',
        default=[0.2],
        validator=lambda x: [float(i) for i in x] if isinstance(x, (list, tuple)) else [float(x)],
        glossary='Increment in Angstroms between van der Waals layers'))


    options.add('RESP', RottenOption(
        keyword='VDW_POINT_DENSITY',
        default=[1.0],
        validator=lambda x: [float(i) for i in x] if isinstance(x, (list, tuple)) else [float(x)],
        glossary='Point density in 1/Angstroms^2 on the molecular surface'))


    options.add('RESP', RottenOption(
        keyword='RADIUS',
        default=[{}],
        validator=lambda x: [dict(i) for i in x] if isinstance(x, (list, tuple)) else [dict(x)],
        glossary='User-defined van der Waals radii'))


    options.add('RESP', RottenOption(
        keyword='WEIGHT',
        default=[1.0],
        validator=lambda x: [float(i) for i in x] if isinstance(x, (list, tuple)) else [float(x)],
        glossary='Weight of molecule in multi-conformational fit'))


    options.add('RESP', RottenOption(
        keyword='RESTRAINT',
        default=True,
        validator=lambda x: bool(x),
        glossary='Apply hyperbolic restraint'))


    options.add('RESP', RottenOption(
        keyword='RESP_A',
        default=[0.0005],
        validator=lambda x: [float(i) for i in x] if isinstance(x, (list, tuple)) else [float(x)],
        glossary='The restraint a factor'))


    options.add('RESP', RottenOption(
        keyword='RESP_B',
        default=[0.1],
        validator=lambda x: [float(i) for i in x] if isinstance(x, (list, tuple)) else [float(x)],
        glossary='The restraint b factor'))


    options.add('RESP', RottenOption(
        keyword='IHFREE',
        default=[True],
        validator=lambda x: [bool(i) for i in x] if isinstance(x, (list, tuple)) else [bool(x)],
        glossary='Do not restrain hydrogens'))


    options.add('RESP', RottenOption(
        keyword='TOLER',
        default=1e-5,
        validator=lambda x: float(x),
        glossary='Tolerance of the charges in the fit'))


    options.add('RESP', RottenOption(
        keyword='MAX_IT',
        default=25,
        validator=lambda x: int(x),
        glossary='Maximum number of iterations'))


    options.add('RESP', RottenOption(
        keyword='QM_PACKAGE',
        default=[''],
        validator=lambda x: [str(i).upper() for i in x] if isinstance(x, (list, tuple)) else [str(x).upper()],
        glossary='Quantum chemistry software to use'))


    options.add('RESP', RottenOption(
        keyword='CONSTRAINT_CHARGE',
        default=[[]],
        validator=lambda x: [list(i) for i in x] if isinstance(x, (list, tuple)) else [list(x)],
        glossary='Constrain a group of atoms to a given charge'))


    options.add('RESP', RottenOption(
        keyword='CONSTRAINT_GROUP',
        default=[[]],
        validator=lambda x: [list(i) for i in x] if isinstance(x, (list, tuple)) else [list(x)],
        glossary='Make all atoms in a group have the same charge'))


    options.add('RESP', RottenOption(
        keyword='CONSTRAINT_EQUAL',
        default=[[]],
        validator=lambda x: [list(i) for i in x] if isinstance(x, (list, tuple)) else [list(x)],
        glossary='Set charges in two groups equal, one by one'))
