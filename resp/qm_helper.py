import qcdb
import numpy as np
import subprocess
import os

class gamess(object):
    def parse_gamess_esp(output):
        """ A fucntion to parse GAMESS output files and extract
            electrostatic potential. It can also parse concatenated
            GAMESS files.
            input:
                  output: path to output file
            output:
                  esp: numpy array of electrostatic potentials
            Note: electrostatic potential is available in grid_esp.dat
        """
        f = open(output, 'r').readlines()
        i = 0
        esp = []
        while i < len(f):
            if 'ELECTROSTATIC POTENTIAL' in f[i]:
                while True:
                    try:
                        line = f[i].split(' ')
                        esp.append(float(line[-1].split('\n')[0]))
                        i += 1
                    except:
                        if 'END OF PROPERTY EVALUATION' in f[i]:
                            break
                        i += 1
                        continue
                    
            i += 1
        esp = np.array(esp)
        np.savetxt('grid_esp.dat', esp, fmt='%14.10f')
        return esp

    def gamess_esp(molecule, method, basis):
        """ Computation of QM electrostatic potential with GAMESS
            input:
                  molecule: an instance of the QCDB Molecule class
            output:
                  esp: numpy array of electrostatic potentials
            Note: Only support RHF/6-31G*. To run other methods
                  manually modify the input files and then use
                  parse_gamess_esp function to extract output.
        """
        grid_complete = open('grid.dat', 'r').readlines()
        count = 0
        arrays = molecule.to_arrays()
        symbols, atomic_number, coordinates = arrays[2], arrays[3], np.array(arrays[0])
        if molecule.units() == 'Angstrom':
            coordinates *= qcdb.physconst.psi_bohr2angstroms

        for i in range(len(grid_complete)//99+1):
            try:
                grid = grid_complete[count:count+99]
            except:
                grid = grid_complete[count:]
            count += 99
            f = open('%i_' %i + molecule.name() + '_gamess.inp', 'w')
            input_file = """ $CONTRL ICHARG=%i MULT=%i RUNTYP=ENERGY UNITS=%s %s $END
 $BASIS  %s                        $END
 $ELPOT  IEPOT=1 WHERE=POINTS OUTPUT=BOTH                    $END
 $DATA
 MOLECULE
 C1
""" %(molecule.molecular_charge(), molecule.multiplicity(), molecule.units()[:4].upper(), method, basis)
            f.write(input_file)
            for j in range(molecule.natom()):
                f.write(' ' + symbols[j] + '   ' + str(atomic_number[j]) + '   ')
                for k in range(3):
                    f.write('%16.9f' %(coordinates[j, k]) + '   ')
                f.write('\n')
            f.write(' $END\n')
            f.write(""" $POINTS
 %s %i
""" %(molecule.units()[:4].upper(), len(grid)))
            for j in grid:
                f.write(j)
            f.write('$END\n\n')

            f.close()
            os.system('rungms '+'%i_'%i+molecule.name()+'_gamess.inp > '+ \
                  '%i_'%i+molecule.name()+'_gamess.out')

        os.system('cat *_%s_*gamess.inp > %s_gamess.inp' %(molecule.name(), molecule.name()))
        os.system('cat *_%s_*gamess.out > %s_gamess.out' %(molecule.name(), molecule.name()))
        os.system('rm *_%s_*gamess.inp *_%s_*gamess.out' %(molecule.name(), molecule.name()))
        esp = gamess.parse_gamess_esp('%s_gamess.out' %molecule.name())

        return esp

class qchem(object):
    def qchem_esp(molecule, method, basis):
        """ Computation of QM electrostatic potential with Q-CHEM
            input:
                  molecule: an instance of the QCDB Molecule class
            output:
                  esp: numpy array of electrostatic potentials
        """
        grid = np.loadtxt('grid.dat')
        arrays = molecule.to_arrays()
        symbols, coordinates = arrays[2], np.array(arrays[0])
        coordinates *= qcdb.physconst.psi_bohr2angstroms
        if molecule.units() == 'Angstrom':
            grid /= qcdb.physconst.psi_bohr2angstroms
        input_file = """$molecule
    %i %i
    """ %(molecule.molecular_charge(), molecule.multiplicity())
        for i in range(molecule.natom()):
            input_file += ' ' + symbols[i] + '   '
            for j in range(3):
                input_file += '%16.9f' %(coordinates[i, j]) + '   '
            input_file += '\n'
        input_file += '$end\n\n'
        input_file += """$rem
       METHOD     %s
       BASIS      %s
       IGDESP     %i
       IANLTY     200
    $end
     
    $plots
    Compute ESP
    """ %(method, basis, len(grid))
        input_file += '1 0 0\n1 0 0\n1 0 0\n0  0  0  0\n$end\n\n$grid\n'

        for i in range(len(grid)):
            input_file += '%16.9f%16.9f%16.9f\n' %(grid[i, 0], grid[i, 1], grid[i, 2])
        input_file += '$end\n'
        with open(molecule.name() + '_qchem.in', 'w') as f:
            f.write(input_file)
        subprocess.run(['qchem', molecule.name() + '_qchem.in', molecule.name() + '_qchem.out'])
        esp = np.loadtxt('fort.9', skiprows=4, usecols=3)
        np.savetxt('grid_esp.dat', esp, fmt='%14.10f') 

        return esp
