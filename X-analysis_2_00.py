# !/usr/bin/env python
# This program approximates the average electron binding energy (X) 
# from the eigenvalues of molecular orbitals 
# An "Experimental Quantum Chemistry" energy decomposition analysis is also performed. 
# Both closed and open-shell wavefunctions are OK
# Requirements: numpy, cclib, >Python 3.0. 

#USAGE:
# use "python X-analysis.py -h" for the help and the list of options

# Written by Martin Rahm, Chalmers University of Technology 
# Modified and updated by Francesco Sessa, Chalmers University of Technology 
# MIT License
# Copyright (c) [2018] [Martin Rahm]

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

# Citation is greatfully appreciated: 
# M. Rahm, R. Hoffmann, J. Am. Chem. Soc, 137, 10282-10291, 2015
# M. Rahm, R. Hoffmann, J. Am. Chem. Soc, 138, 3731-3744, 2016  

# and for cclib:
# N. M. O'Boyle, A. L. Tenderholt, K. M. Langner, J. Comput. Chem. 29, 839-845, 2008

#The most recent version of the code is accessible through a link at https://rahmlab.com/energy-decomposition-analysis/ 

## defining some functions
def print_header():
    print('>>{:s}<<'.format('-'*70))
    print('>>{:^70}<<'.format('EQC ENERGY DECOMPOSITION ANALYSIS'))
    print('>>{:^70}<<'.format('https://rahmlab.com/energy-decomposition-analysis/'))
    print('>>{:^70}<<'.format('Version: 2.04'))
    print('>>{:^70}<<'.format('Last revision: Jul 23, 2020'))
    print('>>{:^70}<<'.format(''))
    print('>>{:^70}<<'.format('All energy value are given in electron volt (eV)'))
    print('>>{:^70}<<'.format(''))
    print('>>{:^70}<<'.format('Citation is gratefully appreciated:'))
    print('>>{:^70}<<'.format('M. Rahm, R. Hoffmann, J. Am. Chem. Soc, 137, 10282-10291, 2015'))
    print('>>{:^70}<<'.format('M. Rahm, R. Hoffmann, J. Am. Chem. Soc, 138, 3731-3744, 2016'))
    print('>>{:^70}<<'.format(''))
    print('>>{:s}<<'.format('-'*70))
    print('')


def custom_round(number,prec=0): 
    # improved rounding to handle floating point errors. 
    # 'prec' is the number of decimal digits
    
    true_prec = 10**prec
    num_out = number * true_prec
    num_out = (num_out)//1 + (((num_out*2)%2)//1)
    num_out /= true_prec
    return num_out

def guess_atom_core(atom_num): 
    # returns a guess of the number of core electrons in an atom, based on the atomic number. 
    
    if (atom_num<=2):    # no core orbital for H and He
        return 2*0 
    elif (atom_num<=10): # core=1s for Li to Ne
        return 2*1
    elif (atom_num<=18): # core=1s,2s,2p for Na to Ar
        return 2*5
    elif (atom_num<=30): # core=1s,2s,2p,3s,3p for K to Zn
        return 2*9
    elif (atom_num<=36): # core=1s,2s,2p,3s,3p,3d for Ga to Kr
        return 2*14
    elif (atom_num<=48): # core=1s,2s,2p,3s,3p,3d,4s,4p for Rb to Cd
        return 2*18
    elif (atom_num<=54): # core=1s,2s,2p,3s,3p,3d,4s,4p,4d for In to Xe
        return 2*23
    elif (atom_num<=71): # core=1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p for Cs to Lu
        return 2*27
    elif (atom_num<=80): # core=1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f for Hf to Hg
        return 2*34
    elif (atom_num<=86): # core=1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d for Tl to Rn
        return 2*39
    elif (atom_num<=103):
        return 2*43
    elif (atom_num<=112):
        return 2*50
    elif (atom_num<=118):
        return 2*55
    else:
        print('WARNING: ATOM NUMBER ABOVE 118 DETECTED!!!')
        sys.exit()

def guess_core (data): 
    # returns a guess of the number of core orbitals in the system. 
    # 'data' is a cclib-parsed output of a QM calculation. 
    
    core_orb_theo=0
    max_occ_orb= -data.charge
    for x in range(data.natom):
        core_orb_theo += guess_atom_core(data.atomnos[x])//2
        max_occ_orb += data.atomnos[x]
    max_occ_orb = max_occ_orb//2 + data.mult//2 
    core_orb = core_orb_theo - (max_occ_orb - (data.homos[0] + 1)) #removes core orbitals in case of pseudo-potentials
    ## In case of partial (non integer) occupation numbers for alpha spin-orbitals, the algorithm guesses too many core orbitals. A correction is introduced here
    sum_part = 0.0
    if hasattr(data,'moocc'):
        for x in data.moocc[0]:
            if x>0.01 and x<0.99:
                sum_part += 1.0 - x
    core_orb = core_orb - int(custom_round(sum_part))
    ##
    return (core_orb)

def nuclear_repulsion_energy(data): 
    # Returns the value of the nuclear repulsion in a system based on the optimized geometry.
    # 'data' is a cclib-parsed output of a QM calculation.
    
    nre = 0.0
    for i in range(data.natom):
        ri = data.atomcoords[-1][i]
        zi = data.atomnos[i]
        for j in range(i+1, data.natom):
            rj = data.atomcoords[-1][j]
            zj = data.atomnos[j]
            d = numpy.linalg.norm(ri-rj)
            nre += zi*zj/d
    nre = cclib.parser.utils.convertor(nre,"bohr", "Angstrom")
    nre = cclib.parser.utils.convertor(nre,"hartree","eV")
    return nre

def Xanalysis(data,core_orbs): 
    # perform the EQC analysis on a single QM calculation. 
    # 'data' is a cclib-parsed output of a QM calculation. 
    # 'core_orbs' is the number of core orbitals of the system (if 0 the analysis makes no core/valence distinction). 
    # Returns a dictionary of the following properties:
    # 'openshell' : boolean. True if 'data' is an openshell calculation. 
    # 'HOMO' and 'LUMO' : floats. eigenvalues of HOMO and LUMO orbitals in eV.
    # 'X_core', 'X_val' and 'X_tot' : floats. eigenvalues sum for core, valence and all orbitals, respectively.
    # 'core_electrons', 'valence_electrons' and 'electrons' : integers. core, valence and total number of electrons in system.
    # 'energy' : float. total energy of the system.
    # 'VNN', 'multi_elec' and 'VNN_EEE' : floats. N-N and e-e interactions, and their difference.
    
    results = {}
    
    openshell = (len(data.homos)==2) # check if it is an openshell calculation
    results['openshell'] = openshell
    
    ## evaluating HOMO and LUMO

    HOMO = data.moenergies[0][data.homos[0]]
    LUMO = data.moenergies[0][data.homos[0] + 1]
    if openshell:
        if (data.moenergies[0][data.homos[0]] <  data.moenergies[1][data.homos[1]]) and (data.homos[1]>=0):
            HOMO = data.moenergies[1][data.homos[1]]
        if data.moenergies[0][data.homos[0] + 1] >  data.moenergies[1][data.homos[1] + 1]:
            LUMO = data.moenergies[1][data.homos[1] + 1]
            
    (results['HOMO'],results['LUMO']) = (HOMO,LUMO)
    

    E = data.scfenergies[-1] # total energies is the last SCF energy
    
    ## Calculating n*X-bar (with valence/core separation as well, if requested)
    
    sum_val=0.0
    sum_core=0.0
    sum_tot=0.0
        
    for spin in range(len(data.homos)):
        for m in range(core_orbs):
            sum_core += data.moenergies[spin][m] * data.moocc[spin][m]
            sum_tot += data.moenergies[spin][m] * data.moocc[spin][m]
        for m in range(core_orbs,data.homos[spin] + 1):
            sum_val += data.moenergies[spin][m] * data.moocc[spin][m]
            sum_tot += data.moenergies[spin][m] * data.moocc[spin][m]
            
    ##


    ## Calculating multi-electron term (Eee)
    VNN_EEE = E - sum_tot
    VNN = nuclear_repulsion_energy(data)
    multi_elec = VNN - VNN_EEE
    
    ##
    
    ## returning results as a dictionary
    
    (results['X_tot'] ,results['energy'] ,results['VNN'] ,results['multi_elec'] ,results['VNN_EEE'] ,
     results['electrons']) = (sum_tot,E, VNN, multi_elec, VNN_EEE,data.nelectrons)
    
    if core_orbs > 0:
        (results['X_core'] ,results['X_val'] ,results['core_electrons'] ,results['valence_electrons'] ) = (sum_core, sum_val,2*core_orbs,data.nelectrons-2*core_orbs)


    return(results)
#

def print_Xanalysis(data,core,Xan_results):
    # prints out the results of the EQC analysis performed on a single QM calculation. 
    # 'data' is a cclib-parsed output of a QM calculation. 
    # 'core' is a boolean flag. If TRUE, valence and core orbitals are treated separately.
    # 'Xan_results' is the dictionary returned by the Xanalysis function.
    
    ## Printing initial metadata
    
    if Xan_results['openshell']:
        print('This wavefunction is open-shell')
        print("There are {:d} alpha and {:d} beta occupied spin-orbitals".format(data.homos[0] + 1,data.homos[1] + 1) )
        if data.metadata['package'] == 'ADF':
            print("Calculations performed with ADF might end up in partial occupation numbers")
    else:
        print('This wavefunction is closed-shell')
        print("There are {:d} occupied orbitals".format(data.homos[0] + 1))

    print("System contains {:d} electrons".format(Xan_results['electrons']))
    
    ##

    ## Listing Orbital Energies
    
    print('')
    print("----ORBITALS--------------------------")
    
    if Xan_results['openshell']:
        print("Energies of Occupied Alpha Orbitals:")
        for m in range(data.homos[0] + 1):
            if core:
                if (m<(Xan_results['core_electrons']//2)):
                    label='core'
                else:
                    label='valence'
            else:
                label=''
            print('{:d}   {:^7}{:10.2f} eV'.format(m,label,custom_round(data.moenergies[0][m],2)) )
        print("Energies of Occupied Beta Orbitals:")
        for m in range(data.homos[1] + 1):
            if core:
                if (m<(Xan_results['core_electrons']//2)):
                    label='core'
                else:
                    label='valence'
            else:
                label=''
            print('{:d}   {:^7}{:10.2f} eV'.format(m,label,custom_round(data.moenergies[1][m],2)) )
            
    else:
        print("Energies of Occupied Orbitals:")
        for m in range(data.homos[0] + 1):
            if core:
                if (m<(Xan_results['core_electrons']//2)):
                    label='core'
                else:
                    label='valence'
            else:
                label=''
            print('{:d}   {:^7}{:10.2f} eV'.format(m,label,custom_round(data.moenergies[0][m],2)) )
            
    print('')
            
    print('LUMO:      {:10.2f} eV'.format(Xan_results['LUMO']))
    print('HOMO:      {:10.2f} eV'.format(Xan_results['HOMO']))
    print('GAP:       {:10.2f} eV'.format(Xan_results['HOMO'] - Xan_results['LUMO']))
    print('')    
            
    ## Printing X_Analysis totals
    
    print("----TOTAL ENERGIES-------------------")

    print("E = n*X-bar + Vnn - Eee")
    print('')
    print("Total SCF energy (E):                                   {:10.2f} eV".format(custom_round(Xan_results['energy'],2))) 
    print("Total energy of occupied orbitals (n*X-bar):            {:10.2f} eV".format(custom_round(Xan_results['X_tot'],2)))
    if core:
        print("Total energy of valence electrons (n*X-bar_valence):    {:10.2f} eV".format(custom_round(Xan_results['X_val'],2)))
        print("Total energy of core orbitals (n*X-bar_core):           {:10.2f} eV".format(custom_round(Xan_results['X_core'],2)))
    print("Nuclear repulsion (Vnn):                                {:10.2f} eV".format(custom_round(Xan_results['VNN'],2)))
    print("Multi-electron term (Eee):                              {:10.2f} eV".format(custom_round(Xan_results['multi_elec'],2)))
    print("Vnn - Eee:                                              {:10.2f} eV".format(custom_round(Xan_results['VNN_EEE'],2)))
    print('')
    
    ##
    
    ## Printing X_Analysis (per electron)
    
    print("----TOTAL ENERGIES--(PER ELECTRON)----")
    
    print("E/n = X-bar + (Vnn - Eee)/n")
    print('')
    print("E/n:           {:40} {:10.2f} eV".format('',custom_round(Xan_results['energy']/Xan_results['electrons'],2))) 
    print("X-bar:         {:40} {:10.2f} eV".format('',custom_round(Xan_results['X_tot']/Xan_results['electrons'],2)))
    if core:
        print("X-bar_valence: {:40} {:10.2f} eV".format('',custom_round(Xan_results['X_val']/Xan_results['valence_electrons'],2)))
        print("X-bar_core:    {:40} {:10.2f} eV".format('',custom_round(Xan_results['X_core']/Xan_results['core_electrons'],2)))
    print("Vnn/n:         {:40} {:10.2f} eV".format('',custom_round(Xan_results['VNN']/Xan_results['electrons'],2)))
    print("Eee/n:         {:40} {:10.2f} eV".format('',custom_round(Xan_results['multi_elec']/Xan_results['electrons'],2)))
    print("(Vnn - Eee)/n: {:40} {:10.2f} eV".format('',custom_round(Xan_results['VNN_EEE']/Xan_results['electrons'],2)))
    print('')
    
    print("-------MULTIELECTRON ANALYSIS----------")
    print("100*(|Eee/E|):      {:6.2f}%".format(abs(Xan_results['multi_elec']/Xan_results['energy'])*100))
    print('')
    print('')
    print('')
    
    return

def exp_values_Xanalysis(atomic_number=1): # only  1 for hydrogen so far is tabulated. atomic_number=0 is a special case, returns a dummy Xan_result for an isolated electron (used for ionizations)
    
    results = {}
    if atomic_number==0:
        results['openshell'] = True
        (results['HOMO'],results['LUMO']) = (0,0)
        (results['X_tot'] ,results['energy'] ,results['VNN'] ,results['multi_elec'] ,results['VNN_EEE'] , results['electrons']) = (0, 0, 0, 0, 0, 1)
        (results['X_core'] ,results['X_val'] ,results['core_electrons'] ,results['valence_electrons'] ) = (0.0,0.0,0,1)
    elif atomic_number==1:
        results['openshell'] = True
        (results['HOMO'],results['LUMO']) = (-13.6,'')
        (results['X_tot'] ,results['energy'] ,results['VNN'] ,results['multi_elec'] ,results['VNN_EEE'] , results['electrons']) = (-13.6,-13.6, 0, 0, 0,1)
        (results['X_core'] ,results['X_val'] ,results['core_electrons'] ,results['valence_electrons'] ) = (0.0,-13.6,0,1)
    
    
    return results

def reaction_Xanalysis(stoich_reactants,Xan_reactants,stoich_products,Xan_products,core=False):
    
    elec_reac = 0
    elec_prod = 0
    for i in range(len(stoich_products)):
        elec_prod += stoich_products[i]*Xan_products[i]['electrons']
    for i in range(len(stoich_reactants)):
        elec_reac += stoich_reactants[i]*Xan_reactants[i]['electrons']
    nelec = elec_prod
    
    ndX = 0
    for i in range(len(stoich_products)):
        ndX += stoich_products[i]*Xan_products[i]['X_tot']
    for i in range(len(stoich_reactants)):
        ndX -= stoich_reactants[i]*Xan_reactants[i]['X_tot']
        
    if core:
        dX_reac = 0
        dX_prod = 0
        elec_reac = 0
        elec_prod = 0
        for i in range(len(stoich_products)):
            dX_prod += stoich_products[i]*Xan_products[i]['X_core']
            elec_prod += stoich_products[i]*Xan_products[i]['core_electrons']
        dX_prod /= elec_prod
        for i in range(len(stoich_reactants)):
            dX_reac += stoich_reactants[i]*Xan_reactants[i]['X_core']
            elec_reac += stoich_reactants[i]*Xan_reactants[i]['core_electrons']
        dX_reac /= elec_reac
        dXcore = dX_prod - dX_reac
        
        dX_reac = 0
        dX_prod = 0
        elec_reac = 0
        elec_prod = 0
        for i in range(len(stoich_products)):
            dX_prod += stoich_products[i]*Xan_products[i]['X_val']
            elec_prod += stoich_products[i]*Xan_products[i]['valence_electrons']
        dX_prod /= elec_prod
        for i in range(len(stoich_reactants)):
            dX_reac += stoich_reactants[i]*Xan_reactants[i]['X_val']
            elec_reac += stoich_reactants[i]*Xan_reactants[i]['valence_electrons']
        dX_reac /= elec_reac
        dXval = dX_prod - dX_reac
        
    
    dE = 0
    for i in range(len(stoich_products)):
        dE += stoich_products[i]*Xan_products[i]['energy']
    for i in range(len(stoich_reactants)):
        dE -= stoich_reactants[i]*Xan_reactants[i]['energy']
        
    dVnn = 0
    for i in range(len(stoich_products)):
        dVnn += stoich_products[i]*Xan_products[i]['VNN']
    for i in range(len(stoich_reactants)):
        dVnn -= stoich_reactants[i]*Xan_reactants[i]['VNN']
        
    dEee = 0
    for i in range(len(stoich_products)):
        dEee += stoich_products[i]*Xan_products[i]['multi_elec']
    for i in range(len(stoich_reactants)):
        dEee -= stoich_reactants[i]*Xan_reactants[i]['multi_elec']
        
    dVEnn = 0
    for i in range(len(stoich_products)):
        dVEnn += stoich_products[i]*Xan_products[i]['VNN_EEE']
    for i in range(len(stoich_reactants)):
        dVEnn -= stoich_reactants[i]*Xan_reactants[i]['VNN_EEE']
        
    Q = ((2*ndX)/dE) - 1
    
    Covalency = 100 * ( (abs(ndX)-ndX)/(2*(abs(dVEnn)-ndX)) + ((abs(ndX)+ndX)/(2*((dVEnn)-ndX))) + 
                       ((abs(ndX)+ndX)/(2*ndX)) )
    
    
    temp_1 = 0
    temp_2 = 0
    temp_3 = 0
    temp_4 = 0
    for i in range(len(stoich_products)):
        temp_1 += stoich_products[i]*Xan_products[i]['multi_elec']
        temp_2 += stoich_products[i]*Xan_products[i]['energy']
    for i in range(len(stoich_reactants)):
        temp_3 += stoich_reactants[i]*Xan_reactants[i]['multi_elec']
        temp_4 += stoich_reactants[i]*Xan_reactants[i]['energy']
    d_EEE_E = abs(temp_1/temp_2) - abs(temp_3/temp_4)

    results = {}
    (results['X_tot'] ,results['energy'] ,results['electrons'] ) = (ndX,dE,nelec)
    (results['VNN'] ,results['multi_elec'] ,results['VNN_EEE'] ) = (dVnn,dEee,dVEnn)
    (results['Q'] ,results['Covalency'] ,results['EEE_E'] ) = (Q,Covalency,d_EEE_E)
    
    if core:
        (results['X_core'] ,results['X_val'] ) = (dXcore,dXval)


    return(results)

def print_reaction_Xanalysis(Xan_results,core=False):
    
    print("E/n = X-bar + (Vnn - Eee)/n")    
    print("-------------------------------------------------")
    print("Electrons (n):           {:10d}".format(Xan_results['electrons']))
    print("Delta(E):                {:10.2f} eV".format(custom_round(Xan_results['energy'],2)))
    print("Delta(n*X-bar):          {:10.2f} eV".format(custom_round(Xan_results['X_tot'],2)))
    print("Delta(Vnn - Eee):        {:10.2f} eV".format(custom_round(Xan_results['VNN_EEE'],2))  )
    print('')
    print("Delta(E/n):              {:10.2f} eV".format(custom_round(Xan_results['energy']/Xan_results['electrons'],2)))
    print("Delta(X-bar):            {:10.2f} eV".format(custom_round(Xan_results['X_tot']/Xan_results['electrons'],2)))
    if core:
        print("Delta(X-bar)_core:       {:10.2f} eV".format(custom_round(Xan_results['X_core'],2)))
        print("Delta(X-bar)_valence:    {:10.2f} eV".format(custom_round(Xan_results['X_val'],2)))
    print("Delta(Vnn/n):            {:10.2f} eV".format(custom_round(Xan_results['VNN']/Xan_results['electrons'],2)))
    print("Delta(Eee/n):            {:10.2f} eV".format(custom_round(Xan_results['multi_elec']/Xan_results['electrons'],2)))
    print("Delta(Vnn-Eee)/n:        {:10.2f} eV".format(custom_round(Xan_results['VNN_EEE']/Xan_results['electrons'],2)))
    print('')
    print("Q:                       {:10.2f}".format(custom_round(Xan_results['Q'],2)) )
    print("Delta(|Eee/E|)           {:10.2f}%".format(custom_round(100*( Xan_results['EEE_E'] )),2))
    
    
    print("Covalency Index:         {:10.2f}%".format(custom_round(Xan_results['Covalency'],2)))
    print("Ionicity Index:          {:10.2f}%".format(custom_round(100-Xan_results['Covalency'],2)))
    
    return

def write_xyz_cclib(data,output_name,step=-1,append=False): 
    # write a .xyz file of a system from a QM calculation. 
    # 'data' is a cclib-parsed output of a QM calculation. 
    
    if append:
        flag='a'
    else:
        flag='w'
    output_file = open(output_name,flag)
    output_file.write('{}\n'.format(data.natom))
    output_file.write('\n')
    for atom in range(data.natom):
        output_file.write('{:2d} '.format(data.atomnos[atom]))
        for i in range(3):
            output_file.write('{:8.4f} '.format(data.atomcoords[step][atom][i]))
        output_file.write('\n')
        
    return

## end of functions

import cclib
from cclib.io import ccopen 
import sys
import os 
import numpy
import argparse

## building the argument parser
class SmartFormatter(argparse.HelpFormatter):

    #def _split_lines(self, text, width):
    #    if text.startswith('R|'):
    #        return text[2:].splitlines() # this is the RawTextHelpFormatter._split_lines
    #    return argparse.HelpFormatter._split_lines(self, text, width)
    def _split_lines(self, text, width):
        return text.splitlines() # this is the RawTextHelpFormatter._split_lines

arg_parser = argparse.ArgumentParser(usage=('%(prog)s [options]'),description=('This program approximates the average electron binding energy (X) from the eigenvalues of molecular orbitals. An "Experimental Quantum Chemistry" energy decomposition analysis is also performed.'),formatter_class=SmartFormatter)

arg_parser.add_argument("-i","--inputfile",required=True,nargs='*',help=("List of output files of QM calculations to be analysed.\n"
                                                                         " "))

arg_parser.add_argument("-m","--mode",type=int,choices=range(4),default=0,help=("Define how the EQC Analysis will be performed:\n"
                                                                                "(0) EQC analysis on each QM calculation provided.\n"
                                                                                "(1) EQC analysis for the transformation A -> B.\n"
                                                                                "      Example: isomerization or ionization\n"
                                                                                "      Requires 2 QM calculations (A and B) as input(-i flag).\n"
                                                                                "(2) EQC analysis for the reaction A + B -> C.\n"
                                                                                "      Requires 3 QM calculations (A, B and C) as input(-i flag).\n"
                                                                                "(3) EQC analysis on a general chemical reaction.\n"
                                                                                "      This is a 'free format reaction' mode.\n"
                                                                                "      Stoichiometry flag (-s) is mandatory.\n"
                                                                                "      Positive stoichiometric coefficients mark species as reactants,\n" 
                                                                                "      negative stoichiometric coefficients mark species as products.\n"
                                                                                "  "))

arg_parser.add_argument("-s","--stoich",type=int,nargs='*',default=[],help=("Specify a list of stoichiometric coefficients (integers) in any reaction mode.\n"
                                                                           "  Each element of the list is associated with the respective QM calculation in the input list.\n"
                                                                           "      Example: if 2 calculation are provided (e.g. "-i A.log B.log"),\n"
                                                                           "      two coefficients must be provided (e.g. "-s 1 1").\n"
                                                                           "  Mandatory if a 'free format reaction' mode is selected.\n"
                                                                           "    In 'free format', positive coefficients mark species as reactants,\n"
                                                                           "    negative coefficients mark species as products.\n"
                                                                           " "))

arg_parser.add_argument("-c","--core",action='store_true',help=("Include a separate treatment for valence and core electrons.\n"
                                                                "By default, the number of core electrons is guessed.\n"
                                                                "The guess is based on atomic numbers, charge and multiplicity.\n"
                                                                " "))

arg_parser.add_argument("-cm","--core-manual",action='store_true',help=("This flag allows for manual input of the number of core electrons.\n"
                                                                        "Only even numbers are supported.\n"
                                                                        "This flag overrides the -c flag.\n"
                                                                        " "))

arg_parser.add_argument("-o","--outputfile",default=argparse.SUPPRESS,help=("Path and name of the file to store the output of the analysis.\n"
                                                                            "If omitted the output is printed on stardard output.\n"
                                                                            " "))

arg_parser.add_argument("-v","--verbose",action='store_true',help=("Print single file EQC analyses in reaction mode.\n"
                                                                   " "))

##

args = arg_parser.parse_args() # parsing the command option

## safety control for reaction mode:
if args.core_manual:
    args.core=True
if args.mode == 1: # A -> B
    if len(args.inputfile) < 2:
        print('WARNING! mode (A -> B) requires 2 input files. Aborting.')
        sys.exit()
    elif len(args.inputfile) > 2:
        print('WARNING! mode (A -> B) requires 2 input files. All files beyond the second will be ignored.')
elif args.mode == 2: # A + B -> C
    if len(args.inputfile) < 3:
        print('WARNING! mode (A + B -> C) requires 3 input files. Aborting.')
        sys.exit()
    elif len(args.inputfile) > 3:
        print('WARNING! mode (A + B -> C) requires 3 input files. All files beyond the third will be ignored.')
##

# Starting the main program
## printing header and preparing outputfile, if requested

print_header()
if args.core:
        print('splitting valence and core molecular orbitals')

if hasattr(args,'outputfile'):
    oldstdout = sys.stdout
    outputfile = open(args.outputfile,'w')
    sys.stdout = outputfile
    print_header()
    outputfile.flush()
##

## cycling on input files. Each one is parsed and appended in data_list
if hasattr(args,'outputfile'):
    sys.stdout = oldstdout

print('Number of input file detected: ',len(args.inputfile))

if hasattr(args,'outputfile'):
    sys.stdout = outputfile

input_path = []
input_names = []
input_root = []
data_list = [] # empty list. Input QM calculations are going to be cclib-parsed here.

for y in range (len(args.inputfile)): 
    
    # removing input file extension to get input root name. It will be used for naming output files
    input_path_split = args.inputfile[y].split('/')[:-1]
    input_names.append(args.inputfile[y].split('/')[-1])
    input_split = input_names[y].split(".")[:-1]
    if len(input_path_split)>0:
        input_path.append( input_path_split[0]+'/' )
        for z in input_path_split:
            input_path[y] += z+'/'
    else:
        input_path.append( '' )
    input_root.append( 'opt_'+input_split[0] )
    for z in input_split:
        input_root[y] += '.'+z
        
        
    if hasattr(args,'outputfile'):
        sys.stdout = oldstdout
        print('parsing file: ' + args.inputfile[y])
        sys.stdout = outputfile

    data_list.append( ccopen(args.inputfile[y]).parse() ) # parsing y-th input file using cclib
        
    if data_list[y].metadata['package'] == 'Gaussian': # if it is a Gaussian file, deep core levels are checked
        
        with open("{0}".format(args.inputfile[y]), 'r') as fo:
            lines = fo.readlines()
                
        line_count=0
        file_replace = False
        for line in lines:
            if "**********-" in line:
                file_replace = True
                print("WARNING!! Deep core levels are not printed correctly in the Gaussian output. These will be ignored")
                lines[line_count] = line.replace("**********-", "0.0000000 -")
            if "Error termination" in line: #avoids output that ends with an error
                print("WARNING!! Error termination in Gaussian output detected. Aborting execution")
                sys.exit()
            line_count += 1
            
        if file_replace: # if necessary, a modified version of the input is parsed to avoid errors
            temp_input = "TEMP_INPUT_FILE"
            print('replacing lines in ' + temp_input)
            with open(temp_input,'w') as temp_data:
                temp_data.writelines(lines)
            data_list[y] = ccopen(temp_input).parse() # replacing the parsed data with the new one, if necessary
            os.remove(temp_input)
            
    if not hasattr(data_list[y],'moocc'):
        if len(data_list[y].homos) == 1:
            data_list[y].moocc = [numpy.zeros(data_list[y].nbasis, "d")]
        else:
            data_list[y].moocc = [numpy.zeros(data_list[y].nbasis, "d"),numpy.zeros(data_list[y].nbasis, "d")]
        for spin in range(len(data_list[y].homos)):
            for m in range(data_list[y].homos[spin] + 1):
                data_list[y].moocc[spin][m] = 1.0
        if len(data_list[y].homos) == 1:
            data_list[y].moocc = numpy.multiply(data_list[y].moocc,2)
##

## the number of core orbitals is handled here, either by guessing, manual input or setting it to zero (if no option is provided)
y=0
ncore = []
for data in data_list:
    if args.core:
        if args.core_manual:
            ncore_temp=int(input('number of core electrons for: {} '.format(input_names[y])))
            while (ncore_temp//2 != ncore_temp/2):
                print('')
                print('>>{}<<'.format('-'*80))
                print('>>{:^80}<<'.format('ERROR: Odd number for core electrons detected. Only even number are supported!!'))
                print('>>{}<<'.format('-'*80))
                print('')
                ncore_temp=int(input('number of core electrons for: {} '.format(input_names[y])))
            ncore.append(ncore_temp//2)
        else:
            ncore.append(guess_core(data))
    else:
        ncore.append(0)
    y+=1
##

## handling the different calculation modes

### single file EQC analysis (multiple files are possible, each one analyzed separately)
if args.mode == 0: 
    
    
    y=0
    for data in data_list:
        if hasattr(args,'outputfile'):
            sys.stdout = oldstdout
            print('analyzing file: '+args.inputfile[y])
            sys.stdout = outputfile
        print('>> {} <<'.format(input_names[y]))
        results = Xanalysis(data,ncore[y])
        print_Xanalysis(data,args.core,results)
        if hasattr(args,'outputfile'):
            outputfile.flush()
        y+=1
###
        
### Analyze the reaction A -> B (for isomers or ionizations) or nA -> nB (homonuclear molecules) (2 files only required)
if args.mode == 1:
    
    if len(args.stoich) < 2:
        args.stoich = [1,-1]
    else:
        if args.stoich[0] < 0:
            args.stoich[0] = -args.stoich[0]
        if args.stoich[1] > 0:
            args.stoich[1] = -args.stoich[1]
    args.mode = 3
    
###

### Analyzing the reaction A + B -> C (3 files required)
if args.mode == 2:
    
    if len(args.stoich) < 3:
        args.stoich = [1,1,-1]
    else:
        if args.stoich[0] < 0:
            args.stoich[0] = -args.stoich[0]
        if args.stoich[1] < 0:
            args.stoich[1] = -args.stoich[1]
        if args.stoich[2] > 0:
            args.stoich[2] = -args.stoich[2]
    args.mode = 3
    
###

### Analyzing a generic reaction. This mode allows for reactions with any number of reactants ad products, which are provided by the stoichiometry flag (positive coefficients are for reactants, negative for products.)
if args.mode in (3,):
    
    stoich_react = []
    names_react = []
    stoich_prod = []
    names_prod = []
    num_atoms_react = 0
    hydrogens_react = 0
    elec_react = 0
    num_atoms_prod = 0
    hydrogens_prod = 0
    elec_prod = 0
    electron_balance = 0
    isolated_atoms_reactants = True
    isolated_atoms_products = True
    for i in range(len(args.stoich)):
        if args.stoich[i] > 0:
            stoich_react.append(args.stoich[i])
            names_react.append(input_names[i])
            isolated_atoms_reactants = isolated_atoms_reactants and ( data_list[i].natom == 1 )
            num_atoms_react += data_list[i].natom * args.stoich[i]
            hydrogens_react += numpy.count_nonzero(data_list[i].atomnos == 1) * args.stoich[i]
            elec_react += args.stoich[i] * data_list[i].nelectrons
            electron_balance += data_list[i].charge * args.stoich[i]
        elif args.stoich[i] < 0:
            names_prod.append(input_names[i])
            stoich_prod.append(-args.stoich[i])
            isolated_atoms_products = isolated_atoms_products and ( data_list[i].natom == 1 )
            num_atoms_prod += data_list[i].natom * (-args.stoich[i])
            hydrogens_prod += numpy.count_nonzero(data_list[i].atomnos == 1)  * (-args.stoich[i])
            elec_prod += (-args.stoich[i]) * data_list[i].nelectrons
            electron_balance += data_list[i].charge * args.stoich[i]
            
        
    if len(stoich_react)==0:
        print('')
        print('>>{}<<'.format('-'*80))
        print('>>{:^80}<<'.format('ERROR: NO REACTANT DETECTED!!'))
        print('>>{:^80}<<'.format('Please refer to the help (-h flag) for the usage of free format reaction mode '))
        print('>>{}<<'.format('-'*80))
        print('')
        sys.exit()
    if len(stoich_prod)==0:
        print('')
        print('>>{}<<'.format('-'*80))
        print('>>{:^80}<<'.format('ERROR: NO PRODUCT DETECTED!!'))
        print('>>{:^80}<<'.format('Please refer to the help (-h flag) for the usage of free format reaction mode '))
        print('>>{}<<'.format('-'*80))
        print('')
        sys.exit()
        
    
    atom_num_mismatch = False
    electron_mismatch = False
    add_experimental_hydrogen = 0
    if ( num_atoms_react != num_atoms_prod ):
        if ( num_atoms_react - hydrogens_react == num_atoms_prod - hydrogens_prod ):
            add_experimental_hydrogen = (hydrogens_prod - hydrogens_react)
        else:
            atom_num_mismatch = True
    if ( add_experimental_hydrogen > 0 ):
        stoich_react.append(add_experimental_hydrogen)
        names_react.append('H(experimental)')
    elif ( add_experimental_hydrogen < 0 ):
        stoich_react.append(-add_experimental_hydrogen)
        names_react.append('H(experimental)')
        
    if ( elec_react+add_experimental_hydrogen != elec_prod ):
        if electron_balance > 0:
            stoich_react.append(electron_balance)
            names_react.append('e')
        elif electron_balance < 0:
            stoich_prod.append(-electron_balance)
            names_prod.append('e')
        else:
            electron_mismatch = True
        
    print('>> REACTION: <<')
    for i in range(len(stoich_react)):
        print('{} * {} '.format(stoich_react[i],names_react[i]) ,end='')
        if i < len(stoich_react)-1:
            print('+ ',end='')
        else:
            print('-> ',end='')
    for i in range(len(stoich_prod)):
        print('{} * {} '.format(stoich_prod[i],names_prod[i]) ,end='')
        if i < len(stoich_prod)-1:
            print('+ ',end='')
        else:
            print('')
    print('')
    if hasattr(args,'outputfile'):
        sys.stdout = outputfile
    print('')
    
    if atom_num_mismatch:
        print('>>{}<<'.format('-'*80))
        print('>>{:^80}<<'.format('ERROR: MISMATCH IN NUMBER OF ATOMS!!'))
        print('>>{}<<'.format('-'*80))
        sys.exit()
    elif electron_mismatch:
        print('')
        print('>>{}<<'.format('-'*80))
        print('>>{:^80}<<'.format('WARNING: initial and final number of electrons are different,'))
        print('>>{:^80}<<'.format('with no change in the system charge'))
        print('>>{}<<'.format('-'*80))
        print('')
            
    
    if args.verbose:
        print('verbose mode: printing single file EQC analysis first')
        print('')
    y=0
    results=[]
    for data in data_list:
        if args.verbose:
            print('analyzing file: '+args.inputfile[y])
            print('')
        results.append( Xanalysis(data,ncore[y]))
        if args.verbose:
            print_Xanalysis(data,args.core,results[y])
        y+=1
        
    reactants = []
    products = []
    for i in range(len(args.stoich)):
        if args.stoich[i] > 0:
            reactants.append(results[i])
        elif args.stoich[i] < 0:
            products.append(results[i])
            
    if ( add_experimental_hydrogen > 0 ):
        reactants.append(exp_values_Xanalysis(1))
    elif ( add_experimental_hydrogen < 0 ):
        products.append(exp_values_Xanalysis(1))
        
    if ( electron_balance > 0 ):
        reactants.append(exp_values_Xanalysis(0))
    elif ( electron_balance < 0 ):
        products.append(exp_values_Xanalysis(0))
        
        
    if (hasattr(args,'outputfile') or args.verbose):
        print('>> REACTION: <<')
        for i in range(len(stoich_react)):
            print('{} * {} '.format(stoich_react[i],names_react[i]) ,end='')
            if i < len(stoich_react)-1:
                print('+ ',end='')
            else:
                print('-> ',end='')
        for i in range(len(stoich_prod)):
            print('{} * {} '.format(stoich_prod[i],names_prod[i]) ,end='')
            if i < len(stoich_prod)-1:
                print('+ ',end='')
            else:
                print('')
        print('')
        
    
    reaction_results = reaction_Xanalysis(stoich_react,reactants,stoich_prod,products,args.core)
    
    print_reaction_Xanalysis(reaction_results,args.core)
       
###

##

## closing output file, if the option was requested
if hasattr(args,'outputfile'):
    outputfile.flush()
    outputfile.close()
##
### end of program 