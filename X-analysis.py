#!/usr/bin/env python
# This program approximates the average electron binding energy (X) 
# from the eigenvalues of molecular orbitals 
# An "Experimental Quantum Chemistry" energy decomposition analysis is also performed. 
# Both closed and open-shell wavefunctions are OK
# Requirements: numpy, cclib, >Python 3.0. 

#USAGE:
# Analysis on single file: 
# ./X-analysis.py QM-outputfile (tested on Gaussian or ORCA, but should work on all programs readable by cclib)

#Analysis of a step A-->B, where A and B are two different structures or electronic states: 
# ./X-analysis.py QM-outputfile_A QM-outputfile_B 

#Analysis of a bonding process A + B --> C:  
# ./X-analysis.py QM-outputfile_A QM-outputfile_B QM-outputfile_C 

#OPTIONS:
# If core-valence separation is desired add the flag --core
# To specify stochiometry X and Y in XA --> 1B or XA + YB --> 1AB add the flag --stoch

# Written by Martin Rahm, Chalmers University of Technology 
#MIT License
#Copyright (c) [2018] [Martin Rahm]

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

print('')
print('>  EQC ENERGY DECOMPOSITION ANALYSIS ')
print('>  Energies in electron volt (eV)   ')
print('>  Revision Jul 20, 2018            ')
print('')
print('>  Citation is greatfully appreciated:           ')
print('>  https://rahmlab.com/energy-decomposition-analysis/ ')
print('>  M. Rahm, R. Hoffmann, J. Am. Chem. Soc, 137, 10282-10291, 2015')
print('>  M. Rahm, R. Hoffmann, J. Am. Chem. Soc, 138, 3731-3744, 2016')
print('')

import cclib
from cclib.io import ccopen 
import sys
import os 
import numpy

inputfile = sys.argv

def Xanalysis(inputfile,core): #Main program code, does analysis on one QM output file
    fo = open("{0}".format(inputfile), 'r') #read in input
    lines = fo.readlines()
    fo.close()

    if core % 2 != 0:
        print("WARNING! Only even number of core electrons are supported. Core averaging will be incorrect!")
    core=int(core/2) # translates number of electrons to number of orbitals

    replace=False   
    for line in lines:
        if "**********-" in line:
            print("Deep core levels are not printed correctly in the Gaussian output. These will be ignored")
            replace=True
            break

    if replace == True: 
        fi = open("TEMP-{0}".format(inputfile), 'w') #create new temporary output file

        for line in lines:
            if "**********-" in line:
                a = line.replace("**********-", "0.0000000 -")
                fi.write(a)
            else:
                fi.write(line)
            if "Error termination" in line: #avoids output that ends with an error
                print("Error termination in Gaussian output detected")
                break
        fi.close()
           
    if replace == False:
        data = ccopen(inputfile)
    else:
        data = ccopen("TEMP-{0}".format(inputfile))

    data = data.parse()
    electrons = 0

    if replace == True:
        os.remove("TEMP-{0}".format(inputfile)) # removes temporary file needed in case deep core levels could not be printed in Gaussian output

    try:
        orbitals = data.homos[0] +1 + data.homos[1]+1 
        orbitals_alpha = data.homos[0]+1
        orbitals_beta = data.homos[1]+1
        print('This wavefunction is open-shell')
        open_shell=True
        electrons=orbitals
    except IndexError:
        orbitals = data.homos[0] +1
        open_shell=False
        print('This wavefunction is closed-shell')
        electrons=2*orbitals 

    print("System contains %d electrons" % electrons)
    print("There are %d occupied orbitals" % orbitals)

    if open_shell==True:
        energies_alpha=data.moenergies[0]
        energies_beta=data.moenergies[1]
    else:
        energies_alpha=data.moenergies[0]
        orbitals_alpha = orbitals 
    sum_val=0
    sum_core=0
    sum_tot=0

    if open_shell==True:
        valence_alpha = energies_alpha[core:orbitals_alpha] 
        valence_beta = energies_beta[core:orbitals_beta] 
        core_alpha = energies_alpha[:core] 
        core_beta = energies_beta[:core] 
        total_alpha = energies_alpha[:orbitals_alpha]
        total_beta = energies_beta[:orbitals_beta] 

    else:
        valence_alpha = energies_alpha[core:orbitals] 
        core_alpha = energies_alpha[:core] 
        total_alpha = energies_alpha[:orbitals]
        valence_beta = [0]
        total_beta = [0]
        core_beta = [0]

    # List Orbital Energies 
    m=0
    LUMO=0.0
    HOMO=0.0
    print("")
    print("----ORBITALS--------------------------")
    if open_shell==False:
        print("Energies of Occupied Orbitals:")
        for e in data.moenergies[0]:
            if m<orbitals:
                m+=1
                print(e)
        m=-1
        for e in data.moenergies[0]: #Search for LUMO
            if m<orbitals:
                m+=1
            if m==orbitals:
                m+=1
                LUMO=e #save LUMO
                print("LUMO:", e)        
        HOMO=valence_alpha[len(valence_alpha)-1]

    if open_shell==True:
        print("Energies of Occupied Alpha Orbitals:")
        for e in data.moenergies[0]:
            if m<orbitals_alpha:
                m+=1
                print(e)
        m=0
        print("Energies of Occupied Beta Orbitals:")
        for e in data.moenergies[1]:
            if m<orbitals_beta:
                m+=1
                print(e)
        m=-1
        for e in data.moenergies[1]: #Search for beta-LUMO 
            if m<orbitals_beta:
                m+=1
            if m==orbitals_beta:
                m+=1
                LUMO=e #save beta-LUMO  
                print("LUMO:", e) 
        try:
            if valence_alpha[len(valence_alpha)-1]>valence_beta[len(valence_beta)-1]:
                HOMO=valence_alpha[len(valence_alpha)-1]
            else:
                HOMO=valence_beta[len(valence_beta)-1]
        except IndexError:
            if len(valence_beta)==0:
                HOMO=valence_alpha[len(valence_alpha)-1]
            if len(valence_alpha)==0:
                HOMO=valence_beta[len(valence_alpha)-1]

    for MO in range(len(valence_alpha)): #Sums the orbital energies
        sum_val += valence_alpha[MO]
        MO += 1
    for MO in range(len(valence_beta)):
        sum_val += valence_beta[MO]
        MO += 1
    for MO in range(len(total_alpha)):
        sum_tot += total_alpha[MO]
        MO += 1
    for MO in range(len(total_beta)):
        sum_tot += total_beta[MO]
        MO += 1
    for MO in range(len(core_alpha)):
        sum_core += core_alpha[MO]
        MO += 1
    for MO in range(len(core_beta)):
        sum_core += core_beta[MO]
        MO += 1

    if open_shell==True:
        valence = len(valence_alpha) + len(valence_beta)
        total = len(total_alpha) + len(total_beta)
        core = len(core_alpha) + len(core_beta)  #NOTE: Open-shell cores are not yet treated correct! Equal number of core alpha and beta is assumed
        average_val = sum_val / valence 
        if core == 0: #If no core exists (e.g. for hydrogen), the valence is also shown as core
            average_core = average_val
        else:
            average_core = sum_core / core

        average_tot = sum_tot / total

    else: #Both unrestricted and restricted count per electron, so *2 here. 
        valence = 2*len(valence_alpha)
        total = 2*len(total_alpha)
        core = 2*len(core_alpha)
        average_val = (2*sum_val) / valence 
        if core == 0: #If no core exists (e.g. for hydrogen), the valence is also shown as core
            average_core = average_val
        else:
            average_core = (sum_core*2) / core
        average_tot = (sum_tot*2) / total

    print("HOMO:", HOMO)
    print("GAP:", HOMO - LUMO)
    print("")    
    print("----TOTAL ENERGIES-------------------")    
    Etot = data.scfenergies
    print("E = n*X-bar + Vnn + Eee")    
    print("Total SCF energy (E):                                   ", Etot[-1]) #print last energy
    if core != 0:
        print("Total energy of valence electrons (n*X-bar_valence):    ", (electrons-core)*average_val)
        print("Total energy of core orbitals (n*X-bar_core):           ", core*average_core)
    print("Total energy of occupied orbitals (n*X-bar):            ", electrons*average_tot) 

    steps = data.atomcoords.shape[0]-1 #Determines number of optimization steps

    def nuclear_repulsion_energy(data): #Defines function for calculating N-N repulsion (written by Karol Langner)
        nre = 0.0
        for i in range(data.natom):
            ri = data.atomcoords[steps][i]
            zi = data.atomnos[i]
            for j in range(i+1, data.natom):
                rj = data.atomcoords[steps][j]
                zj = data.atomnos[j]
                d = numpy.linalg.norm(ri-rj)
                nre += zi*zj/d
        return nre
            
    a2b = cclib.parser.utils.convertor(1.0, "Angstrom", "bohr")
    VNN = 27.211396132*(nuclear_repulsion_energy(data)) / a2b
    print("Nuclear repulsion (Vnn):                                ", VNN)

    W=0 # Calculate multi-electron term (Eee)

    if open_shell==True:
        W = Etot[-1] - total*average_tot-VNN
        print("Multi-electron term (Eee):                              ", -W) # Eee = -W

    if open_shell==False:
        W = Etot[-1] - total*average_tot-VNN
        print("Multi-electron term (Eee):                              ", -W)
    print("Vnn - Eee:                                              ", (VNN+W))

    print("")
    print("----TOTAL ENERGIES PER ELECTRON-------")
    print("E/n = X-bar + (Vnn - Eee)/n")    
    print("E/n:             ", Etot[-1]/electrons) #print last energy  
    if core != 0:
        print("X-bar_valence:   ", average_val)
        print("X-bar_core:      ", average_core)
    print("X-bar:           ", average_tot)
    print("Vnn/n:           ", VNN/electrons)
    print("Eee/n:           ", -W/electrons)
    print("(Vnn - Eee)/n:   ", (VNN/electrons+W/electrons))

    print("")
    print("-------MULTIELECTRON ANALYSIS----------")
    print("100*(|Eee/Etot|):  ",format((W/Etot[-1])*100,'.2f') + " %") #prints the Eee-percentage with two decimal points  

    E = Etot[-1]
    En = Etot[-1]/electrons
    Vnnn = VNN/electrons
    Wn = W/electrons
    WnVnnn = Vnnn + Wn

    return(average_core, average_val, average_tot, E, En, Vnnn, Wn, WnVnnn, electrons)

class NullWriter(object):
    def write(self, arg):
        pass
oldstdout = sys.stdout
nullwrite = NullWriter()
#sys.stdout = nullwrite # disable output 
#sys.stdout = oldstdout # enable output 

# Determine type of calculation based on the number of files to be analyzed: 
valencesplit = False 

# Analysis of single file
if (len(inputfile) == 2) or (len(inputfile) == 3) and (("--core" in inputfile)):
    STATECALC = True 
    print("* "+"{0}".format(inputfile[1]) + " *")
    file1=inputfile[1]
    if "--core" in inputfile:
        core = eval(input ("Please enter number of core electrons in " + file1 + " (0 = all electrons included): "))
        valencesplit = True
    else:
        core = 0
    results=Xanalysis(file1,core)

#Analysis of a step A-->B, where A and B are isomers 
if ((len(inputfile) == 3 ) and ("--core" not in inputfile)) or ((len(inputfile) == 4) and ("--core" in inputfile)) or ((len(inputfile) == 4) and ("--stoch" in inputfile)) or ((len(inputfile) == 5) and ("--core" in inputfile) and ("--stoch" in inputfile)):
    STEPCALC = True
    stoch = 1.0 #stoch is stochiometry of A, e.g. 2A --> B. Default is 1. 
    print("* "+"{0}".format(inputfile[1]) + " --> " + "{0}".format(inputfile[2]) + " *")
    file1=inputfile[1]
    file2=inputfile[2]
    if "--core" in inputfile:
        core = eval(input ("Please enter number of core electrons in " + file1 + " (0 = all electrons included): "))
        valencesplit = True 
        if core % 2 != 0:
            print("WARNING! Only even number of core electrons are supported. Core averaging will be incorrect!")
    else:
        core = 0
    if "--stoch" in inputfile:
        stoch = eval(input ("Please enter the stochiometry of " + file1 + " (1 is default): "))
        
    sys.stdout = nullwrite # disable output 
    results1=Xanalysis(file1,core)
    results2=Xanalysis(file2,core*stoch)
    sys.stdout = oldstdout # enable output 

    print(results1[0], results2[1])

    dX = results2[2] - results1[2]
    dE = results2[3] - stoch*results1[3]
    dVnn = results2[5] - results1[5]
    Q = ((2*results2[8]*dX)/dE)-1
    dEee = results2[6] - results1[6]
    dEeeVnnn = results2[7] - results1[7]
    Covalency = (100*(abs(dX)-dX)/2)*(1/(abs(dEeeVnnn)-dX)) + (100*(abs(dX)+dX)/2)*(1/(dEeeVnnn-dX)+1/dX)  #from SI of JACS,137,10282-10291,2015. 

    print("E/n = X-bar + (Vnn - Eee)/n")    
    print("-------------------------------------------------")
    print("Electrons (n):           " , results2[8] )
    print("Delta(E):                " , format(dE,'.3f') + " eV")
    print("Delta(E/n):              " , format(dE/results2[8],'.3f') + " eV/e")
    print("Delta(X-bar):            " , format(dX,'.3f') + " eV/e")
    print("Delta(Vnn/n):            " , format(dVnn,'.3f') + " eV/e")
    print("Delta(Eee/n):            " , format(-dEee,'.3f') + " eV/e")
    print("Delta(Vnn-Eee)/n:        " , format(dEeeVnnn,'.3f') + " eV/e")
    print("Q:                       " , format(Q,'.3f') )
    print("Delta(|Eee/E|)           " , format(((dEee/dE)*100),'.2f') + " %")
    if valencesplit == True:
        print(results2[0])
        print(results1[0])
        dXcore = results2[0] - results1[0]
        dXval = results2[1] - results1[1]
        dXapprox = (dXval*(results2[8]-core*stoch))/results2[8]
        print("Delta(X-bar)_valence:    " , format(dXval,'.3f') + " eV/e")
        print("Delta(X-bar)_core:       " , format(dXcore,'.3f') + " eV/e")
        print("Delta(X-bar)_val-approx: " , format(dXapprox,'.3f') + " eV/e")
    print("Covalency Index:         " , format(Covalency,'.1f') + " %")
    print("Ionicity Index:          " , format(100-Covalency,'.1f') + " %")

#Analysis of bond formation, A + B --> AB
if ((len(inputfile) == 4 ) and ("--core" not in inputfile) and ("--stoch" not in inputfile)) or ((len(inputfile) == 5) and ("--core" in inputfile) and ("--stoch" not in inputfile)) or ((len(inputfile) == 5) and ("--core" not in inputfile) and ("--stoch" in inputfile)) or ((len(inputfile) == 6) and ("--core" in inputfile) and ("--stoch" in inputfile)):
    BONDCALC = True 
    print("* "+"{0}".format(inputfile[1]) + " + " + "{0}".format(inputfile[2]) + " --> " +  "{0}".format(inputfile[3]) + " *")
    file1=inputfile[1]
    file2=inputfile[2]
    file3=inputfile[3]
    stoch1=1.0
    stoch2=1.0
    if "--core" in inputfile:
        core1 = eval(input ("Please enter number of core electrons in " + file1 + " (0 = all electrons included): "))
        core2 = eval(input ("Please enter number of core electrons in " + file2 + " (0 = all electrons included): "))
        valencesplit = True 
    else:
        core1 = 0
        core2 = 0
    if "--stoch" in inputfile:
        stoch1 = eval(input ("Please enter the stochiometry of " + file1 + " (1 is default): "))
    if "--stoch" in inputfile:
        stoch2 = eval(input ("Please enter the stochiometry of " + file2 + " (1 is default): "))

    core3 = core1*stoch1 + core2*stoch2 
    if core3 % 2 != 0:
        print("WARNING! Only even number of core electrons are supported. Core averaging will be incorrect!")

    sys.stdout = nullwrite # disable output  
    results1=Xanalysis(file1,core1)
    results2=Xanalysis(file2,core2)
    results3=Xanalysis(file3,core3)
    sys.stdout = oldstdout # enable output 
    
    #results(average_core, average_val, average_tot, E, En, Vnnn, Wn, WnVnnn, electrons)    
 
    dX = results3[2] - ((stoch1*results1[8]*results1[2] + stoch2*results2[8]*results2[2])/results3[8])
    dE = results3[3] - (stoch1*results1[3]+stoch2*results2[3])
    dVnn = results3[5] - ((results1[8]*results1[5] + results2[8]*results2[5])/results3[8])
    Q = ((2*results3[8]*dX)/dE)-1
    dEee = results3[6] - ((stoch1*results1[8]*results1[6] + stoch2*results2[8]*results2[6])/results3[8]) 
    dEeeVnnn = results3[7] - ((stoch1*results1[8]*results1[7] + stoch2*results2[8]*results2[7])/results3[8])
    Covalency = (100*(abs(dX)-dX)/2)*(1/(abs(dEeeVnnn)-dX)) + (100*(abs(dX)+dX)/2)*(1/(dEeeVnnn-dX)+1/dX)  #from SI of JACS,137,10282-10291,2015. 

#    print (abs*)
 
    print
    print("E/n = X-bar + (Vnn - Eee)/n")    
    print("-------------------------------------------------")
    print("Electrons (n):           " , results3[8] )
    print("Delta(E):                " , format(dE,'.3f') + " eV")
    print("Delta(E/n):              " , format(dE/results3[8],'.3f') + " eV/e")
    print("Delta(X-bar):            " , format(dX,'.3f') + " eV/e")
    print("Delta(Vnn/n):            " , format(dVnn,'.3f') + " eV/e")
    print("Delta(Eee/n):            " , format(-dEee,'.3f') + " eV/e")
    print("Delta(Vnn-Eee)/n:        " , format(dEeeVnnn,'.3f') + " eV/e")
    print("Q:                       " , format(Q,'.3f') )
    print("Delta(|Eee/E|)           " , format(((dEee/dE)*100),'.2f') + " %")
    if valencesplit == True:
        dXcore = results3[0] - ((stoch1*core1*results1[0] + stoch2*core2*results2[0])/core3)
        dXval = results3[1] - (((results1[8]-core1)*stoch1*results1[1] + (results2[8]-core2)*stoch2*results2[1])/(results3[8] - core3))
        dXapprox = (dXval*(results3[8]-core3))/results3[8]
        print("Delta(X-bar)_valence:    " , format(dXval,'.3f') + " eV/e")
        print("Delta(X-bar)_core:       " , format(dXcore,'.3f') + " eV/e")
        print("Delta(X-bar)_val-approx: " , format(dXapprox,'.3f') + " eV/e")
    print("Covalency Index:         " , format(Covalency,'.1f') + " %")
    print("Ionicity Index:          " , format(100-Covalency,'.1f') + " %")



