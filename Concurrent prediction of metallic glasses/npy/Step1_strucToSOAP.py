# This file converts the lammps structure file into ase.Atoms type and then get
#the SOAP descriptor of the entile structure
from ase import Atoms   #Import ase.Atoms type
from dscribe.descriptors import SOAP  #Import SOAP descriptor
import math
import numpy as np
import shelve
import os
import datetime

def All_soap(rcut,nmax,lmax,sigma):        #Get the SOAP descriptors of all structures
    ##################### Loop to read all the structure file #####################
    Global_variable = globals()

    for coolingrate in ['13_1']:
        rootpath = '../Data/Cooling'
        for temperature in [500, 550, 600, 650, 700]:
            stru_num = 0         #Store the number of structures
            for time in range(1, 201):
                #################### Read the structure file  ##################################
                filename = 'ML' + coolingrate + '_' + str(temperature) + '_' + str(time) + '.txt'
                fullpath = rootpath + coolingrate + '/' + str(temperature) + '/Sort/' + filename 
                try:
                    lmp_structure = open(fullpath,'r')    #Name of structure file
                except:
                    print ("Cannot open file %s" % filename)
                    return 0
                
                #Notice the structure file should be ordered, the atom ID should increase !
                #And type 1 (Zr) comse first and type 2 (Cu) comes after
                print ("File %s is open" % filename)
                stru_num += 1
                for i in range(1):
                    lmp_structure.readline()
                atom_num = int(lmp_structure.readline().split()[0])  #Number of atoms in the structure file
                for i in range(2):
                    lmp_structure.readline()
                xbound = lmp_structure.readline().split()
                xlow = float(xbound[0])
                xhigh = float(xbound[1])
                ybound = lmp_structure.readline().split()
                ylow = float(ybound[0])
                yhigh = float(ybound[1])
                zbound = lmp_structure.readline().split()
                zlow = float(zbound[0])
                zhigh = float(zbound[1])
                for i in range(8):
                    lmp_structure.readline()
                    
                structure = []
                for i in range(atom_num):
                    line = lmp_structure.readline().split()
                    structure.append((float(line[2]),float(line[3]),float(line[4])))
                lmp_structure.close()

                ################################### Transfer to ase.Atoms type ##################
                Zr48Cu52 = Atoms('Zr960Cu1040',positions = structure)
                Zr48Cu52.set_pbc((True,True,True))   #Set periodic boundary condition
                xlength = xhigh - xlow
                ylength = yhigh - ylow
                zlength = zhigh - zlow
                Zr48Cu52.set_cell([xlength,ylength,zlength]) #Set unit cell  #Has been Checked

                ############################ Get SOAP descriptor ###########################
                species = ["Zr","Cu"]

                soap = SOAP(         #Set up SOAP constructor
                species=species,
                periodic=True,
                rcut=rcut,
                nmax=nmax,
                lmax=lmax,
                sigma=sigma,
                rbf='polynomial'
                )

                Global_variable['soap_MG' + str(stru_num)] = soap.create(Zr48Cu52)
                np.save(coolingrate + '_' + str(temperature) + '_soap_MG' + str(stru_num),Global_variable['soap_MG' + str(stru_num)])
                print(stru_num)
                
                del Global_variable['soap_MG' + str(stru_num)]

    return 1


################################## Run main part ##############################
rcut = 5.0
nmax = 11
lmax = 12
sigma = 0.5

All_soap(rcut,nmax,lmax,sigma)

currentDT = datetime.datetime.now()
print(str(currentDT))
