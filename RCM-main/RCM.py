# Programed in Python3 but run in python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 15:39:49 2020

@author: carlos
"""


import sys
import subprocess as sp
import os

if __name__=="__main__":

                
    print("\nRCM - Reverse Cuthill Mckee Algorithm\n")    


    os.system("gcc ./readTables.c ./utils.c ./mesh.c ./input.c ./RCM.c ./main.c -o ./RCM -lm -fopenmp -Wall -O3")
    os.system("./RCM")

