#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
----------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          T. Cisneros-Perez
---------------------------------------------------------------
"""
import os
import urllib.request

LTVERS="LoopTools-2.16"

def main(arg1):
    """
    Function to build a shared library for LoopTools:
    - Download the fortran code
    - Compile the code with appropiate flags
    - Create the shared library
    """  
    print("Fetching LoopTools...")
    os.chdir(arg1)
    url = "http://www.feynarts.de/looptools/" + LTVERS + ".tar.gz"
    urllib.request.urlretrieve(url, "looptools.tar.gz")
    os.system('tar xzvf looptools.tar.gz')
    os.system("cp ./config/configure ./" + LTVERS + "/")
    os.chdir(arg1 + LTVERS)
    os.system('./configure CFLAGS="-fomit-frame-pointer -fPIC" FFLAGS="-fPIC -ff2c" CXXFLAGS="-fomit-frame-pointer -fPIC"')
    os.system("make")
    os.chdir("..")
    print("Building shared library...")
    os.system("gfortran -shared -o liblooptools.so " + os.getcwd() + "/" + LTVERS + "/build/*.o")
    print("Cleaning up directory...")
    os.system("rm -rf looptools.tar.gz")
    os.system("rm -rf " + LTVERS + "/")
    
    
if __name__ == "__main__":
    main("./")
