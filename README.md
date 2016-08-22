# Introduction

This repository contains code created by Joelle Barral to process gold standard inversion recovery T1 maps. The original
code was downloaded from [http://www-mrsrl.stanford.edu/~jbarral/t1map.html](http://www-mrsrl.stanford.edu/~jbarral/t1map.html)
and contains some additional Matlab files written by Yuhan Ma.

I created this repository on Github so that it could be used as a submodule in other repositories. Joelle's original code will
not be modified, any modifications will be created by copying the files in separate directories and renaming them.

# Joelle's original README

T<sub>1</sub> Mapping Matlab package
------------------------

This package includes Matlab functions to perform T<sub>1</sub> mapping. 

It accompanies the journal article:

Joelle K. Barral, Erik Gudmundson, Nikola Stikov, Maryam Etezadi-Amoli,Petre Stoica, and Dwight G. Nishimura "*A robust methodology for in vivo T<sub>1</sub> mapping*". Magn. Reson. Med., 64: 1057–1067. (2010) doi:10.1002/mrm.22497

----------------------------------------------------------------------------

All functions mentioned below can be found in mfiles and are documented. 
The user will be asked for input in the command window

A dataset (DICOM images in singleslicedicomfiles) is given as example. GE's high resolution phantom was scanned at 1.5 T with the SE-IR sequence and the following parameters:
TR = 2550 ms, TE = 14 ms, TI = [50, 400, 1100, 2500] ms, BW = 32 kHz, FOV = 20 cm, Slice thickness = 5 mm, Matrix size 512x128. 


>>> To run simulations: 

Run mainSim (change parameters as desired, default takes a while!, decrease MC
for a quick test)
or
Run compareMethodsSim (add methods as desired)

>>> To fit experimental data: 

In getData: 
   - Change T1path appropriately
   - Specify loadpath and savename 
     (default will use the dicom images given as example)
		
In main.m
   - Change T1path appropriately
   - Specify loadpath, filename, and savename 
     (default will use the example)

Run getData
then
Run mainScan

----------------------------------------------------------------------------

Please send any feedback to Joelle Barral at jbarral@mrsrl.stanford.edu.
This package is available online at
http://www-mrsrl.stanford.edu/~jbarral/t1map.html

Last updated June 02, 2010
Written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov
(c) 2009, 2010, Board of Trustees, Leland Stanford Junior University



