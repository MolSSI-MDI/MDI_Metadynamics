## Description

This repository provides a reference implementation of the Metadynamics method 
using MolSSI MDI. The free energy of dissociation of NaCl in water is used as
a test calculation and the input files are provided.

## Overview of the implementation of Metadynamics using MDI

In this Metadynamics implementation, a single MD engine instance 
establishes a connection with
the Metadynamics driver, which takes control of the calculation. After
 the driver recieves the
nuclear coordinates, the cell vectors and the forces from the engine, a history
dependent bias is computed as a function of the distance of the NaCl ions, adds 
this bias to the nuclear forces, and sends the updated nuclear forces to the engine.

## Installation

This driver requires 

1. LAMMPS, installed through the MDI fork 
2. The MDI Library
3. This driver

Edit the LAMMPS and MDI_Metadynamics files located within the
MDI_Metadynamics/tests/locations folder. Edit the file paths to point to the
LAMMPS executable and the driver binary, respectively.

## Run NaCl-water dissociation calculation

The directory test/npt_spcfw_nacl/ contains the relevant files to execute
a NaCl-water dissociation using MDI. The directory test/npt_spcfw_nacl/data
contains the LAMMPS and MDI_Metadynamics input files as well as reference
data from the SSAGES code for a similar system.

The script test/npt_spcfw_nacl/tcp.sh is used to start the execution
of the simulation. The script test/npt_spcfw_nacl/extract_G.py is used
to analyze the output of MDI_Metadynamics and produce a final plot for
the Gibbs free energy as a function of ionic distance.
