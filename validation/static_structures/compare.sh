#!/usr/bin/env bash

# 1 eV in Kcal/mol (NIST)
ECONV=23.060549

# Executables
LMPEXE=${HOME}/bin/lmp_mpi
MCEXE=${HOME}/bin/mw_water_ls

# Directories containing structures to test
for dir in `ls -d */`
do

    # Change and report directory
    cd $dir
    echo $dir
    echo
    
    # Clear stale output
    rm -f  *log checkpoint*.dat.? ice???_therm.dat *.psf *.dcd log.lammps
    
    # LAMMPS energy of structure
    LMPEN=`${LMPEXE} -in ../in.solid | grep -A1 Step | grep -v Step | awk '{print $2}'`
    
    # mw_water_ls energy of structure converted to LAMMPS kcal/mol
    ${MCEXE} ice.input 
    MCEN=`grep 'Computed energy' node000.log  | awk '{print $5}'`
    
    # Convert LAMMPS energy to eV
    LMPEN=`echo $LMPEN | awk -v c=${ECONV} '{print $1/c}'`
    
    printf "  LAMMPS (eV)               MW_WATER_LS (eV)    \n" 
    printf "  %10.6f                 %10.6f                 \n" $LMPEN $MCEN
    echo
    
    cd ../

done
