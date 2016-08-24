#!/bin/bash
#
# script runs decomposition,database generation and solver
# using this example setup
#
# prior to running this script, you must create the mesh files
# in directory MESH/
# (see section 3.1 "Meshing with CUBIT" in user guide)
#

###################################################

# number of processes
NPROC=16

##################################################

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 5 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo


rm -f OUTPUT_FILES/*


# (re)stores setup
cd OUTPUT_FILES/BAK_FROM_MESHER
cp values_from_mesher.h output_mesher.txt constants.h surface_from_mesher.h ../
cd ../../
cp DATA/Par_file OUTPUT_FILES/
cp DATA/Par_file_ASKI OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/FORCESOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp process.sh OUTPUT_FILES/
cp process_solver_only.sh OUTPUT_FILES/
cp run_specfem3dCartesianForASKI_simulations.py OUTPUT_FILES/

sleep 20

# runs simulation
echo
echo "  running solver..."
echo
mpirun -np $NPROC ./bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`

sleep 20

