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

echo
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
rm -rf OUTPUT_FILES/BAK_FROM_MESHER
mkdir OUTPUT_FILES/BAK_FROM_MESHER
rm -f OUTPUT_FILES/DATABASES_MPI/*


# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/Par_file_ASKI OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/FORCESOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp ../../setup/constants.h OUTPUT_FILES/
cp process.sh OUTPUT_FILES/
cp process_solver_only.sh OUTPUT_FILES/
cp run_specfem3dCartesianForASKI_simulations.py OUTPUT_FILES/

sleep 10

# decomposes mesh
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC ./MESH ./OUTPUT_FILES/DATABASES_MPI/

sleep 5

# runs database generation
echo
echo "  running database generation..."
echo
mpirun -np $NPROC ./bin/xgenerate_databases

cp OUTPUT_FILES/constants.h OUTPUT_FILES/BAK_FROM_MESHER/
cp OUTPUT_FILES/output_mesher.txt OUTPUT_FILES/BAK_FROM_MESHER/
cp OUTPUT_FILES/surface_from_mesher.h OUTPUT_FILES/BAK_FROM_MESHER/
cp OUTPUT_FILES/values_from_mesher.h OUTPUT_FILES/BAK_FROM_MESHER/
cp OUTPUT_FILES/LOG_ASKI_model_external.txt OUTPUT_FILES/BAK_FROM_MESHER/

sleep 10

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
