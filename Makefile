#----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
################################################################
#  This is the Makefile for the extension package SPECFEM3D_Cartesian for ASKI (for GNU Make)
################################################################
#
#-----------------------------------------------------------------------
#  set the compiler
#
COMPILER = gfortran
#
#-----------------------------------------------------------------------
#  General definitions
#
bindir = ../ASKI/bin
obsdir = ../ASKI/obj
#
FFLAGS = -O3 -J$(obsdir) -I/usr/include -Wunused-variable -Wuninitialized -fimplicit-none -ffixed-line-length-132 -fbounds-check -fbacktrace
#
#-----------------------------------------------------------------------
#  Direcories where to search for files to compile to .o by implicit rules below, and dependencies defined in rules.mk
#
vpath %.o $(obsdir)
vpath %.f90 ../ASKI/f90
#
#-----------------------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.f90
	$(COMPILER) -c $(FFLAGS) $< -o $(obsdir)/$@
#
#-----------------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
#
#-----------------------------------------------------------------------
#  Library paths
#
BLAS = /usr/lib/libblas.so
LAPACK = /usr/lib/liblapack.so
#
#-------------------------------------------------------------
#
.PHONY:
#
#----------------------------------------------------------------
#  Include dependencies:
#  rules_SPECFEM3D_Cartesian.mk and ../ASKI/rules.mk are Makefiles because they are included. They
#  contain all dependencies of the .o files. If you change any such dependencies (e.g. by using an
#  additional module in some program/module), please update files rules_SPECFEM3D_Cartesian.mk , 
#  ../ASKI/rules.mk accordingly.
#
-include rules_SPECFEM3D_Cartesian.mk
-include ../ASKI/rules.mk
#
#---------------------------------------------------------------
#
clean:
	-rm -f $(bindir)/transformSpecfem3dCartesianSyntheticData
	-rm -f $(bindir)/transformSpecfem3dCartesianMeasuredData
	-rm -f $(bindir)/convolveWithStf
	-rm -f $(obsdir)/*
#
#----------------------------------------------------------------
# Rules for all programs of SPECFEM3D for ASKI:
#
transformSpecfem3dCartesianSyntheticData: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o seismicStation.o asciiDataIO.o \
	seismicEventList.o inversionBasics.o discreteFourierTransform.o componentTransformation.o mathConstants.o \
	seismicEvent.o argumentParser.o string.o realloc.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o flexibleType.o dateTime.o \
	parameterCorrelation.o readEventStationFile.o modelParametrization.o geminiKernelReferenceModel.o \
	specfem3dKernelReferenceModel.o nexdKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o nexdWavefieldPoints.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o vectorPointer.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o timeUtils.o locatePoint.o streamAccess.o chunkCubedSphere.o externalRadialNodes.o \
	scart2dGrid.o specfem3dForASKI_mod.o dataModelSpaceInfo.o fourierTransform.o specfem3dForASKIFiles.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

transformSpecfem3dCartesianMeasuredData: %: %.o specfem3dForASKI_mod.o errorMessage.o seismicNetwork.o fileUnitHandler.o seismicStation.o asciiDataIO.o \
	seismicEventList.o inversionBasics.o discreteFourierTransform.o componentTransformation.o complexKernelFrequency.o mathConstants.o \
	seismicEvent.o argumentParser.o string.o fourierTransform.o realloc.o flexibleType.o dateTime.o parameterCorrelation.o \
	readEventStationFile.o inputParameter.o modelParametrization.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

convolveWithStf: %: %.o specfem3dForASKI_mod.o inversionBasics.o seismicEvent.o seismicEventList.o seismicStation.o seismicNetwork.o \
	componentTransformation.o asciiDataIO.o fileUnitHandler.o argumentParser.o string.o errorMessage.o inputParameter.o fourierTransform.o \
	parameterCorrelation.o readEventStationFile.o modelParametrization.o mathConstants.o flexibleType.o dateTime.o realloc.o \
	primitiveTypeEncoding.o simpleString.o kindDefinitions.o timeUtils.o
	$(COMPILER) -o $(bindir)/$@ $(obstring)

all: transformSpecfem3dCartesianSyntheticData transformSpecfem3dCartesianMeasuredData convolveWithStf
