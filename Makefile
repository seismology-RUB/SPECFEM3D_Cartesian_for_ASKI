###############################################################
#  This is the Makefile for SPECFEM3D_Cartesian 3.0 for ASKI 1.0
###############################################################
#
#-----------------------------------------------------------------------
#  set the compiler
#
COMPILER = gfortran
#
#-----------------------------------------------------------------------
#  General definitions
#
bindir = ../bin
obsdir = ../obj
#
FFLAGS = -O3 -J$(obsdir) -I/usr/include -Wunused-variable -Wuninitialized -fimplicit-none -ffixed-line-length-132 -fbounds-check -fbacktrace
#
#-----------------------------------------------------------------------
#  Direcories where to search for files to compile to .o by implicit rules below, and dependencies defined in rules.mk
#
vpath %.o $(obsdir)
vpath %.f90 ../f90
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
BLAS = /usr/lib/libblas.so.3gf
LAPACK = /usr/lib/liblapack.so.3gf
#
#-------------------------------------------------------------
#
.PHONY:
#
#----------------------------------------------------------------
#  Include dependencies:
#  rules.mk and ../rules.mk are a Makefile because it is included. They contain all dependencies of
#  the .o files. If you change any such dependencies (e.g. by using an additional module
#  in some program/module), please update files rules.mk , ../rules.mk accordingly.
#
-include rules_SPECFEM3D_Cartesian.mk
-include ../rules.mk
#
#---------------------------------------------------------------
#
clean:
	-rm -f $(bindir)/transformSpecfem3dCartesianSyntheticData
	-rm -f $(bindir)/transformSpecfem3dCartesianMeasuredData
	-rm -f $(bindir)/createSpecfem3dFilters
	-rm -f $(obsdir)/*
#
#----------------------------------------------------------------
# Rules for all programs of SPECFEM3D for ASKI:
#
transformSpecfem3dCartesianSyntheticData: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o seismicStation.o asciiDataIO.o \
	seismicEventList.o inversionBasics.o discreteFourierTransform.o componentTransformation.o mathConstants.o \
	seismicEvent.o argumentParser.o string.o realloc.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o flexibleType.o dateTime.o \
	parameterCorrelation.o readEventStationFile.o modelParametrization.o geminiEarthModel.o \
	specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o vectorPointer.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o timeUtils.o locatePoint.o streamAccess.o chunkCubedSphere.o externalRadialNodes.o \
	scart2dGrid.o specfem3dForASKI_mod.o dataModelSpaceInfo.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

transformSpecfem3dCartesianMeasuredData: %: %.o errorMessage.o seismicNetwork.o fileUnitHandler.o seismicStation.o asciiDataIO.o \
	seismicEventList.o inversionBasics.o discreteFourierTransform.o componentTransformation.o mathConstants.o \
	seismicEvent.o argumentParser.o string.o realloc.o invgridVtkFile.o eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o \
	inversionGrid.o kernelInvertedModel.o wpVtkFile.o inputParameter.o integrationWeights.o flexibleType.o dateTime.o \
	parameterCorrelation.o readEventStationFile.o modelParametrization.o geminiEarthModel.o \
	specfem3dKernelReferenceModel.o geminiWavefieldPoints.o specfem3dWavefieldPoints.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o schunkInversionGrid.o chunksInversionGrid.o vectorPointer.o primitiveTypeEncoding.o \
	simpleString.o kindDefinitions.o timeUtils.o locatePoint.o streamAccess.o chunkCubedSphere.o externalRadialNodes.o \
	scart2dGrid.o specfem3dForASKI_mod.o dataModelSpaceInfo.o specfem3dForASKIFiles.o nexdWavefieldPoints.o nexdKernelReferenceModel.o
	$(COMPILER) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

createSpecfem3dFilters: %: %.o discreteFourierTransform.o asciiDataIO.o errorMessage.o mathConstants.o realloc.o
	$(COMPILER) -o $(bindir)/$@ $(obstring)

all: transformSpecfem3dCartesianSyntheticData transformSpecfem3dCartesianMeasuredData
