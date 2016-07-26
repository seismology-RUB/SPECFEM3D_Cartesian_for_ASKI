#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  This is the generic part of the Makefile. Use this
#  template within each project.
#-----------------------------------------------------------------------
#  set the SHELL
#
SHELL = /usr/bin/tcsh
#----------------------------------------------------------------
#  General definitions
#
bindir = ./bin
obsdir = ./obj
moduledir = ./mod
#
ifeq ($(notdir $(F95)),g95)
	FFLAGS = -O3 -Wunused -fmod=$(moduledir)
else
	FFLAGS = -O3 -J$(moduledir) -Wunused-variable -Wuninitialized -fimplicit-none -ffixed-line-length-132 -fbounds-check -fbacktrace
endif
#----------------------------------------------------------------
#  Library paths
#
BLAS = /usr/lib/libblas.so.3gf
LAPACK = /usr/lib/liblapack.so.3gf
#-------------------------------------------------------
#  Direcories where to search for files to compile to .o by implicit rules below, and dependencies defined in make.incdep
#
vpath %.o $(obsdir)
vpath %.f90 ../src
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.f90
	$(F95) -c $(FFLAGS) $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.PHONY:
#----------------------------------------------------------------
#  include dependencies:
#  make.incdep is a Makefile because it is included. It containes all dependencies of
#  the .o files. If you change any such dependencies (e.g. by using an additional module
#  in some program/module), please update file make.incdep accordingly. 
-include make.incdep
-include ../make.incdep
#---------------------------------------------------------------
clean:
	if (! -e $(bindir)) mkdir -p $(bindir)
	if (! -e $(obsdir)) mkdir -p $(obsdir)
	if (! -e $(moduledir)) mkdir -p $(moduledir)
	-rm -f $(bindir)/*
	-rm -f $(obsdir)/*.o
	-rm -f $(moduledir)/*.mod
#----------------------------------------------------------------
#
createSpecfem3dSyntheticData: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o \
	seismicStation.o asciiDataIO.o seismicEventList.o inversionBasics.o discreteFourierTransform.o \
	componentTransformation.o mathConstants.o seismicEvent.o commandLine.o realloc.o invgridVtkFile.o \
	eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o \
	wpVtkFile.o inputParameter.o integrationWeights.o flexibleType.o dateTime.o readEventStationFile.o \
	modelParametrization.o specfem3dKernelReferenceModel.o specfem3dWavefieldPoints.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o vectorPointer.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o timeUtils.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

createSpecfem3dMeasuredData: %: %.o errorMessage.o seismicNetwork.o iterationStepBasics.o fileUnitHandler.o \
	seismicStation.o asciiDataIO.o seismicEventList.o inversionBasics.o discreteFourierTransform.o \
	componentTransformation.o mathConstants.o seismicEvent.o commandLine.o realloc.o invgridVtkFile.o \
	eventStationVtkFile.o kernelReferenceModel.o wavefieldPoints.o inversionGrid.o kernelInvertedModel.o \
	wpVtkFile.o inputParameter.o integrationWeights.o flexibleType.o dateTime.o readEventStationFile.o \
	modelParametrization.o specfem3dKernelReferenceModel.o specfem3dWavefieldPoints.o ecartInversionGrid.o \
	specfem3dInversionGrid.o scartInversionGrid.o vectorPointer.o primitiveTypeEncoding.o simpleString.o \
	kindDefinitions.o timeUtils.o
	$(F95) -o $(bindir)/$@ $(obstring) $(BLAS) $(LAPACK)

all: createSpecfem3dSyntheticData createSpecfem3dMeasuredData