#
# $Id: GNUmakefile,v 1.1.1.1 2009/05/13 01:38:51 katori Exp $
#
# kinga,josh,echurch@fnal.gov
#
#
include SoftRelTools/arch_spec_root.mk

PACKAGE     := EventFinder
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

#LIBLINK    := $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

override LIBLIBS += -L$(NUTOOLS_LIB) -lSimulationBase -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -lMCCheater -lGeometry -lSimulation -lRecoBase -lRawData -lUtilities
