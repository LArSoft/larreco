#
# $Id: GNUmakefile,v 1.4 2009/06/05 18:49:52 t962cvs Exp $
#
# messier@indiana.edu
#
include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

override LIBLIBS += $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -lFilters -lRecoAlg -lUtilities -lGeometry -lRecoBase -lMCCheater
