#
# $Id: GNUmakefile,v 1.1.1.1 2009/05/13 01:38:51 katori Exp $
#
# kinga,josh,echurch@fnal.gov
#
#
include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

# needed for ClusterCrawlerAlg
ROOTLIBS = $(shell root-config --libs) -lMinuit
override LOADLIBES += $(ROOTLIBS)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -lRecoBase -lGeometry -lUtilities -lMCCheater -lSimulation -lRecoObjects
