#
# $Id: GNUmakefile,v 1.1.1.1 
#
# kinga,josh,echurch@fnal.gov
#
#
include SoftRelTools/arch_spec_root.mk

PACKAGE     := Genfit
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.py)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk

override CXXFLAGS += -Wno-overloaded-virtual -Wno-narrowing -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -Wno-reorder
override LIBLIBS += $(LOADLIBES) -lVMC 
