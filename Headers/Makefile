#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the Headers subdirectory)
#
# !DESCRIPTION: This makefile compiles the various GEOS-Chem Header modules,
#  which contain many PARAMETERs and global arrays for GEOS-Chem routines.
#
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the directory   %%%
# %%% just above this one!                                                 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the GEOS-Chem code
# HDR        Specifies the directory where GEOS-Chem include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
#
# !REVISION HISTORY: 
#  23 Aug 2011 - M. Long     - Initial version
#  19 Mar 2012 - M. Payer    - Added EF_MGN20_mod for SOA + semivolatile POA
#                              simulation (H. Pye)
#  09 Apr 2012 - R. Yantosca - Removed CMN_VEL_mod.F; it's obsolete
#  19 Oct 2012 - R. Yantosca - Add modules for Grid-Independent GEOS-Chem
#  01 Nov 2012 - R. Yantosca - Added gigc_input_opt_mod.F90
#  16 Nov 2012 - R. Yantosca - Added more GIGC updates, removed obsolete 
#  15 Jan 2013 - R. Yantosca - Removed CMN_DEP_mod.F, it's obsolete
#  15 Jan 2013 - R. Yantosca - Added dependency for gigc_input_opt_mod.o
#                              to comode_loop_mod.F, to make it compile
#  19 Mar 2014 - R. Yantosca - Add more visible comment section dividers
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
ROOTDIR = ..
HDR     = $(ROOTDIR)/Headers
HELP    = $(ROOTDIR)/help
LIB     = $(ROOTDIR)/lib
MOD     = $(ROOTDIR)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

# List of source files
SRC = $(wildcard *.F) $(wildcard *.F90)

# Replace .f and .f90 extensions with *.o
TMP = $(SRC:.F=.o)
OBJ = $(TMP:.F90=.o)

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean help

lib: $(OBJ)
	$(AR) crs libHeaders.a $(OBJ)
	mv libHeaders.a $(LIB)

clean:
	rm -f *.o *.mod

help:
	@$(MAKE) -C $(HELP)

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

CMN_DIAG_mod.o       : CMN_DIAG_mod.F           \
                       CMN_SIZE_mod.o           \
                       CMN_FJX_MOD.o            \
                       gigc_errcode_mod.o

CMN_GCTM_mod.o       : CMN_GCTM_mod.F

CMN_NOX_mod.o        : CMN_NOX_mod.F            \
                       CMN_SIZE_mod.o           \
                       gigc_errcode_mod.o

CMN_O3_mod.o         : CMN_O3_mod.F             \
                       CMN_SIZE_mod.o           \
                       gigc_errcode_mod.o

CMN_SIZE_mod.o       : CMN_SIZE_mod.F           \
                       gigc_errcode_mod.o

CMN_mod.o            : CMN_mod.F                \
                       CMN_SIZE_mod.o           \
                       gigc_errcode_mod.o

CMN_FJX_MOD.o        : CMN_FJX_MOD.F            \
                       CMN_SIZE_mod.o

commsoil_mod.o       : commsoil_mod.F           \
                       CMN_SIZE_mod.o           \
                       gigc_errcode_mod.o

comode_loop_mod.o    : comode_loop_mod.F        \
                       CMN_SIZE_mod.o           \
                       gigc_errcode_mod.o       \
                       gigc_input_opt_mod.o

EF_MGN20_mod.o       : EF_MGN20_mod.F

gigc_errcode_mod.o   : gigc_errcode_mod.F90

gigc_input_opt_mod.o : gigc_input_opt_mod.F90   \
                       gigc_errcode_mod.o       \
                       CMN_SIZE_mod.o


gigc_state_chm_mod.o : gigc_state_chm_mod.F90   \
                       gigc_errcode_mod.o       \
                       gigc_input_opt_mod.o

gigc_state_met_mod.o : gigc_state_met_mod.F90   \
                       gigc_errcode_mod.o       \
                       CMN_SIZE_mod.o

gigc_state_phy_mod.o : gigc_state_phy_mod.F90   \
                       gigc_errcode_mod.o

smv_dimension_mod.o  : smv_dimension_mod.F

smv_physconst_mod.o  : smv_physconst_mod.F

#EOC
