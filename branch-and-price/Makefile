#@file    Makefile
#@brief   Makefile for k-vertex cut SCIP project
#@author  Fabio Ciccarelli

#-----------------------------------------------------------------------------
# path (MODIFY IF YOU CHANGE MACHINE)
#-----------------------------------------------------------------------------

SCIPDIR         =       $(SCIP_PATH)

#-----------------------------------------------------------------------------
# Gurobi configuration
#-----------------------------------------------------------------------------

GUROBIDIR       =       $(GUROBI_HOME)
GUROBILIBDIR    =       $(GUROBIDIR)/lib
GUROBIINCDIR    =       $(GUROBIDIR)/include

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP (need to do this twice, once to
# find the correct binary, then, after getting the correct flags from the
# binary (which is necessary since the ZIMPL flags differ from the default
# if compiled with the SCIP Optsuite instead of SCIP), we need to set the
# compile flags, e.g., for the ZIMPL library, which is again done in make.project
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project
SCIPVERSION			:=$(shell $(SCIPDIR)/bin/scip.$(BASE).$(LPS).$(TPI)$(EXEEXTENSION) -v | sed -e 's/$$/@/')
override ARCH		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ARCH=\([^@]*\).*/\1/')
override EXPRINT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* EXPRINT=\([^@]*\).*/\1/')
override GMP		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GMP=\([^@]*\).*/\1/')
override SYM		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SYM=\([^@]*\).*/\1/')
override IPOPT		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPT=\([^@]*\).*/\1/')
override IPOPTOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPTOPT=\([^@]*\).*/\1/')
override LPSCHECK	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSCHECK=\([^@]*\).*/\1/')
override LPSOPT 	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSOPT=\([^@]*\).*/\1/')
override NOBLKBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKBUFMEM=\([^@]*\).*/\1/')
override NOBLKMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKMEM=\([^@]*\).*/\1/')
override NOBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBUFMEM=\([^@]*\).*/\1/')
override THREADSAFE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* THREADSAFE=\([^@]*\).*/\1/')
override READLINE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* READLINE=\([^@]*\).*/\1/')
override SANITIZE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SANITIZE=\([^@]*\).*/\1/')
override ZIMPL		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPL=\([^@]*\).*/\1/')
override ZIMPLOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPLOPT=\([^@]*\).*/\1/')
override ZLIB		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZLIB=\([^@]*\).*/\1/')
include $(SCIPDIR)/make/make.project


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	kvertexcut
CMAINOBJ	=	
CXXMAINOBJ	=	main.o \
		        probdata_kvertexcut.o \
		        reader_kvertexcut.o \
		        vardata_kvertexcut.o \
		        pricer_kvertexcut.o \
			graph.o \
			preprocessing.o \
			min_cut_pricer.o \
			plugins.o \
			custom_compute_symmetry.o \
			prop_problemsymmetry.o 

MAINSRC		=	$(addprefix $(SRCDIR)/,$(CMAINOBJ:.o=.c))
MAINSRC		+=	$(addprefix $(SRCDIR)/,$(CXXMAINOBJ:.o=.cpp))

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(CMAINOBJ))
MAINOBJFILES	+=	$(addprefix $(OBJDIR)/,$(CXXMAINOBJ))

#-----------------------------------------------------------------------------
# External libraries
#-----------------------------------------------------------------------------

FLAGS		+=	-I$(GUROBIINCDIR)
CXXFLAGS	+=	-std=c++17
LDFLAGS		+=	-L$(GUROBILIBDIR) -lgurobi_c++ -lgurobi120 -lm -Wl,-rpath,$(GUROBILIBDIR)

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: run
run:		$(MAINFILE)
		@echo "-> running $(MAINFILE) with Gurobi libraries"
		LD_LIBRARY_PATH=$(GUROBILIBDIR):$$LD_LIBRARY_PATH $(MAINFILE)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) -I$(SCIPDIR) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -USCIP_WITH_READLINE -USCIP_ROUNDING_FE $$i; \
			done'
.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-rm -rf $(OBJDIR)
		@echo "-> remove main objective files"
endif
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
		@echo "-> remove binary"

.PHONY: tags
tags:
	rm -f TAGS; ctags -e src/*.cpp src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

-include	$(MAINOBJFILES:.o=.d)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) $(LINKCXXSCIPALL) $(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		@mkdir -p $(dir $@)
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) $(DFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		@mkdir -p $(dir $@)
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) $(DFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
