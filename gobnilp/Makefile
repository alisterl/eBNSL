#-----------------------------------------------------------------------
#-------------------------- GOBNILP MAKE FILE --------------------------
#-----------------------------------------------------------------------

# This file was created by editing the Makefile for the linear ordering
# example in SCIP

# Set the path to the SCIP directory
SCIPDIR = scip

# The following 'trick' from Stefan Vigerske's Makefile for the CallableLibrary example
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
override GAMS		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GAMS=\([^@]*\).*/\1/')
override GMP		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GMP=\([^@]*\).*/\1/')
override SYM		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SYM=\([^@]*\).*/\1/')
override IPOPT		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPT=\([^@]*\).*/\1/')
override IPOPTOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPTOPT=\([^@]*\).*/\1/')
override LPSCHECK	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSCHECK=\([^@]*\).*/\1/')
override LPSOPT 	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSOPT=\([^@]*\).*/\1/')
override NOBLKBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKBUFMEM=\([^@]*\).*/\1/')
override NOBLKMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKMEM=\([^@]*\).*/\1/')
override NOBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBUFMEM=\([^@]*\).*/\1/')
override PARASCIP	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* PARASCIP=\([^@]*\).*/\1/')
override READLINE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* READLINE=\([^@]*\).*/\1/')
override SANITIZE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SANITIZE=\([^@]*\).*/\1/')
override ZIMPL		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPL=\([^@]*\).*/\1/')
override ZIMPLOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPLOPT=\([^@]*\).*/\1/')
override ZLIB		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZLIB=\([^@]*\).*/\1/')
include $(SCIPDIR)/make/make.project


#OPEN BLAS not installed

ifdef BLASDIR
BLASINCDIR = $(BLASDIR)/include
BLASLIBDIR = $(BLASDIR)/lib
endif

ifdef BLASINCDIR
FLAGS += -I$(BLASINCDIR)
endif

ifdef BLAS
FLAGS += -DBLAS
endif

ifdef BLASLIBDIR
LINKBLAS = $(BLASLIBDIR)/libopenblas.a -lpthread
else
LINKBLAS =
endif

# Declare the phony targets
.PHONY: all doxygen manual pedmanual gobnilp clean test

# This stops the intermediate latex files being deleted
.SECONDARY:

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK) 
endif

#-----------------------------------------------------------------------
#----------------------- Command Line Arguments  -----------------------
#-----------------------------------------------------------------------

# HUSH is used to optionally hide output from latex and Doxygen
# (Turned on by default)
HUSH 				=	TRUE
ifeq ($(HUSH), TRUE)
	HUSHOUT = &>/dev/null
endif

#-----------------------------------------------------------------------
#---------------------------- Declarations  ----------------------------
#-----------------------------------------------------------------------

VERSION         =       development
MAINNAME	=	gobnilp
ALLMAINSRC	=	$(wildcard $(SRCDIR)/*.c)
BGESRC          =       $(wildcard $(SRCDIR)/bge_*.c)
ifdef BLASINCDIR
MAINSRC         =       $(ALLMAINSRC)
else
MAINSRC         =       $(filter-out $(BGESRC),$(ALLMAINSRC))
endif
MAINOBJ		=	$(notdir $(MAINSRC:.c=.o))

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

DOXYGENFILE		=	Doxyfile
DOXYGENDIR		=	docs

MANUALNAME		= 	manual
PEDMANUALNAME	=	pedigree_manual
LATEXCMD			=	@pdflatex -halt-on-error -interaction=errorstopmode

#Needs to appear first so that "make" makes gobnilp
gobnilp: $(SCIP) $(MAINFILE) $(MAINSHORTLINK)


#-----------------------------------------------------------------------
#----------------------------      All      ----------------------------
#-----------------------------------------------------------------------

all: gobnilp doxygen manual pedmanual


#-----------------------------------------------------------------------
#----------------------------    Doxygen    ----------------------------
#-----------------------------------------------------------------------

doxygen: $(DOXYGENDIR)

$(DOXYGENDIR): $(DOXYGENFILE) $(MAINSRC)
		@echo "Making Doxygen files"
		@doxygen $(HUSHOUT)

#-----------------------------------------------------------------------
#----------------------------    Manuals    ----------------------------
#-----------------------------------------------------------------------

manual: $(MANUALNAME).pdf
pedmanual: $(PEDMANUALNAME).pdf

%.pdf: %.ltx %.bbl
		@echo "Making manual"
		@$(LATEXCMD) $*.ltx $(HUSHOUT)
		@$(LATEXCMD) $*.ltx $(HUSHOUT)
		@$(LATEXCMD) $*.ltx $(HUSHOUT)
		@rm -f $*.aux $*.log $*.toc

%.bbl: %.bib %_tmp.aux
		@bibtex $*_tmp $(HUSHOUT)
		@mv $*_tmp.bbl $*.bbl
		@rm -f $*_tmp.blg

%_tmp.aux: %.ltx
		@$(LATEXCMD) -jobname $*_tmp $*.ltx $(HUSHOUT)
		@rm -f $*_tmp.log $*_tmp.pdf $*_tmp.toc



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
		@-(rm -f $(OBJDIR)/*.o && rmdir $(OBJDIR));
		@echo "Removing object files"
endif
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
		@echo "Removing binaries"
		@echo "Removing documentation"
		@if [ -d $(DOXYGENDIR) ] ; then                                  \
			if [ -d $(DOXYGENDIR)/html ] ; then                           \
				if [ -d $(DOXYGENDIR)/html/search ] ; then                 \
					rm -f $(DOXYGENDIR)/html/search/* ;                     \
					rmdir $(DOXYGENDIR)/html/search ;                       \
				fi ;                                                       \
				rm -f $(DOXYGENDIR)/html/* ;                               \
				rmdir $(DOXYGENDIR)/html ;                                 \
			fi ;                                                          \
			rmdir $(DOXYGENDIR) ;                                         \
		fi

		@echo "Removing manual"
		@rm -f $(MANUALNAME)_tmp.aux
		@rm -f $(MANUALNAME).bbl
		@rm -f $(MANUALNAME).pdf
		@rm -f $(PEDMANUALNAME)_tmp.aux
		@rm -f $(PEDMANUALNAME).bbl
		@rm -f $(PEDMANUALNAME).pdf

#-----------------------------------------------------------------------
#----------------------------     Test      ----------------------------
#-----------------------------------------------------------------------

test:           $(MAINFILE)
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/getlastprob.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_set.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_logfiles.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_tmpfile_setup_scip.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/run.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck_cluster.sh);
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(EXECUTABLE) $(SETTINGS) $(BINID) $(OUTPUTDIR) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) \
		$(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(DEBUGTOOL) $(CLIENTTMPDIR) $(REOPT) $(OPTCOMMAND) $(SETCUTOFF) $(MAXJOBS) $(VISUALIZE) $(PERMUTE) \
                $(SEEDS) $(GLBSEEDSHIFT);

#-----------------------------------------------------------------------
#----------------------------    GOBNILP    ----------------------------
#-----------------------------------------------------------------------


# main target
$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) $(LINKCXXSCIPALL) $(LDFLAGS) \
		$(LINKBLAS) \
		$(LINKCXX_o)$@



$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c -g $< $(CC_o)$@

#---- EOF --------------------------------------------------------------------
