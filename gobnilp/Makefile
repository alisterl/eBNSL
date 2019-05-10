#-----------------------------------------------------------------------
#-------------------------- GOBNILP MAKE FILE --------------------------
#-----------------------------------------------------------------------

# This file was created by editing the Makefile for the linear ordering
# example in SCIP

# Set the path the SCIP directory
SCIPDIR = scip

# Include default project Makefile from SCIP
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
