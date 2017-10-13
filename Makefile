#############################################################
# Generic Makefile for C/C++ Program
#
# Description:
# ------------
# This is an easily customizable makefile template. The purpose is to
# provide an instant building environment for C/C++ programs.
#
# It searches all the C/C++ source files in the specified directories,
# makes dependencies, compiles and links to form an executable.
#
# Besides its default ability to build C/C++ programs which use only
# standard C/C++ libraries, you can customize the Makefile to build
# those using other libraries. Once done, without any changes you can
# then build programs using the same or less libraries, even if source
# files are renamed, added or removed. Therefore, it is particularly
# convenient to use it to build codes for experimental or study use.
#
# GNU make is expected to use the Makefile. Other versions of makes
# may or may not work.
#
# Usage:
# ------
# 1. Copy the Makefile to your program directory.
# 2. Customize in the "Customizable Section" only if necessary:
#    * to use non-standard C/C++ libraries, set pre-processor or compiler
#      options to <MY_CFLAGS> and linker ones to <MY_LIBS>
#      (See Makefile.gtk+-2.0 for an example)
#    * to search sources in more directories, set to <SRCDIRS>
#    * to specify your favorite program name, set to <PROGRAM>
# 3. Type make to start building your program.
#
# Make Target:
# ------------
# The Makefile provides the following targets to make:
# http://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
#   $ make           compile and link
#   $ make NODEP=yes compile and link without generating dependencies
#   $ make objs      compile only (no linking)
#   $ make tags      create tags for Emacs editor
#   $ make ctags     create ctags for VI editor
#   $ make clean     clean objects and the executable file
#   $ make distclean clean objects, the executable and dependencies
#   $ make help      get the usage of the makefile
#   $ make run		 execute the main program
#
#===========================================================================

## Customizable Section: adapt those variables to suit your program.
##==========================================================================

# The pre-processor and compiler options.
# -g - turn on debugging (so GDB gives more friendly output)
# -Wall - turns on most warnings
# -O or -O2 - turn on optimizations
# -Ofast Disregard strict standards compliance. -Ofast enables all -O3 optimizations.
# It also enables optimizations that are not valid for all standard-compliant programs.
# It turns on -ffast-math and the Fortran-specific -fno-protect-parens and -fstack-arrays.
# -o <name> - name of the output file
# -c - output an object file (.o)
# -I<include path> - specify an include directory
# -L<library path> - specify a lib directory
# -l<library> - link with library lib<library>.a

# Cantera makefile template
# include /usr/include/cantera/Cantera.mak
include /opt/cantera/v2.4.0/include/cantera/Cantera.mak

#num of processors
P=3

MY_CFLAGS =

#MPI_ROOT_DIR, mpirun, mpiextc, $(MPI_ROOT_DIR)mpic++ directory
MPI_ROOT_DIR= /usr/bin/

# The linker options.
BOOST_ROOT_INCLUDE= -I/usr/local/include
BOOST_ROOT_LIB= -L/usr/local/lib

OPENMPI_ROOT_INCLUDE= -I/usr/lib/openmpi/include
OPENMPI_ROOT_LIB= -L/usr/lib/openmpi/lib

BOOST_MPI_OPTIONS= -lboost_mpi -lboost_serialization -lboost_program_options -lboost_filesystem -lboost_random -lboost_system -lboost_regex

BOOST_REGEX_DIR= ./include/regex/src

# The link options, directory and libs
MY_LIBS   = $(BOOST_ROOT_LIB) $(OPENMPI_ROOT_LIB) $(CANTERA_LIBS) -lgfortran $(BOOST_MPI_OPTIONS)

# The pre-processor options used by the cpp (man cpp for more)
CPPFLAGS  = -Wall -fopenmp $(BOOST_ROOT_INCLUDE) $(OPENMPI_ROOT_INCLUDE) $(CANTERA_INCLUDES) -std=c++11

# The options used in linking as well as in any direct use of ld.
LDFLAGS   = 

# The directories in which source files reside.
# If not specified, only the current directory will be serached.

# Fortran Source directory
FSRC_ROOT_DIR = ./fortran_lib
FSRCDIRS   = $(FSRC_ROOT_DIR) $(FSRC_ROOT_DIR)/dlsode $(FSRC_ROOT_DIR)/chemkin $(FSRC_ROOT_DIR)/cantera

# CXX Source directory
SRCDIRS   = . ./src ./src/drivers ./src/tools/debug ./src/relationshipParser ./src/srkin ./src/time ./src/mechanism ./src/odeSolver ./src/statistics ./src/cubicSpline ./src/fileIO/commandLineConfigFileReader ./src/fileIO/CSVFileReader ./src/fileIO/fileIO ./src/pathwayHandler ./src/propagator/superPropagator ./src/propagator/dlsodePropagator ./src/propagator/ssaPropagator ./src/propagator/SOHRPropagator ./src/tools/map_reduce ./src/tools/block_decomposition ./src/tools/sr_libs ./src/tools/matrix ./src/tools/union_find ./src/tools/chattering ./src/random ./src/reactionNetwork/superReactionNetwork ./src/reactionNetwork/concreteReactionNetwork ./src/reactionNetwork/reactionNetworkODESolver ./src/search_algorithm/eppstein_algorithm $(FSRC_ROOT_DIR)/cantera $(BOOST_REGEX_DIR)

# The executable file name.
# If not specified, current directory name or `a.out' will be used.
PROGRAM   = RxnNetwork_PI

## Implicit Section: change the following only when necessary.
##==========================================================================

# The source file types (headers excluded).
# .c indicates C source files, and others C++ ones.
SRCEXTS = .c .C .cc .cpp .CPP .c++ .cxx .cp
# .f indicates Fortran source files
FSRCEXTS = .f .f90 .f95


# The header file types.
HDREXTS = .h .H .hh .hpp .HPP .h++ .hxx .hp

# The pre-processor and compiler options.
# Users can override those variables from the command line.
CFLAGS  = -g -Ofast
CXXFLAGS= -g -Ofast
# -w: ignore waring, add cantera libs
FFLAGS	= -g -cpp -ffpe-trap=invalid,zero,overflow,underflow -w -Ofast $(CANTERA_INCLUDES) $(CANTERA_FORTRAN_LIBS)

# The C program compiler.
#CC     = gcc
#CC		= icc
CC		= $(MPI_ROOT_DIR)mpic++
# The C++ program compiler.
#CXX    = g++
#CXX		= icpc
CXX		= $(MPI_ROOT_DIR)mpic++
# Un-comment the following line to compile C programs as C++ ones.
#CC     = $(CXX)
# The Fortran compiler
#FC	= g++
#FC	= ifort
FC	= $(MPI_ROOT_DIR)mpic++
# The Fortran C/C++ mixed compiler
#FCC	= g++
#FCC	= ifort
FCC	= $(MPI_ROOT_DIR)mpic++

# The command used to delete file.
RM     = rm -f

ETAGS = etags
ETAGSFLAGS =

CTAGS = ctags
CTAGSFLAGS =

## Stable Section: usually no need to be changed. But you can add more.
##==========================================================================
SHELL   = /bin/sh
EMPTY   =
SPACE   = $(EMPTY) $(EMPTY)
ifeq ($(PROGRAM),)
  CUR_PATH_NAMES = $(subst /,$(SPACE),$(subst $(SPACE),_,$(CURDIR)))
  PROGRAM = $(word $(words $(CUR_PATH_NAMES)),$(CUR_PATH_NAMES))
  ifeq ($(PROGRAM),)
    PROGRAM = a.out
  endif
endif
ifeq ($(SRCDIRS),)
  SRCDIRS = .
endif

SOURCES = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))
HEADERS = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(HDREXTS))))
SRC_CXX = $(filter-out %.c,$(SOURCES))

# Don't compile ckstrt.f
FSOURCES_ALL = $(foreach d,$(FSRCDIRS),$(wildcard $(addprefix $(d)/*,$(FSRCEXTS))))
FSOURCES = $(filter-out %ckstrt.f %ckvariables.f %opkddemos.f %math.f %test.f, $(FSOURCES_ALL))


OBJS    = $(addsuffix .o, $(basename $(SOURCES)))
FOBJS    = $(addsuffix .o, $(basename $(FSOURCES)))
DEPS    = $(OBJS:.o=.d)

## Define some useful variables.
DEP_OPT = $(shell if `$(CC) --version | grep "GCC" >/dev/null`; then \
                  echo "-MM -MP"; else echo "-M"; fi )
DEPEND      = $(CC)  $(DEP_OPT)  $(MY_CFLAGS) $(CFLAGS) $(CPPFLAGS)
DEPEND.d    = $(subst -g ,,$(DEPEND))
COMPILE.c   = $(CC)  $(MY_CFLAGS) $(CFLAGS)   $(CPPFLAGS) -c
COMPILE.cxx = $(CXX) $(MY_CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -c
LINK.c      = $(CC)  $(MY_CFLAGS) $(CFLAGS)   $(CPPFLAGS) $(LDFLAGS)
LINK.cxx    = $(CXX) $(MY_CFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS)
COMPILE.f   = $(FC)  $(FFLAGS)
COMPILE.f.c   = $(FCC)  $(FCFLAGS)

.PHONY: all objs tags ctags clean distclean help show

# Delete the default suffixes
.SUFFIXES:

all: $(PROGRAM)

# Rules for creating dependency files (.d).
#------------------------------------------

%.d:%.c
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

%.d:%.C
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

%.d:%.cc
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

%.d:%.cpp
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

%.d:%.CPP
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

%.d:%.c++
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

%.d:%.cp
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

%.d:%.cxx
	@echo -n $(dir $<) > $@
	@$(DEPEND.d) $< >> $@

# Rules for generating object files (.o).
#----------------------------------------
objs:$(OBJS)

%.o:%.c
	$(COMPILE.c) $< -o $@

%.o:%.C
	$(COMPILE.cxx) $< -o $@

%.o:%.cc
	$(COMPILE.cxx) $< -o $@

%.o:%.cpp
	$(COMPILE.cxx) $< -o $@

%.o:%.CPP
	$(COMPILE.cxx) $< -o $@

%.o:%.c++
	$(COMPILE.cxx) $< -o $@

%.o:%.cp
	$(COMPILE.cxx) $< -o $@

%.o:%.cxx
	$(COMPILE.cxx) $< -o $@

# Rules for generating fortran object files (.o).
#----------------------------------------
fobjs:$(FOBJS)

%.o:%.f
	$(COMPILE.f) -c -o $@ $<

%.o:%.f90
	$(COMPILE.f) -c -o $@ $<

%.o:%.f95
	$(COMPILE.f) -c -o $@ $<


# Rules for generating the tags.
#-------------------------------------
tags: $(HEADERS) $(SOURCES)
	$(ETAGS) $(ETAGSFLAGS) $(HEADERS) $(SOURCES)

ctags: $(HEADERS) $(SOURCES)
	$(CTAGS) $(CTAGSFLAGS) $(HEADERS) $(SOURCES)

# Rules for generating the executable.
#-------------------------------------
$(PROGRAM):$(OBJS) $(FOBJS)
ifeq ($(SRC_CXX),)              # C program
	$(LINK.c)   $(OBJS) $(MY_LIBS) -o $@
	@echo Type ./$@ to execute the program.
else                            # C++ program # The link options have to be at the end of the command
	$(LINK.cxx) $(OBJS) $(FOBJS) $(MY_LIBS) -o $@
	@echo Type ./$@ to execute the program.
endif

ifndef NODEP
ifneq ($(DEPS),)
  sinclude $(DEPS)
endif
endif

clean:
	$(RM) $(FCOBJS) $(FOBJS) $(OBJS) $(PROGRAM) $(PROGRAM).exe

distclean: clean
	$(RM) $(DEPS) TAGS

# Show help.
help:
	@echo 'Generic Makefile for C/C++ Programs (gcmakefile) version 0.5'
	@echo 'Copyright (C) 2007, 2008 whyglinux <whyglinux@hotmail.com>'
	@echo 'Copyright (C) 2013, 2014 Shaun.Bai <Shaun.Bai@hotmail.com>'
	@echo
	@echo 'Usage: make [TARGET]'
	@echo 'TARGETS:'
	@echo '  all       (=make) compile and link.'
	@echo '  NODEP=yes make without generating dependencies.'
	@echo '  objs      compile only (no linking).'
	@echo '  tags      create tags for Emacs editor.'
	@echo '  ctags     create ctags for VI editor.'
	@echo '  clean     clean objects and the executable file.'
	@echo '  distclean clean objects, the executable and dependencies.'
	@echo '  show      show variables (for debug use only).'
	@echo '  run       execute the main program.'
	@echo '  help      print this message.'
	@echo
	@echo 'Report bugs to <bunnysirah AT hotmail DOT com>.'

# Show variables (for debug use only.)
show:
	@echo 'PROGRAM     :' $(PROGRAM)
	@echo 'SRCDIRS     :' $(SRCDIRS)
	@echo 'HEADERS     :' $(HEADERS)
	@echo 'SOURCES     :' $(SOURCES)
	@echo 'SRC_CXX     :' $(SRC_CXX)
	@echo 'OBJS        :' $(OBJS)
	@echo 'DEPS        :' $(DEPS)
	@echo 'DEPEND      :' $(DEPEND)
	@echo 'COMPILE.c   :' $(COMPILE.c)
	@echo 'COMPILE.cxx :' $(COMPILE.cxx)
	@echo 'link.c      :' $(LINK.c)
	@echo 'link.cxx    :' $(LINK.cxx)

# execute the main program
run:
	$(MPI_ROOT_DIR)mpirun -np $(P) ./$(PROGRAM) -c $(PWD)/input/cmd.cfg

## End of the Makefile ##  Suggestions are welcome  ## All rights reserved ##
##############################################################
