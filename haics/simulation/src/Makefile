SHELL = /bin/bash -f
HAICS_HOME = $(shell cd ..; pwd)

# input directory
SDIR = $(HAICS_HOME)/src
# output directory
ODIR =$(HAICS_HOME)/objects
# directory of executables
BDIR =$(HAICS_HOME)/bin

# CC compiler
COMP_PATH = '/usr/bin'
CC     = gcc
# remove -g for productive executables
CFLAGS = -m64 -O3 -Wall -Wno-long-long
CFLAGS += -DNDEBUG
#CFLAGS = -g -m64  -O1 -Wall -Werror
IFLAGS = -I$(SDIR) -I$(ODIR)
#IFLAGS += -I/usr/include
#IFLAGS += -I/usr/local/opt/readline/include
IFLAGS += -I/usr/local/include


#LDFLAGS = -L/usr/local/lib
LDFLAGS = -L/usr/lib
#LDFLAGS += -L/usr/lib/gcc/x86_64-linux-gnu/4.6.3/lto-wrapper
#LDFLAGS += -L/usr/lib/gcc/x86_64-linux-gnu/4.6
#LDFLAGS +=  -L/usr/local/opt/readline/lib
LDFLAGS += -lm
LDFLAGS += -lgsl
LDFLAGS += -llapack
LDFLAGS += -llapacke
#LDFLAGS +=  -lmkl_core
#LDFLAGS += -lopenblas
LDFLAGS += -lgslcblas
#LDFLAGS += -lnsl
LDFLAGS += -lsundials_cvodes
LDFLAGS += -lsundials_nvecserial
LDFLAGS += -lsundials_nvecmanyvector
#LDFLAGS += -L/usr/lib/lapack/liblapack.so
LDFLAGS += -L/usr/lib/x86_64-linux-gnu/libgsl.so
LDFLAGS += -L/usr/lib/x86_64-linux-gnu/libgslcblas.so
LDFLAGS += -L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so
LDFLAGS += -L/usr/lib/x86_64-linux-gnu/liblapacke.so
rm     = rm -f

# Name of the program
LINK_TARGET = haics

# Source files
SRC = 	adapt_integrator.c \
	Globals.c \
	hmc.c \
	integrators.c \
	main.c \
	read_input.c \
	seir_functions.c \
	seir_model_synthetic.c \
	sir_functions.c \
	sir_model_functions.c \
	sir_model_synthetic.c \
	sir_seir_common_functions.c \
	sir_stand_model_synthetic.c \
	sir_stand_model_prevalence_synthetic.c \
	sir_stand_functions.c \
	utils.c

# Objects
objs = $(SRC:.c=.o)
OBJS = $(patsubst %,$(ODIR)/%,$(objs))

# Non sources files
REBUILDABLES = $(OBJS) $(BDIR)/$(LINK_TARGET)

# The rule for "all" is used to incrementally build the system.
# It is the default entry point for make command
all : $(BDIR)/$(LINK_TARGET)
	@echo -e "\n\nAll done. To run the program: $(BDIR)/$(LINK_TARGET) <ID>\n"

# Dependency Rules
adapt_integrator.o : $(SDIR)/Globals.h Definitions.h
Globals.o : $(SDIR)/Globals.h
hmc.o : ${SDIR}/Globals.h integrators.h Definitions.h utils.h sir_seir_common_functions.h sir_model_functions.h sir_functions.h seir_functions.h sir_stand_model_synthetic.h sir_stand_functions.h
integrators.o : ${SDIR}/adapt_integrator.h Globals.h utils.h
main.o : $(SDIR)/read_input.h hmc.h Globals.h
read_input.o : ${SDIR}/Globals.h Definitions.h utils.h integrators.h adapt_integrator.h sir_model.h sir_model_synthetic.h seir_model_synthetic.h sir_stand_model_synthetic.h
seir_functions.o : ${SDIR}/Globals.h utils.h read_input.h sir_seir_common_functions.h
seir_model_synthetic.o : ${SDIR}/seir_model_synthetic.h Globals.h utils.h seir_functions.h read_input.h sir_seir_common_functions.h
sir_functions.o : ${SDIR}/Globals.h utils.h sir_seir_common_functions.h read_input.h
sir_model_functions.o : ${SDIR}/Globals.h utils.h read_input.h sir_seir_common_functions.h
sir_model_synthetic.o : ${SDIR}/sir_seir_common_functions.h Globals.h utils.h sir_functions.h read_input.h
sir_stand_model_synthetic.o : ${SDIR}/sir_seir_common_functions.h Globals.h utils.h sir_functions.h read_input.h sir_stand_functions.h
sir_stand_model_prevalence_synthetic.o : ${SDIR}/sir_seir_common_functions.h Globals.h utils.h sir_functions.h read_input.h sir_stand_functions.h
sir_stand_functions.o : ${SDIR}/Globals.h utils.h read_input.h sir_seir_common_functions.h
sir_seir_common_functions.o : ${SDIR}/Globals.h utils.h read_input.h
utils.o : ${SDIR}/Definitions.h

# $@ for the pattern-matched target
# $< for the pattern-matched dependency
$(ODIR)/%.o : $(SDIR)/%.c
	$(COMP_PATH)/$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

# $@ expands to the rule's target, in this case "hmcsam".
# $^ expands to the rule's dependencies; the filenames of all the prerequisites, separated by spaces, discard duplicates.
$(BDIR)/$(LINK_TARGET) : $(OBJS)
	$(COMP_PATH)/$(CC) $(CFLAGS) $(IFLAGS) $^ -o $@ $(LDFLAGS)

# Headers
# HDRS = hmc.h Globals.h ranlxd.h read_input.h func_caller.h integrators.h sv_input.h

# cleaning your build environment
clean :
	$(rm) $(REBUILDABLES)
	@echo -e "\n\nCleaning done.\n"
