CAVEMAN_VERSION=1.5.0

#Compiler
CC = gcc -O3 -DCAVEMAN_VERSION='"$(CAVEMAN_VERSION)"' -g

#CC = gcc -O3 -DCAVEMAN_VERSION='"$(CAVEMAN_VERSION)"' -g

#compiler flags
# -g adds debug info to the executable file
# -Wall turns on most warnings from compiler
CFLAGS = -Wall

HTSLOC?=$(HTSLIB)
SAMTOOLSLOC?=$(SAMTOOLS)

HTSTMP?=./caveman_tmp

#Define locations of header files
OPTINC?= -I$(SAMTOOLSLOC)/ -I$(HTSLOC)/ -I$(HTSLOC)/htslib
INCLUDES= -Isrc $(OPTINC) -rdynamic

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS?= -L$(SAMTOOLSLOC) -L$(HTSTMP)

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#   option, something like (this will link in libmylib.so and libm.so:
LIBS =-lhts -lbam -lpthread -lz -lm

# define the C source files
SRCS = ./src/file_tests.c ./src/List.c ./src/List_algos.c ./src/bam_access.c ./src/config_file_access.c ./src/fai_access.c ./src/ignore_reg_access.c ./src/alg_bean.c ./src/split_access.c ./src/covs_access.c ./src/cn_access.c ./src/genotype.c ./src/algos.c ./src/output.c ./src/setup.c ./src/split.c ./src/mstep.c ./src/merge.c ./src/estep.c
#Define test sources
TEST_SRC=$(wildcard ./tests/*_tests.c)
TESTS=$(patsubst %.c,%,$(TEST_SRC))

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

MD := mkdir

#Build target executable
CAVEMAN_TARGET=./bin/caveman
UMNORM_TARGET=./bin/generateCavemanUMNormVCF

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean coverage copyscript test make_htslib_tmp remove_htslib_tmp

.NOTPARALLEL: test

all: clean make_bin make_htslib_tmp $(CAVEMAN_TARGET) $(UMNORM_TARGET) copyscript test remove_htslib_tmp
	@echo  Binaries have been compiled.

$(UMNORM_TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(UMNORM_TARGET) $(OBJS) $(LFLAGS) $(LIBS) ./src/generateCavemanVCFUnmatchedNormalPanel.c

$(CAVEMAN_TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(CAVEMAN_TARGET) $(OBJS) $(LFLAGS) $(LIBS) ./src/caveman.c

#Unit Tests
test: $(CAVEMAN_TARGET)
test: CFLAGS += $(INCLUDES) $(OBJS) $(LFLAGS) $(LIBS)
test: $(TESTS)
	sh ./tests/runtests.sh

#Unit tests with coverage
coverage: CFLAGS += --coverage
coverage: test

make_bin:
	$(MD) ./bin

make_htslib_tmp:
	$(MD) $(HTSTMP)
	#Do some magic to ensure we compile CaVEMan with the static libhts.a rather than libhts.so
	ln -s $(HTSLOC)/libhts.a $(HTSTMP)/libhts.a

remove_htslib_tmp:
	rm -rf $(HTSTMP)

copyscript:
	rsync -uE ./scripts/* ./bin/
	chmod u+x ./bin/setupCaveman ./bin/splitCaveman ./bin/mstepCaveman ./bin/mergeCaveman ./bin/estepCaveman ./bin/mergeCavemanResults

valgrind:
	VALGRIND="valgrind --log-file=/tmp/valgrind-%p.log" $(MAKE)


# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) ./src/*.o *~ $(CAVEMAN_TARGET) $(UMNORM_TARGET) ./bin/* ./tests/tests_log $(TESTS) ./src/*.gcda ./src/*.gcov ./src/*.gcno *.gcda *.gcov *.gcno ./tests/*.gcda ./tests/*.gcov ./tests/*.gcno
	rm -rf ./bin

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
