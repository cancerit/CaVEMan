# what flags you want to pass to the C compiler & linker
CFLAGS = -Wall -g -O3 -DCAVEMAN_VERSION='"${VERSION}"' -DTEST_REF='"$(TEST_REF)"' -rdynamic

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = caveman generateCavemanUMNormVCF
caveman_SOURCES = caveman.c file_tests.c List.c List_algos.c bam_access.c config_file_access.c \
                fai_access.c ignore_reg_access.c alg_bean.c split_access.c covs_access.c \
                cn_access.c genotype.c algos.c output.c setup.c \
                split.c mstep.c merge.c estep.c
generateCavemanUMNormVCF_SOURCES = generateCavemanVCFUnmatchedNormalPanel.c file_tests.c List.c List_algos.c bam_access.c \
                                config_file_access.c fai_access.c ignore_reg_access.c \
                                alg_bean.c split_access.c covs_access.c cn_access.c \
                                genotype.c algos.c output.c setup.c split.c \
                                mstep.c merge.c estep.c

caveman_LDADD=-lhts -lpthread -lz -lm -ldl -lrt
caveman_CPPFLAGS=-I./
generateCavemanUMNormVCF_LDADD=-lhts -lpthread -lz -lm -ldl -lrt
generateCavemanUMNormVCF_CPPFLAGS=-I./

# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run_all_tests
run_all_tests_SOURCES = ../tests/run_all_tests.c \
        ../tests/check_alg_bean_tests.c \
        ../tests/check_algos_tests.c \
        ../tests/check_bam_access_tests.c \
        ../tests/check_cn_access_tests.c \
        ../tests/check_config_file_access_tests.c \
        ../tests/check_covs_access_tests.c \
        ../tests/check_fai_access_tests.c \
        ../tests/check_genotype_tests.c \
        ../tests/check_ign_region_access_tests.c \
        ../tests/check_output_tests.c \
        ../tests/check_split_access_tests.c \
        file_tests.c List.c List_algos.c bam_access.c \
        config_file_access.c fai_access.c ignore_reg_access.c \
        alg_bean.c split_access.c covs_access.c cn_access.c \
        genotype.c algos.c output.c setup.c split.c \
        mstep.c merge.c estep.c
run_all_tests_CFLAGS = -I../tests $(PTHREAD_CFLAGS)
run_all_tests_LDADD = -lcheck -lhts -lpthread -lz -lm -ldl $(PTHREAD_LDFLAGS)

if HOST_LINUX
run_all_tests_LDADD += -lrt -lsubunit
endif

