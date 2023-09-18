# SUBDIRS = 

check_PROGRAMS = abc_problem learn_reg learn_reg2_L2 learn_reg2_noise

CLEANFILES = 

abc_problem_SOURCES = abc_problem.cc abc_problem.h abc_problem_solve.cc abc_problem_solve.h abc_problem_algo.cc abc_problem_algo.h abc_problem_grad.cc abc_problem_grad.h abc_problem_func.cc abc_problem_func.h abc_problem_export.cc abc_problem_export.h abc_problem_neural.cc abc_problem_neural.h

learn_reg_SOURCES = learn_reg.cc

learn_reg2_L2_SOURCES = learn_reg2_L2.cc

learn_reg2_noise_SOURCES = learn_reg2_noise.cc


MUMPS_LIBS = @MUMPS_LIBS@
# MUMPS_CFLAGS = @MUMPS_CFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src # $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) @SUPLDFLAGS@ @BOOST_LIBS@

#TESTS = shallow_water.pl

EXTRA_DIST = \
        abc_problem.param		\
        abc_problem2.param		\
        learn_reg.param                \
        learn_reg2_L2.param            \
        learn_reg2_noise.param         

LOG_COMPILER = perl
