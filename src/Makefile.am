ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = unimodularity-test
lib_LIBRARIES = libunimodularity.a
check_PROGRAMS = tu-compare matrix-gen

noinst_HEADERS = \
  algorithm.hpp \
	binary_linear_space.hpp \
	bipartite_graph_bfs.hpp \
	logger.hpp \
	comparators.hpp \
	enumeration.hpp \
	find_minor_sequence.hpp \
	find_wheel_minor.hpp \
	graphicness.hpp \
	graph_utils.hpp \
	matrix.hpp \
	matrix_modified.hpp \
	matrix_permuted.hpp \
	matrix_reorder.hpp \
	matrix_transposed.hpp \
	matrix_transposed_permuted.hpp \
	matroid.hpp \
	matroid_permuted.hpp \
	matroid_reorder.hpp \
	matroid_separate.hpp \
	matroid_transposed.hpp \
	nested_minor_sequence.hpp \
	signing.hpp \
	partition.hpp \
	permutations.hpp \
	r10.hpp \
	separation.hpp \
	vector_three_connectivity.hpp \
	gen_generic.hpp \
  gen_random.hpp \
  gen_network.hpp \
  gen_cycle_violator.hpp \
  combinations.hpp \
  smith_normal_form.hpp \
  linear_algebra.hpp \
  violator_search.hpp \
  gcd.hpp

libunimodularity_includedir = $(includedir)/tu/

libunimodularity_include_HEADERS = \
	common.hpp \
	matroid_graph.hpp \
	total_unimodularity.hpp \
	matroid_decomposition.hpp \
	unimodularity.hpp

libunimodularity_a_SOURCES = \
	determinant.cpp \
	ghouila_houri.cpp \
	zero_plus_minus_one.cpp \
	separation.cpp \
	logger.cpp \
	nested_minor_sequence.cpp \
	total_unimodularity.cpp \
	matroid_graph.cpp \
	matroid_decomposition.cpp \
	unimodularity.cpp

unimodularity_test_SOURCES = \
	unimodularity_test_main.cpp

unimodularity_test_LDADD = \
  libunimodularity.a

matrix_gen_SOURCES = \
  matrix_gen_main.cpp

matrix_gen_LDADD = \
  libunimodularity.a

tu_compare_SOURCES = \
	tu_compare_main.cpp

tu_compare_LDADD = \
  libunimodularity.a

