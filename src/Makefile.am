ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = unimodularity-test matrix-gen
lib_LIBRARIES = libunimodularity.a

noinst_headers = \
  algorithm.hpp \
	binary_linear_space.hpp \
	bipartite_graph_bfs.hpp \
	logger.hpp \
	comparators.hpp \
	enumeration.hpp \
	find_minor_sequence.hpp \
	find_wheel_minor.hpp \
	generator_options.hpp \
	graphicness.hpp \
	graph_utils.hpp \
	matrix.hpp \
	matrix_modified.hpp \
	matrix_permuted.hpp \
	matrix_reorder.hpp \
	matrix_transposed.hpp \
	matroid.hpp \
	matroid_permuted.hpp \
	matroid_reorder.hpp \
	matroid_separate.hpp \
	matroid_transposed.hpp \
	nested_minor_sequence.hpp \
	signing.cpp \
	partition.hpp \
	permutations.hpp \
	r10.hpp \
	row_column_index.hpp \
	separation.hpp \
	vector_three_connectivity.hpp \
	gen_generic.hpp \
  gen_random.hpp \
  gen_network.hpp \
  gen_cycle_violator.hpp \
  smith_normal_form.hpp

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
	unimodularity.cpp \
  smith_normal_form.cpp

unimodularity_test_SOURCES = \
	unimodularity_test_main.cpp

unimodularity_test_LDADD = \
  libunimodularity.a

matrix_gen_SOURCES = \
  matrix_gen_main.cpp

matrix_gen_LDADD = \
  libunimodularity.a
