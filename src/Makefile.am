ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = tu
lib_LIBRARIES = libtotalunimodularity.a

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
	vector_three_connectivity.hpp

libtotalunimodularity_includedir = $(includedir)/tu/

libtotalunimodularity_include_HEADERS = \
	matroid_graph.hpp \
	total_unimodularity.hpp \
	matroid_decomposition.hpp

libtotalunimodularity_a_SOURCES = \
	determinant.cpp \
	ghouila_houri.cpp \
	zero_plus_minus_one.cpp \
	separation.cpp \
	logger.cpp \
	nested_minor_sequence.cpp \
	total_unimodularity.cpp \
	matroid_graph.cpp \
	matroid_decomposition.cpp

tu_SOURCES = \
	tu_main.cpp

tu_LDADD = \
  libtotalunimodularity.a

