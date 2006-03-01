AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = mincfit
dist_man1_MANS = mincfit.man1

EXTRA_DIST = $(m4_files) \
	epm-header.in \
   $(man_MANS)

m4_files = m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	m4/mni_REQUIRE_OPENINVENTOR.m4 \
	m4/mni_cxx_have_koenig_lookup.m4 \
	m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_WITH_BUILD_PATH.m4

mincfit_SOURCES = mincfit.c \
	minc_vector_io.h minc_vector_io.c \
	gamma_func.h gamma_func.c \
	t1_func.h t1_func.c \
	common_func.h
