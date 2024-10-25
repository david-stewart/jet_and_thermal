# The libraries for pythia, fastjet, and root ( "LIB_TRI" )
PY8PATH=/home/davidstewart/software/pythia8310
LIB_PYTH=-I${PY8PATH}/include -L${PY8PATH}/lib -lpythia8
LIB_FASTJET=`/home/davidstewart/software/fastjet-install/bin/fastjet-config --cxxflags --libs`
LIB_ROOT=`root-config --cflags --glibs`
LIB_TRI=${LIB_PYTH} ${LIB_FASTJET} ${LIB_ROOT} -lRecursiveTools

# compilation option
CC=g++
CFLAGS=-std=c++17 -O3 -Wno-deprecated

bin/IP_match_trees: \
	obj/IP_match_trees.o \
	obj/IP_geom_matcher.o \
	obj/BkgGen.o \
	obj/JetClusterer.o \
	obj/RhoMedianBkgClusterer.o \
	obj/JetTree.o \
	obj/JetBranch.o \
	obj/JetscapeHadronReader.o
	${CC} ${CFLAGS} -o $@ $^ ${LIB_TRI}

obj/IP_match_trees.o: src/IP_match_trees.cc src/BkgGen.h \
	src/P8Gen.h src/JetClusterer.h src/JetIndicesMatcher.h \
	src/IP_geom_matcher.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/JetscapeHadronReader.o: src/JetscapeHadronReader.cc \
	src/JetscapeHadronReader.h \
	src/JetClusterer.h 
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/BkgGen.o: src/BkgGen.cc src/BkgGen.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/P8Gen.o: src/P8Gen.cc src/P8Gen.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/JetClusterer.o: src/JetClusterer.cc src/JetClusterer.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/JetIndicesMatcher.o: src/JetIndicesMatcher.cc src/JetIndicesMatcher.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/JetTree.o: src/JetTree.cc src/JetTree.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/RhoMedianBkgClusterer.o: src/RhoMedianBkgClusterer.cc src/RhoMedianBkgClusterer.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/JetBranch.o: src/JetBranch.cc src/JetBranch.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@

obj/IP_geom_matcher.o: src/IP_geom_matcher.cc src/IP_geom_matcher.h
	${CC} ${CFLAGS} ${LIB_TRI} -c $< -o $@
