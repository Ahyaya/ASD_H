# Makefile for test purpose
#
# The project ASD_H is greatly embodied from FINUFFT
# Author email: c.houyuan@mail.scut.edu.cn

CXX = g++
CC = gcc
FC = gfortran
CLINK = -lstdc++
FLINK = $(CLINK)
CFLAGS := -O3 -funroll-loops -march=native -fcx-limited-range
CXXFLAGS := $(CFLAGS)
CXXFLAGS := $(CFLAGS)

FFTWNAME = fftw3
# linux default is fftw3_omp, since 10% faster than fftw3_threads...
FFTWOMPSUFFIX = omp
LIBS := -lm
# multithreading for GCC: C++/C/Fortran, MATLAB, and octave (ICC differs)...
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

# absolute path of this makefile, ie FINUFFT's top-level directory...
FINUFFT = $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

# for macos compatible
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
	CXX=clang++
	CC=clang
	CFLAGS = -O3
	CFLAGS += -I include -I/usr/local/include -I/usr/local/opt/libomp/include -I/opt/homebrew/include
	FFLAGS   = $(CFLAGS)
	CXXFLAGS = $(CFLAGS)
	LIBS += -L/usr/local/lib -L/opt/homebrew/lib
	# OpenMP with clang needs following...
	OMPFLAGS = -Xpreprocessor -fopenmp
	OMPLIBS = -L/usr/local/lib -L/usr/local/opt/libomp/lib -lomp
	# since fftw3_omp doesn't work in OSX, we need...
	FFTWOMPSUFFIX=threads
endif

INCL = -Iinclude
CXXFLAGS := $(CXXFLAGS) $(INCL) -fPIC -std=c++14
CFLAGS := $(CFLAGS) $(INCL) -fPIC
FFLAGS := $(FFLAGS) $(INCL) -I/usr/include -fPIC

LIBSFFT := -l$(FFTWNAME) -l$(FFTWNAME)f $(LIBS)

# multi-threaded libs & flags, and req'd flags (OO for new interface)...
ifneq ($(OMP),OFF)
  CXXFLAGS += $(OMPFLAGS)
  CFLAGS += $(OMPFLAGS)
  FFLAGS += $(OMPFLAGS)
  LIBS += $(OMPLIBS)
  ifneq ($(MINGW),ON)
    ifneq ($(MSYS),ON)
# omp override for total list of math and FFTW libs (now both precisions)...
      LIBSFFT := -l$(FFTWNAME) -l$(FFTWNAME)_$(FFTWOMPSUFFIX) -l$(FFTWNAME)f -l$(FFTWNAME)f_$(FFTWOMPSUFFIX) $(LIBS)
    endif
  endif
endif

# name & location of library we're building...
LIBNAME = libfinufft
DYNLIB = lib/$(LIBNAME).so
DYNLIB_ASDH = lib/libASD_H.so

STATICLIB = lib/$(LIBNAME).a
# absolute path to the .so, useful for linking so executables portable...
ABSDYNLIB = $(FINUFFT)$(DYNLIB)
ABSDYNLIB_ASDH = $(FINUFFT)$(DYNLIB_ASDH)

# spreader is subset of the library with self-contained testing, hence own objs:
# double-prec spreader object files that also need single precision...
SOBJS = src/spreadinterp.o src/utils.o
# their single-prec versions
SOBJSF = $(SOBJS:%.o=%_32.o)
# precision-dependent spreader object files (compiled & linked only once)...
SOBJS_PI = src/utils_precindep.o
# spreader dual-precision objs
SOBJSD = $(SOBJS) $(SOBJSF) $(SOBJS_PI)

# double-prec library object files that also need single precision...
OBJS = $(SOBJS) src/finufft.o src/simpleinterfaces.o
# their single-prec versions
OBJSF = $(OBJS:%.o=%_32.o)
# precision-dependent library object files (compiled & linked only once)...
OBJS_PI = $(SOBJS_PI) contrib/legendre_rule_fast.o
# all lib dual-precision objs
OBJSD = $(OBJS) $(OBJSF) $(OBJS_PI)

# collect headers for implicit depends
HEADERS = $(wildcard include/*.h)

# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $< -o $@
%_32.o: %.cpp $(HEADERS)
	$(CXX) -DSINGLE -c $(CXXFLAGS) $< -o $@
%.o: %.c $(HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@
%_32.o: %.c $(HEADERS)
	$(CC) -DSINGLE -c $(CFLAGS) $< -o $@

# included auto-generated code dependency...
src/spreadinterp.o: src/ker_horner_allw_loop.c src/ker_lowupsampfac_horner_allw_loop.c

default: $(DYNLIB_ASDH)

# build share libs
lib: $(STATICLIB) $(DYNLIB) $(DYNLIB_ASDH)
	@echo "lib: native build success!"

$(STATICLIB): $(OBJSD)
	ar rcs $(STATICLIB) $(OBJSD)
ifeq ($(OMP),OFF)
	@echo "$(STATICLIB) built, single-thread version"
else
	@echo "$(STATICLIB) built, multithreaded version"
endif
$(DYNLIB): $(OBJSD)
# using *absolute* path in the -o here is needed to make portable executables
# when compiled against it, in mac OSX, strangely...
	$(CXX) -shared $(OMPFLAGS) $(OBJSD) -o $(DYNLIB) $(LIBSFFT)
ifeq ($(OMP),OFF)
	@echo "$(DYNLIB) built, single-thread version"
else
	@echo "$(DYNLIB) built, multithreaded version"
endif

$(DYNLIB_ASDH): $(DYNLIB) src/ASD_H.c
	$(CC) -shared -fPIC src/ASD_H.c $(DYNLIB) -lm -o $(DYNLIB_ASDH)
	 
# here $(OMPFLAGS) and $(LIBSFFT) is even needed for linking under mac osx.
# see: http://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html
# Also note -l libs come after objects, as per modern GCC requirement.

test: main/demo.out main/exampleSignals/X.dat main/lib
	@echo "test start"
	cd main && ./demo.out && cd ../
	@echo "test complete"

main/demo.out: $(DYNLIB) $(DYNLIB_ASDH) main/demo.c
	$(CC) main/demo.c $(ABSDYNLIB) $(ABSDYNLIB_ASDH) -lm -o main/demo.out

main/exampleSignals/X.dat: main/signalGen.out
	@echo "Generating example signal"
	cd main && ./signalGen.out && cd ../

main/signalGen.out: main/signalGen.c
	$(CC) main/signalGen.c -lm -o main/signalGen.out

main/lib: $(DYNLIB_ASDH)
	cd main && ln -s $(FINUFFT)lib ./ && cd ../

clean:
	rm -f main/*.test
	rm -f main/*_flt.csv

objclean:
	rm -f src/*.o
	rm -f contrib/*.o

cleanall: clean objclean
	rm -f main/*.out

web: $(DYNLIB) $(DYNLIB_ASDH) main/main.c main/lib
	$(CC) main/main.c $(ABSDYNLIB) $(ABSDYNLIB_ASDH) -lpthread -o main/fftwebcall.out
