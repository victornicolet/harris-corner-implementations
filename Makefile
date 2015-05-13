#Set CXX = icpc or CXX=gcc/g++

ifeq ($(CXX),icpc)
	CXX_FLAGS=-openmp -ipo -g -O3 -xhost -fPIC -shared -debug parallel
else
	CXX_FLAGS=-fopenmp -g -O3 -fPIC -shared
endif

RM=rm -f

SRCS=bench_source/$(wildcard *.cpp)

NO_OVERLAP=_no_overlap
DYN=_dyn
OPT=_opt
NAIVE=_naive
SEQ = _seq
OPT_R=_opt_r
ALIGNED=_aligned
LIB_PREFIX=lib


# Main test file compilation flags
_CXX_FLAGS =-O3 -Lbench_source/ -fPIC -Wall -g -openmp

ifeq ($(CXX),icpc)
	_CXX_FLAGS +=-openmp
else
	_CXX_FLAGS += -fopenmp
endif
ifdef CHECK_FINAL_RESULT
	_CXX_FLAGS += -D CHECK_FINAL_RESULT
endif
ifdef CHECK_LOADING_DATA
	_CXX_FLAGS += -D CHECK_LOADING_DATA
endif
ifeq ($(T),aligned)
	_CXX_FLAGS+= -D VERSION_ALIGNED
endif
RUNPATH = -Wl,-rpath=bench_source/
OPCV_FLAGS =`pkg-config --cflags opencv`
LDFLAGS = `pkg-config --libs opencv` -lharris

#Default test values
IMAGE=images/grand_canyon2.jpg
RUNS=10
APP=harris

#Listing all valid implementations
#HARRIS_IMPLEMS = harris_dyntile.so
ifeq ($(T),aligned)
	W_mesg = "\n **** When compiling test file, ensure you use the aligned \
		version with T=aligned ! **** \n"
 HARRIS_IMPLEMS = harris_n_overlap_aligned.so
 CXX_FLAGS+= -D VERSION_ALIGNED
else
	ifdef T
		HARRIS_IMPLEMS = harris_$(T).so
	else
	W_mesg = "Compiled libs ! \n If you want to use aligned version, use \
		 make T=aligned libs"
	HARRIS_IMPLEMS = harris_larger_noverlap.so
	HARRIS_IMPLEMS += harris_polymage_naive.so
	HARRIS_IMPLEMS += harris_polymage_rewritten.so
	HARRIS_IMPLEMS += harris_sequential.so
	endif
endif


.PHONY : all

all : mesg mtest harris_libs

mesg :
	@echo "Parameters : "
	@echo "__________________________________________"


run: mtest
	./mtest $(IMAGE) $(RUNS)
# Test file
mtest : main.cpp
	rm -f ./mtest
	ln -f -s libs/harris_$(T).so bench_source/libharris.so
	$(CXX) $(_CXX_FLAGS) $(RUNPATH) $(OPCV_FLAGS) $(LDFLAGS) $< -o $@

libs : $(HARRIS_IMPLEMS)
	@echo $(W_mesg)

$(HARRIS_IMPLEMS): %.so : bench_source/%.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/libs/$@ 

clean:
	$(RM) ./mtest bench_source/*.pyc bench_source/*.so bench_source/graph.png

tar:
	tar -cf harris-corner-implementations.tar bench_source/ images/ backups/ \
		Makefile main.cpp
