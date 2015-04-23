CXX=icpc
CXX_FLAGS=-openmp -ipo -g -O3 -xhost -fPIC -shared -debug parallel

#CXX=g++
#CXX_FLAGS=-fopenmp -g -O3 -fPIC -shared

SRCS=bench_source/$(wildcard *.cpp)

NO_OVERLAP_SUFFIX=_no_overlap
DYN_SUFFIX=_dyn
OPT_SUFFIX=_opt
LIB_PREFIX=lib

_CXX_FLAGS =-O -Lbench_source/ -fPIC -lharris -Wall -Werror -g -openmp
RUNPATH = -Wl,-rpath=bench_source/
OPCV_FLAGS =`pkg-config --cflags opencv`
LDFLAGS+=`pkg-config --libs opencv`


APP=harris

.PHONY : all

all : mesg mtest harris_libs

mesg :
	@echo "Parameters : "
	@echo "- Dyn suffix : " $(DYN_SUFFIX)
	@echo "- NO_OVERLAP suffix : " $(NO_OVERLAP_SUFFIX)
	@echo "- Polymage optimized : " $(OPT_SUFFIX)
	@echo $(MESG)
	@echo "__________________________________________"

help:
	@echo "Run tests :"
	@echo "\t make IMAGE=.. RUNS=.. run_test"


run_test:
	./mtest $(IMAGE) $(RUNS)
# Test file
mtest : main.cpp
	rm -f ./mtest
	$(CXX) $(_CXX_FLAGS) $(RUNPATH) $(OPCV_FLAGS) $(LDFLAGS) $< -o $@

$(APP)_libs : $(APP)$(DYN_SUFFIX) $(APP)$(NO_OVERLAP_SUFFIX) $(APP)$(OPT_SUFFIX)

$(APP)$(NO_OVERLAP_SUFFIX): $(LIB_PREFIX)$(APP)$(NO_OVERLAP_SUFFIX).so

$(APP)$(DYN_SUFFIX): $(LIB_PREFIX)$(APP)$(DYN_SUFFIX).so

$(APP)$(OPT_SUFFIX): $(LIB_PREFIX)$(APP)$(OPT_SUFFIX).so

$(LIB_PREFIX)$(APP)$(NO_OVERLAP_SUFFIX).so: bench_source/$(APP)_larger_noverlap.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@

$(LIB_PREFIX)$(APP)$(DYN_SUFFIX).so: bench_source/$(APP)_dyntile.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@

$(LIB_PREFIX)$(APP)$(OPT_SUFFIX).so: bench_source/$(APP)_polymage.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@

clean:
	rm -f ./mtest bench_source/*.pyc bench_source/*.so bench_source/graph.png
