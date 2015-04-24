CXX=icpc
CXX_FLAGS=-openmp -ipo -g -O0 -xhost -fPIC -shared -debug parallel

#CXX=g++
#CXX_FLAGS=-fopenmp -g -O3 -fPIC -shared

SRCS=bench_source/$(wildcard *.cpp)

NO_OVERLAP=_no_overlap
DYN=_dyn
OPT=_opt
NAIVE=_naive
SEQ = _seq
OPT_R=_opt_r
LIB_PREFIX=lib

_CXX_FLAGS =-O0 -Lbench_source/ -fPIC -lharris -Wall -Werror -g -openmp
RUNPATH = -Wl,-rpath=bench_source/
OPCV_FLAGS =`pkg-config --cflags opencv`
LDFLAGS+=`pkg-config --libs opencv`


APP=harris

.PHONY : all

all : mesg mtest harris_libs

mesg :
	@echo "Parameters : "
	@echo "- Dyn suffix : " $(DYN)
	@echo "- NO_OVERLAP suffix : " $(NO_OVERLAP)
	@echo "- Polymage optimized : " $(OPT)
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

$(APP)_libs : $(APP)$(DYN) $(APP)$(NO_OVERLAP) $(APP)$(OPT)

$(APP)$(NO_OVERLAP): $(LIB_PREFIX)$(APP)$(NO_OVERLAP).so

$(APP)$(DYN): $(LIB_PREFIX)$(APP)$(DYN).so

$(APP)$(NAIVE): $(LIB_PREFIX)$(APP)$(NAIVE).so
$(APP)$(SEQ): $(LIB_PREFIX)$(APP)$(SEQ).so

$(APP)$(OPT): $(LIB_PREFIX)$(APP)$(OPT).so
$(APP)$(OPT_R): $(LIB_PREFIX)$(APP)$(OPT_R).so

$(LIB_PREFIX)$(APP)$(NO_OVERLAP).so: bench_source/$(APP)_larger_noverlap.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@

$(LIB_PREFIX)$(APP)$(DYN).so: bench_source/$(APP)_dyntile.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@

$(LIB_PREFIX)$(APP)$(NAIVE).so: bench_source/$(APP)_polymage_naive.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@
$(LIB_PREFIX)$(APP)$(SEQ).so: bench_source/$(APP)_sequential.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@

$(LIB_PREFIX)$(APP)$(OPT).so: bench_source/$(APP)_polymage.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@
$(LIB_PREFIX)$(APP)$(OPT_R).so: bench_source/$(APP)_polymage_rewritten.cpp
	$(CXX) $(CXX_FLAGS) $< -o bench_source/$@

clean:
	rm -f ./mtest bench_source/*.pyc bench_source/*.so bench_source/graph.png
