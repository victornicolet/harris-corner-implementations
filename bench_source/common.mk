#CXX=icpc
#CXX_FLAGS=-openmp -ipo -O3 -xhost -fPIC -shared

CXX=g++
CXX_FLAGS=-fopenmp -g -O3 -fPIC -shared

all: polymage naive dyntile nopar larger

polymage: $(APP)_opt.so
naive: $(APP)_naive.so
dyntile: $(APP)_dyn.so
nopar: $(APP)_nopar.so
larger: $(APP)_larger

$(APP)_opt.so: $(APP)_polymage.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_naive.so: $(APP)_polymage_naive.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_dyn.so: $(APP)_dyntile.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_nopar.so: $(APP)_nopar.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_larger.so: $(APP)_larger_noverlap.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@



clean:
	rm -f *.pyc *.so
