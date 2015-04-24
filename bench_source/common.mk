CXX=icpc
CXX_FLAGS=-openmp -ipo -O3 -xhost -fPIC -shared

#CXX=g++
#CXX_FLAGS=-fopenmp -g -O3 -fPIC -shared

all: polymage naive dyntile sequential larger

polymage: $(APP)_opt.so
naive: $(APP)_naive.so
dyntile: $(APP)_dyn.so
sequential: $(APP)_sequential.so
larger: $(APP)_larger
opt_r:$(APP)_opt_r.so

$(APP)_opt.so: $(APP)_polymage.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_opt_r.so: $(APP)_polymage_rewritten.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_naive.so: $(APP)_polymage_naive.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_dyn.so: $(APP)_dyntile.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_sequential.so: $(APP)_sequential.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@

$(APP)_larger.so: $(APP)_larger_noverlap.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@



clean:
	rm -f *.pyc *.so
