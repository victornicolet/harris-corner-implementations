APP=harris
IMAGE=../images/grand_canyon2.jpg
NRUNS=11
include ./common.mk
SOURCE_SUFFIXES = opt naive dyn nopar larger

$(APP): $(APP)_opt.so
	python main.py opt $(IMAGE) $(NRUNS)

$(APP)_naive: $(APP)_naive.so
	python main.py naive $(IMAGE) $(NRUNS)

$(APP)_dyntile: $(APP)_dyn.so
	python main.py dyn $(IMAGE) $(NRUNS)

$(APP)_nopar: $(APP)_nopar.so
	python main.py nopar $(IMAGE) $(NRUNS)

$(APP)_larger: $(APP)_larger.so
	python main.py larger $(IMAGE) $(NRUNS)

$(APP)_noverlap: $(APP)_noverlap.so
	python main.py  noverlap $(IMAGE) $(NRUNS)
