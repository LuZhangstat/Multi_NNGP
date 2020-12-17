
all:
	g++ -c -fpic nn.cpp -o nn.o
	g++ -shared -o nn.so nn.o
	mkdir -p ../sim/pics
	mkdir -p ../sim/sim2/results
	mkdir -p ../sim/sim3/rep1/results/K1
	mkdir -p ../sim/sim3/rep1/results/K2
	mkdir -p ../sim/sim3/rep1/results/K3
	mkdir -p ../sim/sim3/rep1/results/K4
	mkdir -p ../sim/sim3/rep1/results/K5
	mkdir -p ../sim/sim3/rep1/results/K6
	mkdir -p ../sim/sim3/rep1/results/K7
	mkdir -p ../sim/sim3/rep1/results/K8
	mkdir -p ../sim/sim3/rep1/results/K9
	mkdir -p ../sim/sim3/rep1/results/K10
	mkdir -p ../sim/sim3/rep1/pics
	mkdir -p ../sim/sim4/results
	mkdir -p ../sim/sim_Factor/results/K3-lowK
	mkdir -p ../RDA/pics
	mkdir -p ../RDA/data/conj_data
	mkdir -p ../RDA/data/LMC_data
	mkdir -p ../RDA/results/BSLMC
	mkdir -p ../RDA/results/BSLMC_Factor
	mkdir -p ../RDA/results/conj
	mkdir -p ../RDA/results/misalign_resp	

clean:
	rm *.o *.so
