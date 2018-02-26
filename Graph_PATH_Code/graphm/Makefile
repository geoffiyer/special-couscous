CXX ?= g++
CFLAGS = -Wno-deprecated -O3 -fPIC
#LDFLAGS= -lgsl -lgslcblas
SHVER = 1
all: ./bin/graphm
./bin/graphm: ./obj/algorithm_ca.o ./obj/algorithm_path.o ./obj/algorithm_ext.o ./obj/algorithm_fsol.o  ./obj/algorithm_iden.o ./obj/algorithm_lp.o ./obj/algorithm.o ./obj/algorithm_path.o ./obj/algorithm_qcv.o ./obj/algorithm_rand.o ./obj/algorithm_rank.o ./obj/algorithm_umeyama.o ./obj/algorithm_unif.o ./obj/experiment.o ./obj/graph.o ./obj/hungarian.o  ./obj/rpc.o main.cpp
	$(CXX) $(CFLAGS)  main.cpp ./obj/*o $(LDFLAGS) -o ./bin/graphm
#interface classes
./obj/experiment.o:   experiment.cpp experiment.h
		$(CXX) $(CFLAGS) -c experiment.cpp  -o ./obj/experiment.o
./obj/rpc.o:   rpc.cpp rpc.h
		$(CXX) $(CFLAGS) -c rpc.cpp  -o ./obj/rpc.o
./obj/graph.o:   graph.cpp graph.h
		$(CXX) $(CFLAGS) -c graph.cpp  -o ./obj/graph.o
./obj/hungarian.o:   hungarian.cpp hungarian.h
		$(CXX) $(CFLAGS) -c hungarian.cpp  -o ./obj/hungarian.o
./obj/algorithm.o:   algorithm.cpp algorithm.h
		$(CXX) $(CFLAGS) -c algorithm.cpp  -o ./obj/algorithm.o
#graph matching algorithms
./obj/algorithm_ca.o: algorithm_ca.cpp algorithm_ca.h
		$(CXX) $(CFLAGS) -c algorithm_ca.cpp -o ./obj/algorithm_ca.o
./obj/algorithm_path.o: algorithm_path.cpp algorithm_path.h
		$(CXX) $(CFLAGS) -c algorithm_path.cpp -o ./obj/algorithm_path.o
./obj/algorithm_ext.o:   algorithm_ext.cpp algorithm_ext.h
		$(CXX) $(CFLAGS) -c algorithm_ext.cpp  -o ./obj/algorithm_ext.o
./obj/algorithm_fsol.o:   algorithm_fsol.cpp algorithm_fsol.h
		$(CXX) $(CFLAGS) -c algorithm_fsol.cpp  -o ./obj/algorithm_fsol.o
./obj/algorithm_iden.o:   algorithm_iden.cpp algorithm_iden.h
		$(CXX) $(CFLAGS) -c algorithm_iden.cpp  -o ./obj/algorithm_iden.o
./obj/algorithm_lp.o:   algorithm_lp.cpp algorithm_lp.h
		$(CXX) $(CFLAGS) -c algorithm_lp.cpp  -o ./obj/algorithm_lp.o
./obj/algorithm_qcv.o:   algorithm_qcv.cpp algorithm_qcv.h
		$(CXX) $(CFLAGS) -c algorithm_qcv.cpp  -o ./obj/algorithm_qcv.o
./obj/algorithm_rand.o:   algorithm_rand.cpp algorithm_rand.h
		$(CXX) $(CFLAGS) -c algorithm_rand.cpp  -o ./obj/algorithm_rand.o
./obj/algorithm_rank.o:   algorithm_rank.cpp algorithm_rank.h
		$(CXX) $(CFLAGS) -c algorithm_rank.cpp  -o ./obj/algorithm_rank.o
./obj/algorithm_umeyama.o:   algorithm_umeyama.cpp algorithm_umeyama.h
		$(CXX) $(CFLAGS) -c algorithm_umeyama.cpp  -o ./obj/algorithm_umeyama.o
./obj/algorithm_unif.o:   algorithm_unif.cpp algorithm_unif.h
		$(CXX) $(CFLAGS) -c algorithm_unif.cpp  -o ./obj/algorithm_unif.o

clean:
	rm -f ./obj/*o
