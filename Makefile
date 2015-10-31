ReadData_C:
	h5c++ -I$(ROOTSYS)/include -o calc_charge charge.cpp -L$(ROOTSYS)/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -l\
Rint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lhdf5_cpp -g3 -std=c++0x -Wall
