CPPFLAGS=-g -I$(ROOTSYS)/include -std=c++14 -I $(RATROOT)/include
LDFLAGS=-g -L$(ROOTSYS)/lib -L $(RATROOT)/lib
LDLIBS=-lRATEvent `root-config --libs --cflags` -lcurl -pthread

all: timing_ntuple timing_root timing_root_digitizer

timing_ntuple: timing_ntuple.cc
	$(CXX) $(LDFLAGS) -o timing_ntuple timing_ntuple.cc $(CPPFLAGS) $(LDLIBS)

timing_root: timing_root.cc
	$(CXX) $(LDFLAGS) -o timing_root timing_root.cc $(CPPFLAGS) $(LDLIBS)

timing_root_digitizer: timing_root_digitizer.cc
	$(CXX) $(LDFLAGS) -o timing_root_digitizer timing_root_digitizer.cc $(CPPFLAGS) $(LDLIBS)

timing_ntuple.o: timing_ntuple.cc

timing_root.o: timing_root.cc

timing_root.o: timing_root_digitizer.cc

clean:
	$(RM) timing_root timing_ntuple timing_root_digitizer
