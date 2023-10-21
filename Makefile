EXECUTABLE   = par-heat-tranfer-2d.x
CXX          = mpic++
CXXFLAGS     = -std=c++11 -O3

OBJS = src/*.cpp  HeatTransfer2D.cpp
$(EXECUTABLE): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ 

.PHONY: clean
clean:
	$(RM) $(EXECUTABLE)
