# Compiler settings
CXX = g++
CXXFLAGS = -O2


# Rules
default : buffergas.o flowfieldentries.o molecule.o vec3d.o simulation \
	  load_and_convert_u_data


buffergas.o : buffergas.cpp
	$(CXX) $(CXXFLAGS) -c buffergas.cpp

flowfieldentries.o : flowfieldentries.cpp
	$(CXX) $(CXXFLAGS) -c flowfieldentries.cpp

molecule.o : molecule.cpp
	$(CXX) $(CXXFLAGS) -c molecule.cpp

vec3d.o : vec3d.cpp
	$(CXX) $(CXXFLAGS) -c vec3d.cpp

simulation: simulation.cpp buffergas.o flowfieldentries.o molecule.o
	$(CXX) $(CXXFLAGS) -o simulation simulation.cpp buffergas.o flowfieldentries.o molecule.o 

load_and_convert_u_data : load_and_convert_u_data.cpp vec3d.o flowfieldentries.o
	$(CXX) $(CXXFLAGS) -o load_and_convert_u_data load_and_convert_u_data.cpp vec3d.o flowfieldentries.o


clean :
	rm *.o simulation load_and_convert_u_data