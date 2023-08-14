# Makefile
# CXX=icpc
# CFLAGS =  -Ilib  -std=c++11
# CXXFLAGS = -Ilib  -std=c++11 -g -O0 -lm# -Wall# -O0 -no-ip -w3 -ip-no-inlining#-no-vec  #-qopt-report=5# -static -I/opt/OpenBLAS/include  -L/opt/OpenBLAS/lib  -lz -lopenblas #-qopt-report=5
# CXX=icpc

CXX=mpiicpc

INC =  -Ilib -#I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/

CXXFLAGS = -Ilib -std=c++11 -O3  -qopenmp  -lsundials_cvode  -g -lz #-xcore-avx2  -axcore-avx2   #-lopenblas # -axcommon-avx512 -g #-lm
# if wanna accuracy: add to CXXFLAGS -fp-model consistent


LDFLAGS=  -lz -lsundials_cvode #-lopenblas #-lsundials_sunlinsolklu   # -lm# -lopenblas

LDFLAGS=  -lz # -lopenblas #-lm


BUILD_DIR=lib/build

INC_PARAMS=-Ilib -Ilib_tissue -I/usr/local/include -I/home/hbni/sundials/include/ -L/home/hbni/sundials/lib/# (if g++)

# CXX = g++
srcs = $(wildcard lib_tissue/*.cpp lib/*.cpp)
objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps = $(srcs:.cpp=.d)

Atria_2D_MPI_Ghost: Atria_2D_MPI_Ghost.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@

ONE_D_MPI_Ghost: ONE_D_MPI_Ghost.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@
ONE_D_MPI_Ghost_VW: ONE_D_MPI_Ghost_VW.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@


$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(INC_PARAMS)   $(CXXFLAGS)  -MMD -MP -c $< -o $@ 
.PHONY: clean

# $(RM) is rm -f by default
clean:
	$(RM) $(objs) $(deps) main_HAM_Signalling_population main_HAM_Signalling_population_level

-include $(deps)

MKDIR_P ?= mkdir -p
