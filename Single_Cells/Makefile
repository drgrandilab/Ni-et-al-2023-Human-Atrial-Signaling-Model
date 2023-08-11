# Makefile


CXX=icpc

# CXX=mpiicpc

INC =  -Ilib -#I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/
CFLAGS =  -Ilib  -std=c++11 -g
CXXFLAGS = -Ilib -std=c++11 -O3 -lsundials_cvode  -xcore-avx2  -axcore-avx2 -g #-lopenblas # -axcommon-avx512 -g #-lm


LDFLAGS=  -lz -lsundials_cvode 
LDLIBS=


LDFLAGS=  -lz 
LDLIBS=-lz  



BUILD_DIR=lib/build

INC_PARAMS=-Ilib -Ilib_tissue -I/usr/local/include

# CXX = g++
srcs = $(wildcard lib_tissue/*.cpp lib/*.cpp)
objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps = $(srcs:.cpp=.d)


main_HAM_Signalling_cvode_new: main_HAM_Signalling_cvode_new.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@

Atria_2D_MPI_Ghost: Atria_2D_MPI_Ghost.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@



$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(INC_PARAMS)   $(CXXFLAGS)  -MMD -MP -c $< -o $@ 
.PHONY: clean

# $(RM) is rm -f by default
clean:
	$(RM) $(objs) $(deps) main_HAM_Signalling_population main_HAM_Signalling_population_level main_HAM_Signalling_cvode_new

-include $(deps)

MKDIR_P ?= mkdir -p