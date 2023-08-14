# Makefile
# CXX=icpc
# CFLAGS =  -Ilib  -std=c++11
# CXXFLAGS = -Ilib  -std=c++11 -g -O0 -lm# -Wall# -O0 -no-ip -w3 -ip-no-inlining#-no-vec  #-qopt-report=5# -static -I/opt/OpenBLAS/include  -L/opt/OpenBLAS/lib  -lz -lopenblas #-qopt-report=5
CXX=icpc

# CXX=mpiicpc

INC =  -Ilib -#I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/
CFLAGS =  -Ilib  -std=c++11 -g
CXXFLAGS = -Ilib -std=c++11 -O3 -lsundials_cvode  -xcore-avx2  -axcore-avx2 -g #-lopenblas # -axcommon-avx512 -g #-lm


# CFLAGS =  -Ilib  -std=c++11 #-I/usr/local/suitesparse/include/ -L/usr/local/suitesparse/lib/
# CXXFLAGS = -Ilib -std=c++11 -g -O3 -ffast-math -lsundials_cvode -I/usr/local/suitesparse/include/ #-rpath=/usr/local/suitesparse/lib/ #-xHost#  -axcommon-avx512

# LDFLAGS=  -lz -lsundials_cvode -lsundials_sunlinsolklu   # -lm# -lopenblas

# CXXFLAGS = -Ilib -std=c++11 -g -O3 -lsundials_cvode  -xcore-avx2  -axcore-avx2#-I/usr/local/suitesparse/include/ #-rpath=/usr/local/suitesparse/lib/ #-xHost#  -axcommon-avx512

LDFLAGS=  -lz -lsundials_cvode #-lopenblas #-lsundials_sunlinsolklu   # -lm# -lopenblas
LDLIBS=#-lz -lsundials_cvode -lsundials_sunlinsolklu -lrt -lklu -lamd -lcolamd -lbtf -lsuitesparseconfig -lsundials_sunlinsollapackdense -lsundials_sunmatrixsparse # -lm #-lopenblas -lpthread
# main_CNZ:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))

LDFLAGS=  -lz # -lopenblas #-lm
LDLIBS=-lz  #-lopenblas -lpthread #-lm
# main_CNZ:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))
# main_HAM_Signalling_population:$(patsubst  %.cpp,%.o,$(wildcard lib/*.cpp $LDFLAGS)) 

BUILD_DIR=lib/build

INC_PARAMS=-Ilib -Ilib_tissue -I/usr/local/include

# CXX = g++
srcs = $(wildcard lib_tissue/*.cpp lib/*.cpp)
objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps = $(srcs:.cpp=.d)

# main_HAM_Signalling_population: main_HAM_Signalling_population.cpp $(objs)  
# main_HAM_Signalling_population_comp: main_HAM_Signalling_population_comp.cpp $(objs)  
# main_HAM_Signalling_population_level: main_HAM_Signalling_population_level.cpp $(objs)  
# 	$(CXX) $(CXXFLAGS) $^ -o $@
# main_HAM_Signalling_population: main_HAM_Signalling_population.cpp $(objs)  
# 	$(CXX) $(CXXFLAGS) $^ -o $@
# main_HAM_Signalling_population: main_HAM_Signalling_population.cpp $(objs)  
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@
# main_HAM_Signalling_population_comp: main_HAM_Signalling_population_comp.cpp $(objs)  
# 	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@


main_HAM_Signalling_cvode_new: main_HAM_Signalling_cvode_new.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@


main_HAM_Signalling_cvode_phos_level: main_HAM_Signalling_cvode_phos_level.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@


main_HAM_Signalling_cvode_phos_comp: main_HAM_Signalling_cvode_phos_comp.cpp $(objs) 
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@

main_HAM_Signalling_cvode_new_Vm_para: main_HAM_Signalling_cvode_new_Vm_para.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@

main_HAM_Signalling_cvode_new_VC: main_HAM_Signalling_cvode_new_VC.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@

main_HAM_Signalling_cvode: main_HAM_Signalling_cvode.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@

Atria_2D_MPI_Ghost: Atria_2D_MPI_Ghost.cpp $(objs)  
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