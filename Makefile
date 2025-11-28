CXX = gcc
# -O5 doesn’t exist (highest is -O3); also fix missing dash before -mavx2 -mavx512f -mavx512dq 
CXXFLAGS = -O3 -march=native -mtune=native -funroll-loops -DNDEBUG 
LDFLAGS = -fopenmp

all: optimized/sw_opt

optimized/sw_opt: optimized/sw_opt.c
	mkdir -p optimized
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f optimized/sw_opt

