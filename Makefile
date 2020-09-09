CC=gcc
CXX=g++
CFLAGS = -O2
CXXFLAGS = $(CFLAGS)
LDFLAGS=
LDLIBS = -lblas -llapack \
	 -I /usr/local/Cellar/openblas/0.3.10_1/include/

EXES= test matinv mmse gemm dgemm_demo

all: $(EXES)

clean:
	rm -fr *.o
	rm $(EXES)

#%:%.c
#	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@
#%:%.cpp
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@
