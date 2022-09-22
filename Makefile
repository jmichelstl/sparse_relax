CC=g++ --std=c++11 -O3 -fexceptions -fPIC -fopenmp
INC=-I/usr/include/suitesparse
SPARSELINK= -lsuitesparseconfig
GPULINK= -lsuitesparseconfig -lSuiteSparse_GPURuntime -lGPUQREngine

all:
	make libflists.a
	make sp_relax
	make grelax

sp_relax: libflists.a
	$(CC) $(INC) src/sparse_relax_no_gpu.cpp -o sp_relax libflists.a \
	-L/usr/local/lib -lcholmod -lspqr -lblas \
	-lboost_filesystem -lboost_regex

grelax:	libflists.a
	$(CC) -DGPU_BLAS $(INC) src/sparse_relax.cpp -o grelax libflists.a \
	-L/usr/local/lib -lcholmod -lspqr $(GPULINK) -lblas \
	-lboost_filesystem -lboost_regex

libflists.a:
	g++ --std=c++11 -c src/file_lists.cpp -o libflists.a


sclean:
	rm sp_relax

lclean:
	rm libflists.a

gclean:
	rm grelax

clean:
	rm sp_relax
	rm libflists.a
	rm grelax
