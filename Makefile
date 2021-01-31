CXX=mpiicpc
CXXFLAGS=-axMIC-AVX512 -qopenmp
OPTRPT=-qopt-report=5

default : app

L.o : L.cc L.h
	${CXX} -c ${OPTRPT} ${CXXFLAGS} -o "$@" "$<"

worker.o : worker.cc L.o
	${CXX} -c ${OPTRPT} ${CXXFLAGS} -o "$@" "$<" L.o

app : main.cc worker.o
	${CXX} ${OPTRPT} ${CXXFLAGS} -o "$@" "$<" worker.o L.o

now:
	mpirun -np 4 ./app 0.2 6 14

clean :
	rm app worker.o L.o *.optrpt
