
.SUFFIXES: .cpp .o	
CXX=g++
CPPFLAGS=-O3
INCLUDE=-I/usr/include/python2.4
OBJS=main.o inputeqdsk.o memalloc.o meshgenerator.o postprocess.o separatrix.o vega.o

all: ${OBJS}
		${CXX} ${CPPFLAGS} -o VEGA2.0 ${OBJS}
		cp VEGA2.0 ../run/

swig: vega.i
		swig -python -c++ vega.i
		${CXX} ${CPPFLAGS} -fPIC -c main.cpp inputeqdsk.cpp memalloc.cpp meshgenerator.cpp postprocess.cpp separatrix.cpp vega.cpp vega_wrap.cxx ${INCLUDE}
		${CXX} ${CPPFLAGS} -shared -o _vega.so ${OBJS} vega_wrap.o
		cp -r vega.py _vega.so ../run/

clean: 
		rm -rf *.o VEGA2.0
		rm -rf *_wrap.* *.so *.py *.pyc 
