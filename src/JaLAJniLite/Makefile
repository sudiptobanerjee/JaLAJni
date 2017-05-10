DEST = .
CC = gcc
CFLAGS = $(INCLUDE) $(LIB) 

JAVA_HOME = ../..
JAVA_PACKAGES = $(JAVA_HOME)/packages
JAVA_LIB = $(JAVA_HOME)/lib

JNI_SRC = $(JAVA_HOME)/src/JaLAJniLite

INCLUDE = -I$(JNI_SRC)/c/include \
	-I$(JAVA_HOME)/include
LIB = -L/usr/local/lib -L$(JAVA_LIB)
C_SRC_LAPACK = $(JNI_SRC)/c/jniLAPACK.c
C_SRC_BLAS = $(JNI_SRC)/c/jniBLAS.c
DYLIB_LAPACK = $(JAVA_HOME)/lib/liblapack_lite.so
DYLIB_BLAS = $(JAVA_HOME)/lib/libblas_lite.so
	

# For package Jama
JC = javac
JH = javah
JFLAGS = -d
JAVA_SRC_LAPACK = $(JNI_SRC)/java/jniLAPACK.java
JAVA_SRC_BLAS = $(JNI_SRC)/java/jniBLAS.java
PACKAGE = $(JAVA_PACKAGES)/JaLAJniLite
HEADER_LAPACK = $(JNI_SRC)/c/include/jniLAPACK.h
HEADER_BLAS = $(JNI_SRC)/c/include/jniBLAS.h

# header file jniLAPACKE.h


all: $(PACKAGE) $(HEADER_LAPACK) $(HEADER_BLAS) $(DYLIB_LAPACK) $(DYLIB_BLAS)

$(PACKAGE): 
	$(JC) $(JFLAGS) $(JAVA_PACKAGES) $(JAVA_SRC_LAPACK) $(JAVA_SRC_BLAS)

$(HEADER_LAPACK): 
	$(JH) -o $(HEADER_LAPACK) -classpath $(JAVA_PACKAGES) JaLAJniLite.jniLAPACK

$(HEADER_BLAS):
	$(JH) -o $(HEADER_BLAS) -classpath $(JAVA_PACKAGES) JaLAJniLite.jniBLAS

$(DYLIB_LAPACK): $(C_SRC_LAPACK)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK) -llapack -o $(DYLIB_LAPACK)

$(DYLIB_BLAS): $(C_SRC_BLAS) 
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_BLAS) -lblas -o $(DYLIB_BLAS)

clean: 
	rm -fr $(JAVA_PACKAGES)/JaLAJniLite
	rm -f $(JNI_SRC)/c/include/*.h
	rm -f $(JAVA_HOME)/lib/liblapack_lite.so
	rm -f $(JAVA_HOME)/lib/libblas_lite.so
