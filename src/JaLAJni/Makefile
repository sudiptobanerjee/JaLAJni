DEST = .
CC = gcc
CFLAGS = $(INCLUDE) $(LIB) 

JAVA_HOME = ../..
JAVA_PACKAGES = $(JAVA_HOME)/packages
JAVA_LIB = $(JAVA_HOME)/lib

JNI_SRC = $(JAVA_HOME)/src/JaLAJni

INCLUDE = -I$(JNI_SRC)/c/include \
	-I$(JAVA_HOME)/include
LIB = -L/usr/local/lib -L$(JAVA_LIB)
C_SRC_LAPACK = $(JNI_SRC)/c/jniLAPACKE.c
C_SRC_BLAS = $(JNI_SRC)/c/jniCBLAS.c
DYLIB_LAPACK = $(JAVA_HOME)/lib/libjni_java_lapacke.so
DYLIB_BLAS = $(JAVA_HOME)/lib/libjni_cblas.so
	

# For package Jama
JC = javac
JH = javah
JFLAGS = -d
JAVA_SRC_LAPACK = $(JNI_SRC)/java/jniLAPACKE.java
JAVA_SRC_BLAS = $(JNI_SRC)/java/jniCBLAS.java
PACKAGE = $(JAVA_PACKAGES)/JaLAJni
HEADER_LAPACK = $(JNI_SRC)/c/include/jniLAPACKE.h
HEADER_BLAS = $(JNI_SRC)/c/include/jniCBLAS.h

# header file jniLAPACKE.h


all: $(PACKAGE) $(HEADER_LAPACK) $(HEADER_BLAS) $(DYLIB_LAPACK) $(DYLIB_BLAS)

$(PACKAGE): 
	$(JC) $(JFLAGS) $(JAVA_PACKAGES) $(JAVA_SRC_LAPACK) $(JAVA_SRC_BLAS)

$(HEADER_LAPACK): 
	$(JH) -o $(HEADER_LAPACK) -classpath $(JAVA_PACKAGES) JaLAJni.jniLAPACKE 

$(HEADER_BLAS):
	$(JH) -o $(HEADER_BLAS) -classpath $(JAVA_PACKAGES) JaLAJni.jniCBLAS

$(DYLIB_LAPACK): $(C_SRC_LAPACK)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK) -llapacke -llapack -lblas -o $(DYLIB_LAPACK)

$(DYLIB_BLAS): $(C_SRC_BLAS) 
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_BLAS) -llapack -lblas -lcblas -o $(DYLIB_BLAS)

clean: 
	rm -fr $(JAVA_PACKAGES)/JaLAJni
	rm -f $(JNI_SRC)/c/include/*.h
	rm -f $(JAVA_HOME)/lib/libjni_java_lapacke.so
	rm -f $(JAVA_HOME)/lib/libjni_cblas.so
