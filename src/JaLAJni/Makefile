DEST = .
CC = gcc
CFLAGS = $(INCLUDE) $(LIB) 

JAVA_HOME = ../..

JNI_SRC = $(JAVA_HOME)/src/JaLAJni

INCLUDE = -I$(JNI_SRC)/c/include \
	-I$(JAVA_HOME)/include
LIB = -L/usr/local/lib  
C_SRC_LAPACK = $(JNI_SRC)/c/jniLAPACKE.c
C_SRC_BLAS = $(JNI_SRC)/c/jniCBLAS.c
DYLIB_LAPACK = $(JAVA_HOME)/lib/libjni_java_lapacke.dylib
DYLIB_BLAS = $(JAVA_HOME)/lib/libjni_cblas.dylib
	

# For package JaLAJni
JC = javac
JH = javah
JFLAGS = -d
JAVA_SRC_LAPACK = $(JNI_SRC)/java/jniLAPACKE.java
JAVA_SRC_BLAS = $(JNI_SRC)/java/jniCBLAS.java
PACKAGE = $(JAVA_HOME)/class/JaLAJni
HEADER_LAPACK = $(JNI_SRC)/c/include/jniLAPACKE.h
HEADER_BLAS = $(JNI_SRC)/c/include/jniCBLAS.h

# header file jniLAPACKE.h


all: $(PACKAGE) $(HEADER_LAPACK) $(HEADER_BLAS) $(DYLIB_LAPACK) $(DYLIB_BLAS)

$(PACKAGE): $(JAVA_SRC_LAPACK) $(JAVA_SRC_BLAS)
	@mkdir -p $(@D)
	$(JC) $(JFLAGS) $(JAVA_HOME)/class $(JAVA_SRC_LAPACK) $(JAVA_SRC_BLAS)

$(HEADER_LAPACK): 
	@mkdir -p $(@D)
	$(JH) -o $(HEADER_LAPACK) -classpath $(JAVA_HOME)/class JaLAJni.jniLAPACKE 

$(HEADER_BLAS):
	@mkdir -p $(@D)
	$(JH) -o $(HEADER_BLAS) -classpath $(JAVA_HOME)/class JaLAJni.jniCBLAS

$(DYLIB_LAPACK): $(C_SRC_LAPACK)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK) -llapacke -llapack -lblas -o $(DYLIB_LAPACK)

$(DYLIB_BLAS): $(C_SRC_BLAS) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_BLAS) -llapack -lblas -lcblas -o $(DYLIB_BLAS)

clean: 
	rm -fr $(JAVA_HOME)/class/JaLAJni
	rm -f $(JNI_SRC)/c/include/*.h
	rm -f $(JAVA_HOME)/lib/libjni_java_lapacke.dylib
	rm -f $(JAVA_HOME)/lib/libjni_cblas.dylib
	rm -f $(JAVA_HOME)/test/JaLAJni_tests/*.class
	rm -f $(JAVA_HOME)/test/JaLAJni_tests/*.txt
