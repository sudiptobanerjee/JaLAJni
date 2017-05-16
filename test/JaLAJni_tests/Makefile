LAPACK_CLASS = lapacketest.class
LAPACK_RUN = lapack_test_result.txt
BLAS_CLASS = cblastest.class
BLAS_RUN = blas_result_result.txt

all: $(LAPACK_CLASS) $(LAPACK_RUN) $(BLAS_CLASS) $(BLAS_RUN)

$(LAPACK_CLASS):
	javac -classpath ../../class/:. lapacketest.java 

$(LAPACK_RUN):
	java -Djava.library.path=../../lib/. -classpath ../../class/:. lapacketest > lapack_test_result.txt

$(BLAS_CLASS):
	javac -classpath ../../class/:. cblastest.java 

$(BLAS_RUN):
	java -Djava.library.path=../../lib/. -classpath ../../class/:. cblastest > blas_test_result.txt

clean: 
	rm -f *.txt *.class
