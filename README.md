
<b>JaLAJni</b> is a JAVA package providing a java interface for lapack and blas library. This package contains two different versions: JaLAJni and JaLAJniLite. Both of them serve as a java interface for lapack and blas library. The main difference is that the former one - JaLAJni - calls an intermediate c interface - cblas and lapacke - which requires previous installation. So if you want to use this version called JaLAJni, you have to install cblas and lapacke first (detailed installation instructions are listed below). The latter one (JaLAJniLite), is a simplified version which requires no installation of any libraries other than lapack and blas. In JaLAJni, some machines might fail to compile the source file if you link both the libblas and libcblas when building the blas interface, while in some other ones you cannot compile without any of them. But the version - JaLAJniLite - does not involve such issues.  



Build Instructions
------------------
#### JaLAJni




#### JaLAJniLite

* JaLAJniLite requires the installation of lapack and blas. 

* To compile the package, go to the folder src/JaLAJniLite, find the Makefile and run "make" on the command line. Notice that you may have to change the extension of generated libraries in the Makefile base on your operating system. On OS X you have to change all the extensions of dynamic library to .dylib while on Linux the corresponding extensions are .so or .a. 

* To clean generated file, type “make clean” on the command line. 


Running the tests
-----------------
For testing, change directory to test/JaLAJni_tests (test/JaLAJniLite_tests), then type “make” to run the test files. If you want to clean testing results and all class files, type "make clean".  


Notes
---------
This package is intended for some basic problems we encounter when solving linear algebra problems. So we only include several most basic and widely used routines of the blas and lapack library.


Source Repository
-----------------
JaLAJni's source-code repository is hosted here on GitHub.


Authors
---------




#### JaLAJniLite
| Lu Zhang (maintainer)|    Department of Biostatistics  UCLA | lu.zhang@ucla.edu|
| LiZhen Nie | Department of Statistics Chicago University |    |                             
| Sudipto Banerjee   |  Department of Biostatistics  UCLA |     |
                             


Licensing
---------
JaLAJni is licensed under the Creative Commons Attribution License. 



