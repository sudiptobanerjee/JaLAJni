//
//  main.c
//  lapacktest.java
//
//  Created by Lu Zhang on 5/15/16.
//  Copyright © 2016 Lu Zhang. All rights reserved.
//

import JaLAJniLite.jni_lapack.*;
public final class lapacktest {
    private lapacktest() {}
    public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
        System.out.println("###   Testfile for JaLAJniLite jni_lapack   ### \n \n");
        int M=3, N=3, K=3;
        int[] info = new int[]{0};
        double[] A = new double[] {12, -51, 4, 6, 167, -68, -4, 24, -41};
        double[] cpA = new double[] {12, -51, 4, 6, 167, -68, -4, 24, -41};
        double[] AT = new double[] {12, 6, -4, -51, 167, 24, 4, -68, -41};
        double[] B = new double[] {4, -2, -6, -2, 10, 9, -6, 9, 14};
        double[] D = new double[9];
        double[] D2 = new double[9];
        int[] ipiv = new int[] {0, 0, 0};
        double[] tau = new double[]{0, 0, 0};
        int[] jpvt = new int[] {1, 2, 3};
        double[] sym = new double[]{1,4,0,4,2,5,0,5,3};
        double[] cpsym = new double[]{1,4,0,4,2,5,0,5,3};
        double[] eq = new double[]{1,2,3};
        double[] I = new double[]{1,0,0,0,1,0,0,0,1};
        
        double[] VL = new double[M * N];
        double[] VR = new double[M * N];
        double[] WR = new double[N];
        double[] WI = new double[N];
        double[] S = new double[N];
        double[] U = new double[M * M];
        double[] VT = new double[N * N];
        double[] superb = new double[N - 1];

        double[] work = new double[1000];
        int lwork = 3;
        int[] iwork = new int[8*M];
        //
        // Set Option
        //
        int matrix_layout = jniLAPACK.LAYOUT.RowMajor;
        char Trans = jniLAPACK.TRANSPOSE.NoTrans;
        char uplo = jniLAPACK.UPLO.Upper;
        char Jobvl = jniLAPACK.JOBV.Compute;
        char Jobvr = jniLAPACK.JOBV.Compute;
        char Jobu = jniLAPACK.JOB.All;
        char Jobvt = jniLAPACK.JOB.All;
        char Jobz = jniLAPACK.JOB.Overwritten;
        int itype = jniLAPACK.ITYPE.first;
        int result;
        //
        /* ---- LU ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   LU   ###");
        //
        // dgetrf
        //
        System.out.println();
        System.out.println("The LU decomposition A=P*L*U (L has unit diagnoal elements):");
        D = A.clone();
	
        int m = M;          // m is the number of rows of matrix A
        int n = N;          // n is the number of columns of matrix A
        double[] a = D;     // a is a general m-by-n matrix
        int lda = N;        // lda is the leading dimension of the array a,
                            // here since it is row-major, lda >= max(1, n)
                            // ipiv is the pivot indices
                            //  info is integer (array)
        printMatrix("Matrix A = ", matrix_layout, a, m, n);
        result = jniLAPACK.dgetrf( matrix_layout, m, n, a, lda, ipiv, info);

        System.out.println("dgetrf");
        printMatrix("The LU matrix of A = ", matrix_layout, a, m, n);
        printIntArray("  The permutation vector ipiv = ", ipiv, m);
        //
        // dgetrs
        //
        System.out.println();
        D = B.clone();
        D2 = A.clone();
        printMatrix("Matrix A = ", matrix_layout, D, M, N);
        printMatrix("Matrix B = ", matrix_layout, D2, N, N);
        m = M;              // m = the number of rows of matrix A
        n = N;              // n = the number of columns of matrix A
        a = D;              // a = a general m-by-n matrix to be factorized
        lda = N;            // lda = the leading dimension of the array a, here since it is row-major, lda >= max(1, n)
                            // ipiv = the pivot indices
                            // info = integer (array)
        result = jniLAPACK.dgetrf( matrix_layout, m, n, a, lda, ipiv, info);

        char trans = Trans; // trans specifies the form of the system of equations
                            // n = the order of matrix a
        int nrhs = N;       // nrhs = the number of columns of matrix b
                            // a is returned by dgetrf
                            // lda = the leading dimension of the array a,
                            // here since it is row_major, lda >= max(1, n)
                            // ipiv = the pivot indices from dgetrf
        double[] b = D2;    // b = the right hand side matrix, dimension(n,nrhs)
        int ldb = nrhs;     // ldb = the leading dimension of the array b, here since it is row-majored, ldb >= max(1, nrhs)
                            //info is integer (array)
        jniLAPACK.dgetrs(matrix_layout, trans, n, nrhs, a, lda, ipiv, b, ldb, info);
        System.out.println("dgetrs");
        printMatrix("The solution  matrix of AX = B by LU is:", matrix_layout, b, ldb, nrhs);
        //
        // dgetri
        //
        System.out.println();
        D = A.clone();
        printMatrix("Matrix A = ", matrix_layout, D, N, N);
        m = M;              // m = the number of rows of matrix A
        n = N;              // n = the number of columns of matrix A
        a = D;              // a = a general m-by-n matrix to be factorized
        lda = N;            // lda = the leading dimension of the array a, here since it is row-major, lda >= max(1, n)
                            //ipiv = the pivot indices
                            //info = integer (array)
        result = jniLAPACK.dgetrf( matrix_layout, m, n, a, lda, ipiv, info);
        n = M;              // n = the order of the matrix a, a is computed by dgetrf
                            //lda = the leading dimension of matrix a
                            //ipiv = pivot indices from dgetrf
                            //work is double precision array, dimension((max(1, lwork))
        lwork = N;          // lwork  = the dimension of the array work, lwork >= max(1,n)
                            //info is integer(array)
        jniLAPACK.dgetri(matrix_layout, n, a, lda, ipiv, work, lwork, info);
        System.out.println("dgetri");
        printMatrix("  The Inverse of  matrix A = ", matrix_layout, a, n, n);
        //
        /* ---- Cholesky ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   Cholesky   ###");
        //
        // dpotrf 
        //
        System.out.println();
        D = B.clone();
                            //uplo specifies which part of a is stored
        n = M;              // n = the order of matrix a
        a = D;              // a = real symmetric positive definite matrix
        lda = M;            // lda = the leading dimension of the array a
                            //info is integer(array)
        printMatrix("Matrix A = ", matrix_layout, a, n, n);
        jniLAPACK.dpotrf(matrix_layout, uplo, n, a, lda, info);
        System.out.println("dpotrf");
        printMatrix("  The Cholesky factorization of matrix A = ", matrix_layout, a, n, n);
        //
        // dpotri 
        //
        System.out.println();
                            //uplo specifies which part of a is stored
        n = M;              // n = the order of matrix a
        a = D;              // a is real symmetric positive definite matrix
        lda = M;            // lda = the leading dimension of the array a
                            //info is integer(array)
        jniLAPACK.dpotri(matrix_layout, uplo, n, a, lda, info);
        System.out.println("dpotri");
        printMatrix("The upper triangular part of Inverse of matrix A is :", matrix_layout, a, n, n);
        //
        //  dpotrs 
        //
        System.out.println();
        D = B.clone();
        D2 = A.clone();
                            //uplo specifies which part of a is stored
        n = M;              // n = the order of matrix a
        a = D;              // a = real symmetric positive definite matrix
        lda = M;            // lda = the leading dimension of the array a
                            //info is integer(array)
        jniLAPACK.dpotrf(matrix_layout, uplo, n, a, lda, info);
                            //uplo specifies which part of a is stored
                            //n = the order of matrix a
        nrhs = N;           // nrhs = the number of columns of matrix b
                            //a is computed by dpotrf
                            //lda = the leading dimension of array a
        b = D2;             // b dimension (n, nrhs)
        ldb = nrhs;         // ldb = the leading dimension of array b
                            //info is integer(array)
        System.out.println("dpotrf");
        printMatrix("Matrix A = ", matrix_layout, a, n, n);
        printMatrix("Matrix B = ", matrix_layout, b, n, nrhs);
        jniLAPACK.dpotrs(matrix_layout, uplo, n, nrhs, a, lda, b, ldb, info);
        printMatrix("  The solution  matrix of AX=B by Cholesky is :", matrix_layout, b, ldb, nrhs);
        //
        /* ---- QR ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   QR decomposition   ###");
        //
        // dgeqrf 
        //
        System.out.println();
        System.out.println("A=QR:");
        D = A.clone();
        m = M;              // m = the number of rows of the matrix a
        n = N;              // n = the number of columns of matrix a
        a = D;              // a dimension (m, n)
        lda = n;            // lda = leading dimension of a
                            //tau = the scalar factors of the elementary
                            //reflectors, dimension (min(m,n))
                            //work dimension (max(1,lwork))
                            //lwork is integer
                            //info is integer(array)
        System.out.println("dgeqrf");
        printMatrix("Matrix A = ", matrix_layout, a, m, n);
        jniLAPACK.dgeqrf(matrix_layout, m, n, a, lda, tau, work, lwork, info);
        printMatrix("  The rewritten matrix A = ", matrix_layout, a, m, n);
        //
        // dorgqr 
        //
        System.out.println();
        m = M;              // m = the number of rows of matrix q
        n = N;              // n = the number of columns of matrix q
        int k = M;          // k = the number of elementary reflectors whose product defines the matrix q
                            // a is returned by dgeqrf
	                        // lda = leading dimension of a
	                        // tau is returned by dgeqrf
	                        // work dimension (max(1,lwork))
	                        // lwork = dimension of array work
	                        // info integer(array)
        System.out.println("dorgqr");
        jniLAPACK.dorgqr(matrix_layout, m, n, k, a, lda, tau, work, lwork, info);
        printMatrix("  The matrix Q is :", matrix_layout, a, m, n);
        //
        // dgeqp3 
        //
        System.out.println();
        D = A.clone();
        lwork = 30;
        m = M;              // m = the number of rows of matrix a
        n = N;              // n = the number of columns of matrix a
        a = D;              // a dimension(m,n)
        lda = n;            // lda = leading dimension of a
                            //jpvt integer array, dimension (n)
                            //tau double precision array, dimension (min(m,n))
                            //work double precision array, dimension (max(1,lwork))
                            //lwork integer, the dimension of array work
                            //info integer(array)
        System.out.println("dgeqp3");
        printMatrix("Matrix A = ", matrix_layout, a, m, n);
        jniLAPACK.dgeqp3(matrix_layout, m, n, a, lda, jpvt, tau, work, lwork, info);
        printMatrix("  The R matrix computed by dgeqp3:", matrix_layout, a, m, n);
        //
        /* ---- Eigenvector and SVD ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   Eigenvector and SVD   ###");
        //
        // dgeev 
        //
        System.out.println();
        System.out.println(" \nGet Eigenvectors of A:");
        D = A.clone();
        char jobvl = Jobvl; //jobvl specifies whether left eigenvectors of a are computed
        char jobvr = Jobvr; //jobvr specifies whether right eigenvectors of a are computed
        n = N;              //n = the order of matrix a
        a = D;              //a double precision array, dimension(lda,n)
        lda = n;            //n = the leading dimension of array a
        double[] wr = WR;   // wr is double precision array, dimension (n)
        double[] wi = WI;   // wi is double precision array, dimension (n)
        double[] vl = VL;   // vl is double precision array, dimension (ldvl, n)
        int ldvl = N;       // ldvl = the leading dimension of array vl
        double[] vr = VR;   // vr is double precision array, dimension (ldvr, n)
        int ldvr = N;       // ldvr = the leading dimension of array vr
                            //work double precision array, dimension(max(1,lwork))
                            //lwork = dimension of array work, integer
                            //info integer(array)
        System.out.println("dgeev");
        printMatrix("Matrix A = ", matrix_layout, a, n, n);
        jniLAPACK.dgeev(matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
        printMatrix("\n  Left eigenvectors (in columns) :", matrix_layout, vl, ldvl, n);
        printMatrix("\n  Right eigenvectors (in columns):", matrix_layout, vr, ldvr, n);
        printMatrix("\n  WR:", matrix_layout, wr, 1, n);
        printMatrix("\n  WI:", matrix_layout, wi, 1, n);
        //
        // dgesvd 
        //
        System.out.println();
        System.out.println(" \nSVD of A: A = U * SIGMA * transpose(V):");
        D = A.clone();
        lwork = 50;
        char jobu = Jobu;   // jobu specifies options for computing all or part of the matrix u
        char jobvt = Jobvt; // jobvt specifies options for computing all or part of the matrix vt
        m = M;              // m = the number of rows of the input matrix a
        n = N;              // n = the number of columns of the input matrix a
        a = D;              // a is double precision array, dimension(lda,n)
        lda = M;            // lda = the leading dimension of the array a
        double[] s = S;     // s = the singular values of a, sorted so that s(i)>=s(i+1)
        double[] u = U;     // u dimension(ldu,ucol)
        int ldu = N;        // ldu = leading dimension of array u
        double[] vt = VT;   // vt dimension(ldvt,n)
        int ldvt = N;       // ldvt = leading dimensio of array vt
                            //work double precision array, dimension(max(1,lwork))
                            //lwork = dimension of array work
                            //info integer(array)
        System.out.println("dgesvd");
        printMatrix("Matrix A = ", matrix_layout, a, m, n);
        jniLAPACK.dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
        printMatrix("\n  D = ", matrix_layout, s, 1, m);
        printMatrix("\n  U = ", matrix_layout, u, m, n);
        printMatrix("\n  VT = ", matrix_layout, vt, n, n);
        //
        // dgesdd
        //
        System.out.println();
        System.out.println(" \nSVD of A calculated by dgesdd:");
        D = A.clone();
        char jobz = Jobu;   // jobz specifies options for computing all or part of matrix u
        m = M;              // number of rows of input matrix a
        n = N;              // number of columns of input matrix a
        a = D;              // a double precision array, dimension (lda,n)
        lda = n;            // lda = leading dimension of array a
        s = S;              //s = singular values of a, sorted so that s(i)>=s(i+1)
        u = U;              //u double precision array, dimension (ldu,ucol)
        ldu = n;            //ldu = leading dimension of array u
        vt = VT;            //vt double precision array, dimension (ldvt, n)
        ldvt = n;           //ldvt = leading dimension of array vt
                            //work double precision array, dimension(max(1,lwork))
                            //lwork = dimension of the array work
                            //iwork, integer array, dimension(8*min(m,n))
                            //info integer
        System.out.println("dgesdd");
        printMatrix("Matrix A = ", matrix_layout, a, m, n);
        jniLAPACK.dgesdd(matrix_layout, jobu, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
        printMatrix("\n  D = ", matrix_layout, s, 1, m);
        printMatrix("\n  U = ", matrix_layout, u, m, n);
        printMatrix("\n  VT = ", matrix_layout, vt, n, n);
        //
        //dsyev
        //
        System.out.println();
        System.out.println
        (" \ndsyev: computes the eigenvalues and right/left eigenvectors for symmetric matrix A");
                            //jobz specifies whether eigenvectors will be computed
                            //uplo specifies which part of a is stored
        n = N;              // the order of matrix a
        a = sym;            // a is symmetric matrix
        lda = n;            // lda = leading dimension of array a
        double[] w = S;     // w = eigenvalues in ascending order, dimension (n)
                            //work double precision array, dimension(max(1,lwork))
                            //lwork = length of array work
                            //info integer(array)
        System.out.println("dsyev");
        printMatrix("matrix A = ", matrix_layout, a, n, n);
        jniLAPACK.dsyev(matrix_layout, jobvl, uplo, n, a, lda, w, work, lwork, info);
        printMatrix("  Eigenvectors of A :", matrix_layout, a, n, n);
        printMatrix("\n  The singular values of A :", matrix_layout, w, 1, n);
        //
        //dgesv
        //
        System.out.println();
        System.out.println(" \ndgesv: computes the solution to linear equation A*X = B for GE matrices");
        n = N;              // n = order of matrix a
        nrhs = 1;           // nrhs = number of columns of matrix b
        a = A;              //a = coefficient matrix a
        lda = n;            // lda = leading dimension of array a
                            //ipiv = pivot indices that define the permutation matrix p
        b = eq;             // b = rith hand side matrix B
        ldb = nrhs;         //ldb = leading dimension of array b
                            //info integer(array)
        System.out.println("dgesv");
        printMatrix("  matrix A :", matrix_layout, a, n, n);
        printMatrix("  matrix B :", matrix_layout, b, n, nrhs);
        jniLAPACK.dgesv(matrix_layout, n, nrhs, a, lda, ipiv, b, ldb, info);
        printMatrix("  X = ", matrix_layout, b, n, nrhs);
        //
        //dggev
        //
        System.out.println();
        System.out.println(" \ndggev: computes the eigenvalues and right/left eigenvectors for symmetric matrix A");
        sym = cpsym.clone();
        jobvl = Jobvl;      //jobvl specifies whether to compute the generalized eigenvectors or not
        jobvr = Jobvr;      //jobvr specifies whether to compute the generalized eigenvectors or not
        n = N;              //the order of matrix a
        a = cpA;            //matrix a in the pair (a,b)
        lda = N;            //lda = leading dimension of a
        b = sym;            //matrix b in the pair (a,b)
        ldb = N;            //ldb = leading dimension of b
        double[] alphar = S;//alphar double precision array, dimension(n)
        double[] alphai = eq;//alhpar double precision array, dimesnion(n)
        double[] beta = tau;//beta double precision array, dimension(n)
        vl = VL;            //vl = left eigenvectors
        ldvl = N;           //ldvl = leading dimension of matrix vl
        vr = VR;            //vr = right eigenvectors
        ldvr = N;           //ldvr = leading dimension of matrix vr
                            //work double precision array, dimension(max(1,lwork))
                            //lwork = dimension of array work
                            //info integer(array)
        System.out.println("dggev");
        printMatrix("  matrix A :", matrix_layout, a, n, n);
        printMatrix("  matrix B :", matrix_layout, b, n, n);
        jniLAPACK.dggev(matrix_layout, jobvl, jobvr, n, a, lda, b,
                        ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
        printMatrix("\n alphar = ", matrix_layout, alphar, 1, n);
        printMatrix("\n alphai = ", matrix_layout, alphai, 1, n);
        printMatrix("\n beta = ", matrix_layout, beta, 1, n);
        printMatrix("\n VL = ", matrix_layout, vl, ldvl, n);
        printMatrix("\n VR = ", matrix_layout, vr, ldvr, n);
        //
        //dsygv
        //
        System.out.println();
        System.out.println(" \ndsygv: computes the eigenvalues and right/left eigenvectors for real-generalized symmetric definite eigenproblem");
        sym = cpsym.clone();
                            //itype specifies the problem type to be solved
        jobz = Jobvl;       //jobvl specifies whether the eigenvectors should be computed
                            //uplo specifies which part of a and b are stored
        n = N;              //n = the order of the matrices a and b
        a = sym;            //a is symmetric matrix a
        lda = N;            //lda = leading dimension of array a
        b = I;              //b is symmetric positive definite matrix b
        ldb = n;            //ldb = leading dimension of array b
        w = S;              //w = eigenvalues in ascending order, dimension (n)
                            //work double precision array, dimension(max(1,lwork))
                            //lwork = length of array work
                            //info integer(array)
        System.out.println("dsygv");
        printMatrix("  matrix A :", matrix_layout, a, n, n);
        printMatrix("  matrix B :", matrix_layout, b, n, n);
        jniLAPACK.dsygv(matrix_layout, itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
        printMatrix("\n Eigenvectors = ", matrix_layout, a, n, n);
        printMatrix("\n Eigenvalues = ", matrix_layout, w, 1, n);
        
    }
    
    /** Print the matrix X. */
    private static void printMatrix(String prompt, int layout, double[] X, int I, int J) {
        System.out.println(prompt);
        if (layout == jniLAPACK.LAYOUT.ColMajor) {
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[j*I+i]));
                System.out.println();
            }
        }
        else if (layout == jniLAPACK.LAYOUT.RowMajor){
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[i*J+j]));
                System.out.println();
            }
        }
        else{System.out.println("** Illegal layout setting");}
    }
    
    private static void printIntArray(String prompt, int[] X, int L) {
        System.out.println(prompt);
        for (int i=0; i<L; i++) {
                System.out.print("\t" + string(X[i]));
        }
        System.out.println();
    }
    
    /** Shorter string for real number. */
    private static String string(double re) {
        String s="";
        if (re == (long)re)
            s += (long)re;
        else
            s += re;
        return s;
    }
}







