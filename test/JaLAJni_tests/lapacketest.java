//
//  main.c
//  lapacketest.java
//
//  Created by Lu Zhang on 5/15/16.
//  Copyright © 2016 Lu Zhang. All rights reserved.
//
import JaLAJni.*;
/* lapacketest.java */
public final class lapacketest {
    private lapacketest() {}
    public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
        System.out.println("###   Testfile for jni_lapacke   ### \n \n");
        System.out.println("###   Parameter Preparation   ###");
        int M=3, N=3, K=3;
        int info = 0;
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
        //
        // Print the parameters
        //
        printMatrix("Matrix A", jniLAPACKE.LAYOUT.RowMajor, A, M, K);
        printMatrix("Matrix B", jniLAPACKE.LAYOUT.RowMajor, B, K, N);
        //
        // Set Option
        //
        int matrix_layout = jniLAPACKE.LAYOUT.RowMajor;
        char trans = jniLAPACKE.TRANSPOSE.NoTrans;
        char Uplo = jniLAPACKE.UPLO.Upper;
        char Jobvl = jniLAPACKE.JOBV.Compute;
        char Jobvr = jniLAPACKE.JOBV.Compute;
        char Jobu = jniLAPACKE.JOB.All;
        char Jobvt = jniLAPACKE.JOB.All;
        char Jobz = jniLAPACKE.JOB.Overwritten;
	int itype = jniLAPACKE.ITYPE.first;
        //
        //
        /* ---- LU ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   LU   ###");
        //
        /* dgetrf */
        //
        System.out.println();
        System.out.println("The LU decomposition A=PLU (L has unit diagnoal elements):");
	D = A.clone();
        int m = M;//m = the number of rows of the matrix a
	int n = N;//n = the number of columns of the matrix a
	double[] a = D;//a double precision, dimension(m,n)
	int lda = N;//lda = leading dimension of array a
	//ipiv = pivot indices       
	jniLAPACKE.dgetrf( matrix_layout, m, n, a, lda, ipiv);
        printMatrix("  The LU matrix of A:",matrix_layout, a, m, n);
        printIntArray("  The permutation vector ipiv: ", ipiv, n);
        //
        /* dgetrs */
        //
        System.out.println();
        D = B.clone();
        D2 = A.clone();
	a = D;
        jniLAPACKE.dgetrf(matrix_layout, m, n, a, lda, ipiv);
	
	//trans specifies the form of the system of equations
	n = M;//n = order of matrix a
	int nrhs = N;//nrhs = number of columns of matrix b 
	//a is computed by dgetrf
	lda = n;//lda = leading dimension of a
	//ipiv = pivot indices from dgetrf
	double[] b = D2;//b double precision, dimension(n,nrhs)
	int ldb = nrhs;//ldb = leading dimension of b
        jniLAPACKE.dgetrs(matrix_layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
        printMatrix("  The solution  matrix of BX=A by LU is:", matrix_layout, b, n, nrhs);
        //
        /* dgetri */
        //
        System.out.println();
        D = A.clone();
	
	m = M;
	n = N;
	a = D;
	lda = M;
        jniLAPACKE.dgetrf(matrix_layout, m, n, a, lda, ipiv);
	
	n = M;//n = order of matrix a
	a = D;//a is computed by dgetrf, double precision
	lda = n;//lda = leading dimension of a
	//ipiv is pivot indices from dgetrf       
	jniLAPACKE.dgetri(matrix_layout, n, a, lda, ipiv);
        printMatrix("  The Inverse of  matrix A is :", matrix_layout, a, n, n);
        //
        /* ---- Cholesky ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   Cholesky   ###");
        //
        /* dpotrf */
        //
        System.out.println();
        D = B.clone();
	
	char uplo = Uplo;//uplo specifies which part of a is stored
	n = M;//n = order of matrix a
	a = D;//a is symmetric matrix
	lda = M;//lda = leading dimension of a
        jniLAPACKE.dpotrf(matrix_layout, uplo, n, a, lda);
        printMatrix("  The Cholesky factorization of matrix B is :", matrix_layout, a, n, n);
        //
        /* dpotri */
        //
        System.out.println();
        jniLAPACKE.dpotri(matrix_layout, Uplo, M, D, M);
        printMatrix("  The upper triangular part of Inverse of matrix B is :", matrix_layout, D, M, N);
        //
        /*  dpotrs */
        //
        System.out.println();
        D = B.clone();
        D2 = A.clone();
	
	//uplo specifies which part of a is stored
	n = M;//n = order of matrix a
	a = D;//a is symmetric matrix
	lda = n;//lda = leading dimension of a
        jniLAPACKE.dpotrf(matrix_layout, uplo, n, a, lda);

	//uplo specifies which part of a is stored
	n = M;//n = order of matrix a
	nrhs = N;//nrhs = number of columns of matrix b
	a = D;//a is computed by dpotrf
	lda = n;//lda = leading dimension of matrix a
	b = D2;//b = right hand matrix, dimension(n,nrhs)
	ldb = nrhs;//ldb = leading dimension of matrix b
        jniLAPACKE.dpotrs(matrix_layout, uplo, n, nrhs, a, lda, b, ldb);
        printMatrix("  The solution  matrix of BX=A by Cholesky is :", matrix_layout, b, n, nrhs);
        //
        /* ---- QR ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   QR decomposition   ###");
        //
        /* dgeqrf */
        //
        System.out.println();
        System.out.println("A=QR:");
        D = A.clone();

	m = M;//m = number of rows of matrix a
	n = N;//n = number of columns of matrix a
	a = D;//a double precision, dimension(m,n)
	lda = n;//lda = leading dimension of matrix a
	//tau = scalar factors of the elementary reflectors
        jniLAPACKE.dgeqrf(matrix_layout, m, n, a, lda, tau);
        printMatrix("  The rewritten matrix A is :", matrix_layout, a, m, n);
        //
        /* dorgqr */
        //
        System.out.println();
	m = M;//m = number of rows of matrix q
	n = N;//n = number of columns of matrix q
	int k = M;//k = number of elementary reflectors whose product defines the matrix q
	a = D;//a is returned by dgeqrf
	lda = M;//lda = leading dimension of a
	//tau is returned by dgeqrf
        jniLAPACKE.dorgqr(matrix_layout, m, n, k, a, lda, tau);
        printMatrix("  The matrix Q is :", matrix_layout, a, m, n);
        //
        /* dgeqp3 */
        //
        System.out.println();
        D = A.clone();
	
	m = M;//m = number of rows of matrix a
	n = N;//n = number of columns of matrix a
	a = D;//a double precision
	lda = M;//lda = leading dimension of a
	//jpvt integer array, dimension(n)
	//tau = scalar factors of the elementary reflectors, dimension(min(m,n))
        jniLAPACKE.dgeqp3(matrix_layout, m, n, a, lda, jpvt, tau);
        printMatrix("  The R matrix computed by dgeqp3:", matrix_layout, a, m, n);
        //
        /* ---- Eigenvector and SVD ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   Eigenvector and SVD   ###");
        //
        /* dgeev */
        //
        System.out.println();
        System.out.println(" \nGet Eigenvectors of A:");
        D = A.clone();

	char jobvl = Jobvl;//jobvl specifies whether left eigenvectors of a are computed
	char jobvr = Jobvr;//jobvr specifies whether right eigenvectors of a are computed
	n = N;//n = order of matrix a
	a = D;//a double precision
	lda = N;//lda = leading dimension of a
	double[] wr = WR;//wr = real part of eigenvalues
	double[] wi = WI;//wi = image part of eigenvalues
	double[] vl = VL;//vl = left eigenvectors
	int ldvl = N;//ldvl = leading dimension of vl
	double[] vr = VR;//vr = right eigenvectors
	int ldvr = N;//ldvr = leading dimension of vr
        jniLAPACKE.dgeev(matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
        printMatrix("  The overwritten matrix A :", matrix_layout, a, n, n);
        printMatrix("\n  Left eigenvectors (in columns) :", matrix_layout, vl, n, ldvl);
        printMatrix("\n  Right eigenvectors (in columns):", matrix_layout, vr, n, ldvr);
        //
        /* dgesvd */
        //
        System.out.println();
        System.out.println(" \nSVD of A: A = U * SIGMA * transpose(V):");
        D = A.clone();
	
	char jobu = Jobu;//jobu specifies options for computing all or part of the matrix u
	char jobvt = Jobvt;//jobvt specifies options for computing all or part of the matrix vt
	m = M;//m = number of rows of input matrix a
	n = N;//n = number of columns of input matrix a
	a = D;//a double precision
	lda = n;//lda = leading dimension of matri a
	double[] s = S;//s = singular values of a
	double[] u = U;//u = U
	int ldu = M;//ldu = leading dimension of u
	double[] vt = VT;//vt = VT;
	int ldvt = N;//ldvt = leading dimension of vt
	jniLAPACKE.dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
        printMatrix("  The overwritten matrix A :", matrix_layout, a, m, n);
        printMatrix("\n  The singular values of A :", matrix_layout, s, 1, m);
        printMatrix("\n  Left singular vectors (in columns) :", matrix_layout, u, m, n);
        printMatrix("\n  Right singular vectors (in rows):", matrix_layout, vt, n, n);
        //
        /* dgesdd */
        System.out.println();
        System.out.println(" \nSVD of A calculated by dgesdd:");
        D = A.clone();

	jobu = Jobu;//jobu specifies options for computing all or part of the matrix u 
	m = M;//m = number of rows of input matrix a
	n = N;//n = number of columns of input matrix a
	a = D;//a double precision
	lda = n;//lda = leading dimension of matrix a
	s = S;//s = singular values of a
	u = U;//u = U
	ldu = M;//ldu = leading dimension of u
	vt = VT;//vt = VT;
	ldvt = N;//ldvt = leading dimension of vt
        jniLAPACKE.dgesdd(matrix_layout, jobu, m, n, a, lda, s, u, ldu, vt, ldvt);
        printMatrix("  The overwritten matrix A :", matrix_layout, a, m, n);
        printMatrix("\n  The singular values of A :", matrix_layout, s, 1, m);
        printMatrix("\n  Left singular vectors (in columns) :", matrix_layout, u, m, n);
        printMatrix("\n  Right singular vectors (in rows):", matrix_layout, vt, n, n);
        //

	//dsyev
	System.out.println();
        System.out.println(" \ndsyev: computes the eigenvalues and right/left eigenvectors for symmetric matrix A");

	char jobz = Jobvl;//jobz specifies whether to compute eigenvectors or not
	//uplo specifies which part of a is stored
	n = N;//n = order of matrix a
	a = sym;//a is symmetric matrix
	lda = N;//lda = leading dimension of matrix a
	double[] w = S;//w = eigenvalues in ascendign order

        printMatrix("  matrix A :", matrix_layout, a, n, n);
        jniLAPACKE.dsyev(matrix_layout, jobz, uplo, n, a, lda, w);
        printMatrix("  The overwritten matrix A :", matrix_layout, a, n, n);
        printMatrix("\n  The singular values of A :", matrix_layout, w, 1, n);

	//dgesv
	System.out.println();
        System.out.println(" \ndgesv: computes the solution to linear equation A*X = B for GE matrices");

	n = N;//n = order of matrix a
	nrhs = 1;//nrhs = number of columns of matrix b
	a = A;//a double precision
	lda = N;//lda = leading dimension of a
	//ipiv = pivot indices that define the permutation matrix P
	b = eq;//b = right hand side matrix, dimension(n,nrhs)
	ldb = 1;//ldb = leading dimension of matrix b

        printMatrix("  matrix A :", matrix_layout, a, n, n);
        printMatrix("  matrix B :", matrix_layout, b, n, nrhs);
        jniLAPACKE.dgesv(matrix_layout, n, nrhs, a, lda, ipiv, b, ldb);
        printMatrix("  X = ", matrix_layout, b, n, nrhs);
        
        //dggev
	System.out.println();
        System.out.println(" \ndggev: computes the eigenvalues and right/left eigenvectors for symmetric matrix A");
	sym = cpsym.clone();
        
	//jobvl specifies whether the left generalized eigenvectors are computed
	//jobvr specifies whether the right generalized wigenvectors are computed
	n = N;//n = order of matrices a, b, vl, vr
	a = cpA;//a is matrix a in the pair(a,b)
	lda = N;//lda = leading dimension of matrix a
	b = sym;//b is matrix b in the pari(a,b)
	ldb = N;//ldb = leading dimension of matrix b
	double[] alphar = S;//alphar dimension(n)
	double[] alphai = eq;//alphai dimension(n)
	double[] beta = tau;//beta double precision, dimension(n)
	vl = VL;//vl = left eigenvectors
	ldvl = M;//ldvl = leading dimension of vl
	vr = VR;//vr = right eigenvectors
	ldvr = M;//ldvr = leading dimension of vr

        printMatrix("  matrix A :", matrix_layout, a, n, n);
	printMatrix("  matrix B :", matrix_layout, b, n, n);
        jniLAPACKE.dggev(matrix_layout, jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr);
        printMatrix("\n alphar = ", matrix_layout, alphar, 1, n);
        printMatrix("\n alphai = ", matrix_layout, alphai, 1, n);
        printMatrix("\n beta = ", matrix_layout, beta, 1, n);
	printMatrix("\n VL = ", matrix_layout, vl, n, n);
	printMatrix("\n VR = ", matrix_layout, vr, n, n);

	//dsygv
	System.out.println();
        System.out.println(" \ndsygv: computes the eigenvalues and right/left eigenvectors for real-generalized symmetric definite eigenproblem");
	sym = cpsym.clone();

	//itype specifies the problem type to be solved
	jobz = Jobvl;//jobz specifies whether to compute eigenvectors or not
	//uplo specifies which part of a and b are stored
	n = N;//n = order the matrices a and b
	a = sym;//a is symmetric matrix
	lda = N;//lda = leading dimension of a
	b = I;//b double precision
	ldb = N;//ldb = leading dimension of matrix b
	w = S;//w double precision, dimension(max(1,lwork))

        printMatrix("  matrix A :", matrix_layout, a, n, n);
        printMatrix("  matrix B :", matrix_layout, b, n, n);
        jniLAPACKE.dsygv(matrix_layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
        printMatrix("\n A = ", matrix_layout, a, n, n);
        printMatrix("\n W = ", matrix_layout, w, 1, n);
        
    }
    
    /** Print the matrix X. */
    private static void printMatrix(String prompt, int layout, double[] X, int I, int J) {
        System.out.println(prompt);
        if (layout == jniLAPACKE.LAYOUT.ColMajor) {
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[j*I+i]));
                System.out.println();
            }
        }
        else if (layout == jniLAPACKE.LAYOUT.RowMajor){
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







