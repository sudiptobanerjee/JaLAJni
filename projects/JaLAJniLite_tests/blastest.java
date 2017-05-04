//import jama.;
import JaLAJniLite.*;
//
//  blastest.java
//
//  Created by Lu Zhang on 5/15/16.
//  Copyright © 2016 Lu Zhang. All rights reserved.
//

/* blastest.java */
public final class blastest {
    private blastest() {}
    public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
        System.out.println("###   Parameter Preparation   ###");
        //
        //Parameter for test
        double alpha=2, beta=-1;
        int L = 3;
        int M = 3, N = 3, K = 3;
        int M2 = 2, N2 = 4, K2 = 3;
        int incx = 1, incy = 1;
        //
        // vector for test
        double[] x = new double[] {1, 2, 1};
	double[] cpx = new double[] {1, 2, 1};
        double[] y = new double[] {1, 2, 3};
        double[] cpy = new double[] {1, 2, 3};
        double[] y2 = new double[] {1, 2};
        double[] cpy2 = new double[] {1, 2};
        //
        // Square matrix in col-major and row-major for test
        double[] As = new double[] {12, -51, 4, 6, 167, -68, -4, 24, -41};
        double[] Bs = new double[] {4, -2, -6, -2, 10, 9, -6, 9, 14};
        double[] Cs = new double[] {2, -4, 3, 7, -3, -7, 0, 18, -20};
        double[] Asc = new double[] {12, 6, -4, -51, 167, 24, 4, -68, -41};
        double[] Csc = new double[] {2, 7, 0, -4, -3, 18, 3, -7, -20};
        //
	double[] cpAs = new double[] {12, -51, 4, 6, 167, -68, -4, 24, -41};
        double[] cpBs = new double[] {4, -2, -6, -2, 10, 9, -6, 9, 14};
        double[] cpCs = new double[] {2, -4, 3, 7, -3, -7, 0, 18, -20};
        double[] cpAsc = new double[] {12, 6, -4, -51, 167, 24, 4, -68, -41};
        double[] cpCsc = new double[] {2, 7, 0, -4, -3, 18, 3, -7, -20};
        
        // rectangular matrix in col-major and row-major for test
        double[] Ac = new double[] {1, 4, 2, 5, 3, 6};
        double[] Bc = new double[] {7, 11, 15, 8, 12, 16, 9, 13, 17, 10, 14, 18};
        double[] Cc = new double[] {2, -3, -4, -7, 3, 0, 7, 18};
        double[] Ec = new double[] {1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12};
        double[] Gc = new double[] {4, 6, 2, 5, 3, 1};
        double[] Ar = new double[] {1, 2, 3, 4, 5, 6};
        double[] Br = new double[] {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
        double[] Cr = new double[] {2, -4, 3, 7, -3, -7, 0, 18};
        double[] Er = new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        double[] Fr = new double[] {1, 2, 3, 4, 5, 6, 7, 8};
        double[] Gr = new double[] {4, 2, 3, 6, 5, 1};
        double[] Br2 = new double[] {7, 8, 9, 10, 11, 12, 13, 14};
        double[] Bc2 = new double[] {7, 11, 8, 12, 9, 13, 10, 14};
        double[] C8 = new double[8];
        double[] C12 = new double[12];
        //
	double[] cpAc = new double[] {1, 4, 2, 5, 3, 6};
        double[] cpBc = new double[] {7, 11, 15, 8, 12, 16, 9, 13, 17, 10, 14, 18};
        double[] cpCc = new double[] {2, -3, -4, -7, 3, 0, 7, 18};
        double[] cpEc = new double[] {1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12};
        double[] cpGc = new double[] {4, 6, 2, 5, 3, 1};
        double[] cpAr = new double[] {1, 2, 3, 4, 5, 6};
        double[] cpBr = new double[] {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
        double[] cpCr = new double[] {2, -4, 3, 7, -3, -7, 0, 18};
        double[] cpEr = new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        double[] cpFr = new double[] {1, 2, 3, 4, 5, 6, 7, 8};
        double[] cpGr = new double[] {4, 2, 3, 6, 5, 1};
        double[] cpBr2 = new double[] {7, 8, 9, 10, 11, 12, 13, 14};
        double[] cpBc2 = new double[] {7, 11, 8, 12, 9, 13, 10, 14};
        double[] cpC8 = new double[8];
        double[] cpC12 = new double[12];
        
        // Working space set for saving results
        // vector
        double[] z = new double[3];
        double[] z2 = new double[2];
        // matrix
        double[] D6 = new double[6];
        double[] D9 = new double[9];
        double[] D8 = new double[8];
        double[] D12 = new double[12];

	double[] cpz = new double[3];
        double[] cpz2 = new double[2];
        // matrix
        double[] cpD6 = new double[6];
        double[] cpD9 = new double[9];
        double[] cpD8 = new double[8];
        double[] cpD12 = new double[12];
        
        //
        // Print the parameters
        //
        printMatrix("Array x", jniBLAS.LAYOUT.RowMajor, x, L, 1);
        printMatrix("Array y", jniBLAS.LAYOUT.RowMajor, y, L, 1);
        printMatrix("Array y2", jniBLAS.LAYOUT.RowMajor, y2, M2, 1);
        System.out.println();
        System.out.println("alpha=" + string(alpha));
        System.out.println("beta=" + string(beta));
        System.out.println();
        System.out.println("Square matrix:");
        printMatrix("Matrix As:", jniBLAS.LAYOUT.RowMajor, As, M, K);
        printMatrix("Matrix Bs:", jniBLAS.LAYOUT.RowMajor, Bs, K, N);
        printMatrix("Matrix Cs:", jniBLAS.LAYOUT.RowMajor, Cs, M, N);
        System.out.println();
        System.out.println("Rectangle Matrix:");
        printMatrix("Matrix A (2 by 3):", jniBLAS.LAYOUT.RowMajor, Ar, M2, K2);
        printMatrix("Matrix B (3 by 4):", jniBLAS.LAYOUT.RowMajor, Br, K2, N2);
        printMatrix("Matrix C (2 by 4):", jniBLAS.LAYOUT.RowMajor, Cr, M2, N2);
        printMatrix("Matrix E (3 by 4):", jniBLAS.LAYOUT.RowMajor, Er, K2, N2);
        printMatrix("Matrix F (4 by 2):", jniBLAS.LAYOUT.RowMajor, Fr, N2, M2);
        printMatrix("Matrix G (2 by 3):", jniBLAS.LAYOUT.RowMajor, Gr, M2, K2);
        //
        // Set Option
        //
        int Layout = jniBLAS.LAYOUT.RowMajor;
        int Layout2 = jniBLAS.LAYOUT.ColMajor;
        int TransA = jniBLAS.TRANSPOSE.NoTrans;
        int TransB = jniBLAS.TRANSPOSE.NoTrans;
        int Uplo = jniBLAS.UPLO.Upper;
        int Diag = jniBLAS.DIAG.NonUnit;
        int Side = jniBLAS.SIDE.Left;
        //
        //
        /* Level 1 */
        //
        System.out.println();
        System.out.println("###   Level 1   ###");
        //
        //
        System.out.println();
        System.out.println("dscal: x = alpha * x");

        int n = L; //n = the length of vector x
	//alpha is constant
	x = cpx.clone();//x is vector
	//incx is increment	
	jniBLAS.dscal(n, alpha, x, incx);

        printMatrix("Resulting x", jniBLAS.LAYOUT.RowMajor, x, n, 1);
        //
        System.out.println();
        System.out.println("daxpy: y = alpha * x + y");
	//n is length of x and y
	//alpha is constant
	x = cpx.clone();//x is vector
	//incx is increment of x
	y = cpy.clone();//y is vector
	//incy is increment of y
        jniBLAS.daxpy(n, alpha, x, incx, y, incy);
        printMatrix("Resulting y", jniBLAS.LAYOUT.RowMajor, y, n, 1);
        //
        /* ddot: A*B */
        //
        System.out.println();
	y = cpy.clone();
        System.out.println("ddot: x dot y = "+ string(jniBLAS.ddot(n, x, incx, y, incy)));
        //
        //
        /* Level 2*/
        //
        System.out.println();
        System.out.println("###   Level 2   ###");
        //
        /* dgemv: y := alpha*A*x + beta*y or y := alpha*A**T*x + beta*y */
        //
        System.out.println();
        System.out.println("dgemv: y2 := alpha*A*x + beta*y2 ");
        z2 = cpy2.clone();
        System.out.println("Row-major: ");
	
	int layout = Layout;//layout specifies layout of matrix
	int trans = TransA;//trans specifies the operation to be performed
	int m = M2;//m = number of rows of matrix a
	n = K2;//n = number of columns of matrix a
	//alpha	is scalar
	double[] a = Ar;//a double precision array, dimension(m,n)
	//x vector
	//incx specifies the increment for the elements of x
	//beta is scalar
	y = z2;//y double precision array
	//incy specifies the increment for the elements of y
        jniBLAS.dgemv(layout, trans, m, n, alpha, a, x, incx, beta, y, incy);
        printMatrix("Resulting y", layout, y, m, 1);
        System.out.println("Col-major: ");
        z2 = cpy2.clone();

	layout = Layout2;
	a = Ac;
	y = z2;
        jniBLAS.dgemv(layout, trans, m, n, alpha, a, x, incx, beta, y, incy);
        printMatrix("Resulting y", layout, y, m, 1);
        //
        /* dtrmv: x := A*x, or x := A**T*x */
        /*A is an n by n unit, or non-unit, upper or lower triangular matrix*/
        //
        System.out.println();
        System.out.println("dtrmv: x := As * x");
        System.out.println("As is treated as a non-unit, upper triangular matrix");
        System.out.println("Row-major: ");
        
	layout = Layout;
	int uplo = Uplo;//specifies whether the matrix is upper or lower triangular part
	//trans specifies the operation to be performed
	int diag = Diag;//diag specifies whether or not a is unit triangular
	n = N;//n = the order of matrix a
	a = As;//a double precision array, dimension(n,n)
	x = cpx.clone();//x double precision array
	//incx specifies the increment for the elements of x
        printMatrix("Initial x", layout, x, n, 1);
        jniBLAS.dtrmv(layout, uplo, trans, diag, n, a, x, incx);
        printMatrix("Resulting x", layout, x, n, 1);
        System.out.println("Col-major: ");
        z = cpx.clone();

	layout = Layout2;
	a = Asc;
	x = cpx.clone();
        jniBLAS.dtrmv(layout, uplo, trans, diag, n, a, x, incx);
        printMatrix("Resulting z", layout, x, n, 1);
        //
        //
        /* dsymv:  y := alpha*A*x + beta*y */
        /* Where A is a n by n symmetric matrix */
        //
        System.out.println();
        System.out.println("dsymv: y := alpha*Bs*x + beta*y; Bs is a symmetric matrix");
        z = cpy.clone();
        System.out.println("Row-major:");

	layout = Layout;
	//uplo specifies which part of a is specified
	//n = the order of matrix a
	//alpha is scalar
	a = Bs;//a double precision array
	x = cpx.clone();//x double precision array
	//incx specifies the increment for the elements of x
	//beta is scalar
	y = z;//y double precision array
	//incy specifies the increment for the elements of y
        jniBLAS.dsymv(layout, uplo, n, alpha, a, x, incx, beta, y, incy);
        printMatrix("Resulting y", layout, y, n, 1);

        z = cpy.clone();
        System.out.println("Col-major:");
	layout = Layout2;
	uplo = Uplo;//uplo specifies which part of a is specified
	n = N;//n = the order of matrix a
	//alpha is scalar
	a = Bs;//a is double precision array
	//x double precision array
	//incx specifies the increment for the elements of x
	y = z;//y double precision array
	//incy specifies the increment for the elements of y
        jniBLAS.dsymv(layout, uplo, n, alpha, a, x, incx, beta, y, incy);
        printMatrix("Resulting y", layout, y, n, 1);
        //
        //
        /* Level 3 */
        //
        System.out.println();
        System.out.println("###   Level 3   ###");
        //
        //
        /* dgemm: C := alpha*op( A )*op( B ) + beta*C */
        /* where op( X ) = X   or   op( X ) = X**T */
        //
        System.out.println("dgemm: C := alpha*op( A )*op( B ) + beta*C \n where op( X ) = X   or   op( X ) = X**T");
        TransA = jniBLAS.TRANSPOSE.NoTrans;
        TransB = jniBLAS.TRANSPOSE.NoTrans;
        System.out.println("\n A NoTrans, B NoTrans");
        System.out.println("dgemm: C := alpha * A * B + beta * C");
        System.out.println("Row-major: ");
        D8 = Cr.clone();
	
	layout = Layout;
	int transA = TransA; //transA specifies the form of op(a) to be used in the matrix multiplication
	int transB = TransB; //transB specifies the form of op(b) to be used in matrix multiplication
	m = M2;//m = number of rows of matrix op(a)
	n = N2;//n = number of columns of matrix op(b)
	int k = K2;//k = number of columns of matrix op(a)
	//alpha is scalar
	a = Ar;//a is double precision array
	double[] b = Br;//b is double precision array
	//beta is scalar
	double[] c = D8;//c is double precision array
        jniBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        System.out.println("Col-major: ");
        D8 = cpCc.clone();

	layout = Layout2;
	a = Ac;
	b = Bc;
	c = D8;
        jniBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        System.out.println();
        System.out.println("Matrix Multiplication in Row-major ");
        D8 = cpCr.clone();

	layout = Layout;
	transA = TransA;
	transB = TransB;
	m = M2;
	n = N2;
	k = K2;
	alpha = 1;
	a = Ar;
	b = Br;
	beta = 0;
	c = D8;
        jniBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting A * B", layout, c, m, n);
	
	m = M2;
	n = M2;
	k = N2;
	a = c;
	b = Fr;
        jniBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting ABF", Layout, c, m, n);
        printMatrix("Result is saved in row-major", Layout2, c, 1, m * n);
        System.out.println();
        
        System.out.println();
        transA = jniBLAS.TRANSPOSE.Trans;
        transB = jniBLAS.TRANSPOSE.NoTrans;
        System.out.println("A Trans, B NoTrans");
        System.out.println("Row-major: dgemm: C := alpha * t(A) * (B2) ");
        D12 = cpC12.clone();
	
	layout = Layout;
	m = K2;
	n = N2;
	k = M2;
	alpha = 2;
	a = Ar;
	b = Br2;
	beta = -1;
	c = D12;
        jniBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        System.out.println();
        System.out.println("Col-major: dgemm: C := alpha * t(A) * (B2) ");

        D12 = cpC12.clone();
	layout = Layout2;
	a = Ac;
	b = Bc2;
	c = D12;
        jniBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        
        System.out.println();
        transA = jniBLAS.TRANSPOSE.Trans;
        transB = jniBLAS.TRANSPOSE.Trans;
        System.out.println("A Trans, B Trans");
        System.out.println("Row-major: dgemm: C := alpha * t(B) * t(A)");
        
	D8 = cpC8.clone();
	layout = Layout;
	m = N2;
	n = M2;
	k = K2;
	a = Br;
	b = Ar;
	c = D8;
        jniBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        System.out.println();
        System.out.println("Col-major: dgemm: C := alpha * t(B) * t(A)");
        
	D8 = cpC8.clone();
	layout = Layout2;
	m = N2;
	n = M2;
	k = K2;
	a = Bc;
	b = Ac;
	c = D8;
        jniBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        //
        /* Options for dtrmm and dsymm */
        TransA = jniBLAS.TRANSPOSE.NoTrans;
        Uplo = jniBLAS.UPLO.Upper;
        //
        /* dtrmm: B := alpha*op( A )*B,   or   B := alpha*B*op( A ) */
        /* A  is a unit, or non-unit,  upper or lower triangular matrix */
        /* op( A ): op( A ) = A   or   op( A ) = A**T*/
        //
        System.out.println();
        System.out.println("dtrmm: B := alpha*op( A )*B,   or   B := alpha*B*op( A ) ");
        System.out.println("A is treated as a non-unit upper triangular matrix");
        System.out.println("dtrmm: B := alpha* As * B");
        System.out.println("Row-major: ");
        
	D12 = cpBr.clone();
	layout = Layout;
	int side = Side;//side specifies whether op(a) multiplies b from the left or right as follows 
	uplo = Uplo;//uplo specifies whether a is an upper or lower triangular matrix
	transA = TransA;//transA specifies the form of op(a)
	//diag specifies whether or not a is unit triangular
	m = M;//m = the number of rows of b
	n = N2;//n = the number of columns of b
	//alpha is scalar
	a = As;//a double precision array
	b = D12;//b double precision array
        jniBLAS.dtrmm(layout, side, uplo, transA, diag, m, n, alpha, a, b);
        printMatrix("Resulting B", layout, b, m, n);
        System.out.println("Col-major: ");
        
	D12 = cpBc.clone();
	layout = Layout2;
	a = Asc;
	b = D12;
        jniBLAS.dtrmm(layout, side, uplo, transA, diag, m, n, alpha, a, b);
        printMatrix("Resulting B", layout, b, m, n);
        System.out.println("\n Test Side = Right: dtrmm: A := alpha * A * As");
        Side = jniBLAS.SIDE.Right;
        System.out.println("Row-major: ");
        
	D6 = cpAr.clone();
	layout = Layout;
	side = Side;
	m = M2;
	n = N;
	a = As;
	b = D6;
        jniBLAS.dtrmm(layout, side, uplo, transA, diag, m, n, alpha, a, b);
        printMatrix("Resulting A", layout, b, m, n);
        System.out.println("Col-major: ");
        
	D6 = cpAc.clone();
	layout = Layout2;
	a = Asc;
	b = D6;
        jniBLAS.dtrmm(layout, side, uplo, transA, diag, m, n, alpha, a, b);
        printMatrix("Resulting A", layout, b, m, n);
        //
        //
        /* dsymm: C := alpha*A*B + beta*C or C := alpha*B*A + beta*C */
        /* A is a symmetric matrix*/
        //
        Side = jniBLAS.SIDE.Left;
        System.out.println();
        System.out.println("dsymm: E := alpha*Bs*B + beta*E; Bs is a symmetric matrix");
        System.out.println("Row-major: ");
        D12 = cpEr.clone();
	layout = Layout;
	side = Side;//side specifies whether the symmetric matrix a appears on the left or right in the operation
	//uplo specifies whether the upper or lower triangular part of the symmetric matrix a is to be referenced	
	m = M;//m = the number of rows of matrix c
	n = N2;//n = number of columns of matrix c
	//alpha is scalar	
	a = Bs;//a double precision array
	b = Br;//b double precision array
	c = D12;//c double precision array
        jniBLAS.dsymm(layout, side, uplo, m, n, alpha, a, b, beta, c);
        printMatrix("Resulting E", layout, c, m, n);
        System.out.println("Col-major: ");
        
	D12 = cpEc.clone();
	layout = Layout2;
	a = Bs;
	b = Bc;
	c = D12;
        jniBLAS.dsymm(layout, side, uplo, m, n, alpha, a, b, beta, c);
        printMatrix("Resulting E", layout, c, m, n);
        
        System.out.println("\n Test Side = Right: dsymm: G := alpha*A*Bs + beta*G;");
        Side = jniBLAS.SIDE.Right;
        System.out.println("Row-major: ");
        
	D6 = cpGr.clone();
	layout = Layout;
	side = Side;
	m = M2;
	n = K2;
	a = Bs;
	b = Ar;
	c = D6;
        jniBLAS.dsymm(layout, side, uplo, m, n, alpha, a, b, beta, c);
        printMatrix("Resulting G", layout, c, m, n);
        System.out.println("Col-major: ");
        
	D6 = cpGc.clone();
	layout = Layout2;
	a = Bs;
	b = Ac;
	c = D6;
        jniBLAS.dsymm(layout, side, uplo, m, n, alpha, a, b, beta, c);
        printMatrix("Resulting G", layout, c, m, n);
    }
    
    /** Print the matrix X. */
    private static void printMatrix(String prompt, int layout, double[] X, int I, int J) {
        System.out.println(prompt);
        if (layout == jniBLAS.LAYOUT.ColMajor) {
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[j*I+i]));
                System.out.println();
            }
        }
        else if (layout == jniBLAS.LAYOUT.RowMajor){
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[i*J+j]));
                System.out.println();
            }
        }
        else{System.out.println("** Illegal layout setting");}
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







