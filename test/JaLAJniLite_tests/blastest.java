
import JaLAJniLite.jni_blas.*;
//
//  blastest.java
//
//  Created by Lu Zhang on 5/15/16.
//  Copyright © 2016 Lu Zhang. All rights reserved.
//

public final class blastest {
    private blastest() {}
    public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
        System.out.println("###   Parameter Preparation   ###");
        //
        // Parameter for test
        //
        double alpha=2, beta=-1;
        int L = 3;
        int M = 3, N = 3, K = 3;
        int M2 = 2, N2 = 4, K2 = 3;
        int incx = 1, incy = 1;
        //
        // vector for test
        //
        double[] x = new double[] {1, 2, 1};
        double[] cpx = new double[] {1, 2, 1};
        double[] y = new double[] {1, 2, 3};
        double[] cpy = new double[] {1, 2, 3};
        double[] y2 = new double[] {1, 2};
        double[] cpy2 = new double[] {1, 2};
        //
        // Square matrix in col-major and row-major for test
        //
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
        //
        // rectangular matrix in col-major and row-major for test
        //
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
        double[] cpAc = new double[] {1, 4, 2, 5, 3, 6};
        double[] cpBc = new double[] {7, 11, 15, 8, 12, 16, 9, 13, 17,
            10, 14, 18};
        double[] cpCc = new double[] {2, -3, -4, -7, 3, 0, 7, 18};
        double[] cpEc = new double[] {1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12};
        double[] cpGc = new double[] {4, 6, 2, 5, 3, 1};
        double[] cpAr = new double[] {1, 2, 3, 4, 5, 6};
        double[] cpBr = new double[] {7, 8, 9, 10, 11, 12, 13, 14, 15,
            16, 17, 18};
        double[] cpCr = new double[] {2, -4, 3, 7, -3, -7, 0, 18};
        double[] cpEr = new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        double[] cpFr = new double[] {1, 2, 3, 4, 5, 6, 7, 8};
        double[] cpGr = new double[] {4, 2, 3, 6, 5, 1};
        double[] cpBr2 = new double[] {7, 8, 9, 10, 11, 12, 13, 14};
        double[] cpBc2 = new double[] {7, 11, 8, 12, 9, 13, 10, 14};
        double[] cpC8 = new double[8];
        double[] cpC12 = new double[12];
        //
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
        double[] cpD6 = new double[6];
        double[] cpD9 = new double[9];
        double[] cpD8 = new double[8];
        double[] cpD12 = new double[12];
        //
        // Print the parameters
        //
        printMatrix("Array x", jniCBLAS.LAYOUT.RowMajor, x, L, 1);
        printMatrix("Array y", jniCBLAS.LAYOUT.RowMajor, y, L, 1);
        printMatrix("Array y2", jniCBLAS.LAYOUT.RowMajor, y2, M2, 1);
        System.out.println();
        System.out.println("alpha=" + string(alpha));
        System.out.println("beta=" + string(beta));
        System.out.println();
        System.out.println("Square matrix:");
        printMatrix("Matrix As:", jniCBLAS.LAYOUT.RowMajor, As, M, K);
        printMatrix("Matrix Bs:", jniCBLAS.LAYOUT.RowMajor, Bs, K, N);
        printMatrix("Matrix Cs:", jniCBLAS.LAYOUT.RowMajor, Cs, M, N);
        System.out.println();
        System.out.println("Rectangle Matrix:");
        printMatrix("Matrix A (2 by 3):", jniCBLAS.LAYOUT.RowMajor, Ar, M2, K2);
        printMatrix("Matrix B (3 by 4):", jniCBLAS.LAYOUT.RowMajor, Br, K2, N2);
        printMatrix("Matrix C (2 by 4):", jniCBLAS.LAYOUT.RowMajor, Cr, M2, N2);
        printMatrix("Matrix E (3 by 4):", jniCBLAS.LAYOUT.RowMajor, Er, K2, N2);
        printMatrix("Matrix F (4 by 2):", jniCBLAS.LAYOUT.RowMajor, Fr, N2, M2);
        printMatrix("Matrix G (2 by 3):", jniCBLAS.LAYOUT.RowMajor, Gr, M2, K2);
        //
        // Set Option
        //
        int Layout = jniCBLAS.LAYOUT.RowMajor;
        int Layout2 = jniCBLAS.LAYOUT.ColMajor;
        int TransA = jniCBLAS.TRANSPOSE.NoTrans;
        int TransB = jniCBLAS.TRANSPOSE.NoTrans;
        int Uplo = jniCBLAS.UPLO.Upper;
        int Diag = jniCBLAS.DIAG.NonUnit;
        int Side = jniCBLAS.SIDE.Left;
        //
        /* Level 1 */
        //
        System.out.println();
        System.out.println("###   Level 1   ###");
        System.out.println();
        //
        // dscal
        //
        System.out.println("dscal: x = alpha * x");
        int n = L;
        x = cpx.clone();
        jniCBLAS.dscal(n, alpha, x, incx);
        printMatrix("Resulting x", jniCBLAS.LAYOUT.RowMajor, x, n, 1);
        //
        // daxpy
        //
        System.out.println();
        System.out.println("daxpy: y = alpha * x + y");
        x = cpx.clone();
        y = cpy.clone();
        jniCBLAS.daxpy(n, alpha, x, incx, y, incy);
        printMatrix("Resulting y", jniCBLAS.LAYOUT.RowMajor, y, n, 1);
        //
        // ddot
        //
        System.out.println();
        y = cpy.clone();
        System.out.println("ddot: x dot y = "+
                           string(jniCBLAS.ddot(n, x, incx, y, incy)));
        //
        /* Level 2*/
        //
        System.out.println();
        System.out.println("###   Level 2   ###");
        //
        // dgemv
        //
        System.out.println();
        System.out.println("dgemv: y2 := alpha*A*x + beta*y2 ");
        z2 = cpy2.clone();
        System.out.println("Row-major: ");
        int layout = Layout;
        int trans = TransA;
        int m = M2;
        n = K2;
        double[] a = Ar;
        y = z2;
        jniCBLAS.dgemv(layout, trans, m, n, alpha, a, x, incx, beta, y, incy);
        printMatrix("Resulting y", layout, y, m, 1);
        System.out.println("Col-major: ");
        z2 = cpy2.clone();
        layout = Layout2;
        a = Ac;
        y = z2;
        jniCBLAS.dgemv(layout, trans, m, n, alpha, a, x, incx, beta, y, incy);
        printMatrix("Resulting y", layout, y, m, 1);
        //
        // dtrmv
        //
        System.out.println();
        System.out.println("dtrmv: x := As * x");
        System.out.println("As is treated as a non-unit, upper triangular matrix");
        System.out.println("Row-major: ");
        layout = Layout;
        int uplo = Uplo;
        int diag = Diag;
        n = N;
        a = As;
        x = cpx.clone();
        printMatrix("Initial x", layout, x, n, 1);
        jniCBLAS.dtrmv(layout, uplo, trans, diag, n, a, x, incx);
        printMatrix("Resulting x", layout, x, n, 1);
        System.out.println("Col-major: ");
        z = cpx.clone();
        layout = Layout2;
        a = Asc;
        x = cpx.clone();
        jniCBLAS.dtrmv(layout, uplo, trans, diag, n, a, x, incx);
        printMatrix("Resulting z", layout, x, n, 1);
        //
        // dsymv
        //
        System.out.println();
        System.out.println("dsymv: y := alpha*Bs*x + beta*y; Bs is a symmetric matrix");
        z = cpy.clone();
        System.out.println("Row-major:");
        layout = Layout;
        a = Bs;
        x = cpx.clone();
        y = z;
        jniCBLAS.dsymv(layout, uplo, n, alpha, a, x, incx, beta, y, incy);
        printMatrix("Resulting y", layout, y, n, 1);
        z = cpy.clone();
        System.out.println("Col-major:");
        layout = Layout2;
        uplo = Uplo;
        n = N;
        a = Bs;
        y = z;
        jniCBLAS.dsymv(layout, uplo, n, alpha, a, x, incx, beta, y, incy);
        printMatrix("Resulting y", layout, y, n, 1);
        //
        /* Level 3 */
        //
        System.out.println();
        System.out.println("###   Level 3   ###");
        //
        // dgemm
        //
        System.out.println
        ("dgemm: C := alpha*op( A )*op( B ) + beta*C \n where op( X ) = X   or   op( X ) = X**T");
        TransA = jniCBLAS.TRANSPOSE.NoTrans;
        TransB = jniCBLAS.TRANSPOSE.NoTrans;
        System.out.println("\n A NoTrans, B NoTrans");
        System.out.println("dgemm: C := alpha * A * B + beta * C");
        
        System.out.println("Row-major: ");
        D8 = Cr.clone();
        layout = Layout;
        int transA = TransA;
        int transB = TransB;
        m = M2;
        n = N2;
        int k = K2;
        a = Ar;
        double[] b = Br;
        double[] c = D8;
        jniCBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        
        System.out.println("Col-major: ");
        D8 = cpCc.clone();
        layout = Layout2;
        a = Ac;
        b = Bc;
        c = D8;
        jniCBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
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
        jniCBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        
        printMatrix("Resulting A * B", layout, c, m, n);
        m = M2;
        n = M2;
        k = N2;
        a = c;
        b = Fr;
        jniCBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting ABF", Layout, c, m, n);
        printMatrix("Result is saved in row-major", Layout2, c, 1, m * n);
        System.out.println();
        
        System.out.println();
        transA = jniCBLAS.TRANSPOSE.Trans;
        transB = jniCBLAS.TRANSPOSE.NoTrans;
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
        jniCBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        System.out.println();
        
        System.out.println("Col-major: dgemm: C := alpha * t(A) * (B2) ");
        D12 = cpC12.clone();
        layout = Layout2;
        a = Ac;
        b = Bc2;
        c = D12;
        jniCBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        System.out.println();
        
        transA = jniCBLAS.TRANSPOSE.Trans;
        transB = jniCBLAS.TRANSPOSE.Trans;
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
        jniCBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
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
        jniCBLAS.dgemm(layout, transA, transB, m, n, k, alpha, a, b, beta, c);
        printMatrix("Resulting C", layout, c, m, n);
        TransA = jniCBLAS.TRANSPOSE.NoTrans;
        Uplo = jniCBLAS.UPLO.Upper;
        //
        // dtrmm
        //
        System.out.println();
        System.out.println("dtrmm: B := alpha*op( A )*B,   or   B := alpha*B*op( A ) ");
        System.out.println("A is treated as a non-unit upper triangular matrix");
        System.out.println("dtrmm: B := alpha* As * B");
        
        System.out.println("Row-major: ");
        D12 = cpBr.clone();
        layout = Layout;
        int side = Side;
        uplo = Uplo;
        transA = TransA;
        m = M;
        n = N2;
        a = As;
        b = D12;
        jniCBLAS.dtrmm(layout, side, uplo, transA, diag, m, n, alpha, a, b);
        printMatrix("Resulting B", layout, b, m, n);
        
        System.out.println("Col-major: ");
        D12 = cpBc.clone();
        layout = Layout2;
        a = Asc;
        b = D12;
        jniCBLAS.dtrmm(layout, side, uplo, transA, diag, m, n, alpha, a, b);
        printMatrix("Resulting B", layout, b, m, n);
        System.out.println("\n Test Side = Right: dtrmm: A := alpha * A * As");
        Side = jniCBLAS.SIDE.Right;
        
        System.out.println("Row-major: ");
        D6 = cpAr.clone();
        layout = Layout;
        side = Side;
        m = M2;
        n = N;
        a = As;
        b = D6;
        jniCBLAS.dtrmm(layout, side, uplo, transA, diag, m, n, alpha, a, b);
        printMatrix("Resulting A", layout, b, m, n);
        
        System.out.println("Col-major: ");
        D6 = cpAc.clone();
        layout = Layout2;
        a = Asc;
        b = D6;
        jniCBLAS.dtrmm(layout, side, uplo, transA, diag, m, n, alpha, a, b);
        printMatrix("Resulting A", layout, b, m, n);
        //
        // dsymm
        //
        Side = jniCBLAS.SIDE.Left;
        System.out.println();
        System.out.println
        ("dsymm: E := alpha*Bs*B + beta*E; Bs is a symmetric matrix");
        
        System.out.println("Row-major: ");
        D12 = cpEr.clone();
        layout = Layout;
        side = Side;
        m = M;
        n = N2;
        a = Bs;
        b = Br;
        c = D12;
        jniCBLAS.dsymm(layout, side, uplo, m, n, alpha, a, b, beta, c);
        printMatrix("Resulting E", layout, c, m, n);
        
        System.out.println("Col-major: ");
        D12 = cpEc.clone();
        layout = Layout2;
        a = Bs;
        b = Bc;
        c = D12;
        jniCBLAS.dsymm(layout, side, uplo, m, n, alpha, a, b, beta, c);
        printMatrix("Resulting E", layout, c, m, n);
        System.out.println
        ("\n Test Side = Right: dsymm: G := alpha*A*Bs + beta*G;");
        Side = jniCBLAS.SIDE.Right;
        
        System.out.println("Row-major: ");
        D6 = cpGr.clone();
        layout = Layout;
        side = Side;
        m = M2;
        n = K2;
        a = Bs;
        b = Ar;
        c = D6;
        jniCBLAS.dsymm(layout, side, uplo, m, n, alpha, a, b, beta, c);
        printMatrix("Resulting G", layout, c, m, n);
        
        System.out.println("Col-major: ");
        D6 = cpGc.clone();
        layout = Layout2;
        a = Bs;
        b = Ac;
        c = D6;
        jniCBLAS.dsymm(layout, side, uplo, m, n, alpha, a, b, beta, c);
        printMatrix("Resulting G", layout, c, m, n);
    }
    /** Print the matrix X. */
    private static void printMatrix(String prompt, int layout,
                                    double[] X, int I, int J) {
        System.out.println(prompt);
        if (layout == jniCBLAS.LAYOUT.ColMajor) {
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[j*I+i]));
                System.out.println();
            }
        }
        else if (layout == jniCBLAS.LAYOUT.RowMajor){
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







