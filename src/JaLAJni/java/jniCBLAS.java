/**JavaLapack.java*/
package JaLAJni;

public final class jniCBLAS {
    private jniCBLAS(){}
    static {
     
    /* load library (which will contain wrapper for cblas function.)*/
    System.loadLibrary("jni_cblas");
 
    }

    public final static class ORDER {
        private ORDER() {}
        /** row-major arrays */
        public final static int RowMajor= 101;
        /** column-major arrays */
        public final static int ColMajor= 102;
    }
    
    public final static class TRANSPOSE {
        private TRANSPOSE() {}
        /** trans='N' */
        public final static int NoTrans = 111;
        /** trans='T' */
        public final static int Trans= 112;
        /** trans='C' */
        public final static int ConjTrans=113;
    }

    public final static class UPLO {
        private UPLO() {}
        /** Upper triangular matrix */
        public final static int Upper = 121;
        /** Lower triangular matrix*/
        public final static int Lower= 122;
    }
    
    public final static class DIAG {
        private DIAG() {}
        /** not assumed to be unit  */
        public final static int NonUnit = 131;
        /** assumed to be unit */
        public final static int Unit = 132;
    }
    
    public final static class SIDE {
        private SIDE() {}
        /** B := alpha*op( A )*B */
        public final static int Left = 141;
        /** B := alpha*B*op( A ) */
        public final static int Right= 142;
    }
    
    
    


    
/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */
public static native float  sdsdot(int N, float alpha, float[] X, int incX, float[] Y, int incY);

public static native double dsdot(int N, float[] X, int incX, float[] Y, int incY);

public static native float  sdot(int N, float[] X, int incX, float[] Y, int incY);

public static native double ddot(int N, double[] X, int incX, double[] Y, int incY);

/*
 * Functions having prefixes Z and C only
 */

public static native void   cdotu_sub(int N, float[] x, int incX, float[] y, int incY, float[] dotu);

public static native void   cdotc_sub(int N, float[] x, int incX, float[] y, int incY, float[] dotc);

public static native void   zdotu_sub(int N, double[] x, int incX, double[] y, int incY, double[] dotu);

public static native void   zdotc_sub(int N, double[] x, int incX, double[] y, int incY, double[] dotc);


/*
 * Functions having prefixes S D SC DZ
 */

public static native float  snrm2(int N, float[] X, int incX);

public static native float  sasum(int N, float[] X, int incX);

public static native double dnrm2(int N, double[] X, int incX);

public static native double dasum(int N, double[] X, int incX);

public static native float  scnrm2(int N, float[] X, int incX);

public static native float  scasum(int N, float[] X, int incX);

public static native double dznrm2(int N, double[] X, int incX);

public static native double dzasum(int N, double[] X, int incX);


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
//CBLAS_INDEX cblas_isamax(const int N, const float  *X, const int incX);
//CBLAS_INDEX cblas_idamax(const int N, const double *X, const int incX);
//CBLAS_INDEX cblas_icamax(const int N, const void   *X, const int incX);
//CBLAS_INDEX cblas_izamax(const int N, const void   *X, const int incX);

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (s, d, c, z)
 */


public static native void sswap(int N, float[] X, int incX, float[] Y, int incY);

public static native void scopy(int N, float[] X, int incX, float[] Y, int incY);

public static native void saxpy(int N, float alpha, float[] X, int incX, float[] Y, int incY);

public static native void dswap(int N, double[] X, int incX, double[] Y, int incY);

public static native void dcopy(int N, double[] X, int incX, double[] Y, int incY);

public static native void daxpy(int N, double alpha, double[] X, int incX, double[] Y, int incY);

public static native void cswap(int N, float[] x, int incX, float[] y, int incY);

public static native void ccopy(int N, float[] x, int incX, float[] y, int incY);

public static native void caxpy(int N, float[] alpha, float[] x, int incX, float[] y, int incY);

public static native void zswap(int N, double[] x, int incX, double[] y, int incY);

public static native void zcopy(int N, double[] x, int incX, double[] y, int incY);

public static native void zaxpy(int N, double[] alpha, double[] x, int incX, double[] y, int incY);


/* 
 * Routines with S and D prefix only
 */

public static native void srotg(float[] a, float[] b, float[] c, float[] s);

public static native void srotmg(float[] d1, float[] d2, float[] b1, float b2, float[] p);

public static native void srot(int N, float[] x, int incX, float[] y, int incY, float c, float s);

public static native void srotm(int N, float[] x, int incX, float[] y, int incY, float[] p);

public static native void drotg(double[] a, double[] b, double[] c, double[] s);

public static native void drotmg(double[] d1, double[] d2, double[] b1, double b2, double[] p);

public static native void drot(int N, double[] x, int incX, double[] y, int incY, double c, double  s);

public static native void drotm(int N, double[] x, int incX, double[] y, int incY, double[] p);


/* 
 * Routines with S D C Z CS and ZD prefixes
 */

public static native void sscal(int N, float alpha, float[] x, int incX);

public static native void dscal(int N, double alpha, double[] x, int incX);

public static native void cscal(int N, float[] alpha, float[] x, int incX);

public static native void zscal(int N, double[] alpha, double[] x, int incX);

public static native void csscal(int N, float alpha, float[] x, int incX);

public static native void zdscal(int N, double alpha, double[] x, int incX);

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */

public static native void sgemv(int order, int TransA, int M, int N, float alpha, float[] a, int lda, float[] x, int incX, float beta, float[] y, int incY);

public static native void sgbmv(int order, int TransA, int M, int N, int KL, int KU, float alpha, float[] a, int lda, float[] x, int incX, float beta, float[] y, int incY);

public static native void strmv(int order, int Uplo, int TransA, int Diag, int N, float[] a, int lda, float[] x, int incX);

public static native void stbmv(int order, int Uplo, int TransA, int Diag, int N, int K, float[] a, int lda, float[] x, int incX);

public static native void stpmv(int order, int Uplo, int TransA, int Diag, int N, float[] Ap, float[] x, int incX);

public static native void strsv(int order, int Uplo, int TransA, int Diag, int N, float[] A, int lda, float[] X, int incX);

public static native void stbsv(int order, int Uplo, int TransA, int Diag, int N, int K, float[] A, int lda, float[] X, int incX);

public static native void stpsv(int order, int Uplo, int TransA, int Diag, int N, float[] Ap, float[] X, int incX);

public static native void dgemv(int order, int TransA, int M, int N, double alpha, double[] A, int lda, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dgbmv(int order, int TransA, int M, int N, int KL, int KU, double alpha, double[] A, int lda, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dtrmv(int order, int Uplo, int TransA, int Diag, int N, double[] A, int lda, double[] X, int incX);

public static native void dtbmv(int order, int Uplo, int TransA, int Diag, int N, int K, double[] A, int lda, double[] X, int incX);

public static native void dtpmv(int order, int Uplo, int TransA, int Diag, int N, double[] Ap, double[] X, int incX);

public static native void dtrsv(int order, int Uplo, int TransA, int Diag, int N, double[] A, int lda, double[] X, int incX);

public static native void dtbsv(int order, int Uplo, int TransA, int Diag, int N, int K, double[] A, int lda, double[] X, int incX);

public static native void dtpsv(int order, int Uplo, int TransA, int Diag, int N, double[] Ap, double[] X, int incX);

public static native void cgemv(int order, int TransA, int M, int N, float[] alpha, float[] A, int lda, float[] X, int incX, float[] beta, float[] Y, int incY);

public static native void cgbmv(int order, int TransA, int M, int N, int KL, int KU, float[] alpha, float[] A, int lda, float[] X, int incX, float[] beta, float[] Y, int incY);

public static native void ctrmv(int order, int Uplo, int TransA, int Diag, int N, float[] A, int lda, float[] X, int incX);

public static native void ctbmv(int order, int Uplo, int TransA, int Diag, int N, int K, float[] A, int lda, float[] X, int incX);

public static native void ctpmv(int order, int Uplo, int TransA, int Diag, int N, float[] Ap, float[] X, int incX);

public static native void ctrsv(int order, int Uplo, int TransA, int Diag, int N, float[] A, int lda, float[] X, int incX);

public static native void ctbsv(int order, int Uplo, int TransA, int Diag, int N, int K, float[] A, int lda, float[] X, int incX);

public static native void ctpsv(int order, int Uplo, int TransA, int Diag, int N, float[] Ap, float[] X, int incX);

public static native void zgemv(int order, int TransA, int M, int N, double[] alpha, double[] A, int lda,double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void zgbmv(int order, int TransA, int M, int N, int KL, int KU, double[] alpha, double[] A, int lda, double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void ztrmv(int order, int Uplo, int TransA, int Diag, int N, double[] A, int lda, double[] X, int incX);

public static native void ztbmv(int order, int Uplo, int TransA, int Diag, int N, int K, double[] A, int lda, double[] X, int incX);

public static native void ztpmv(int order, int Uplo, int TransA, int Diag, int N, double[] Ap, double[] X, int incX);

public static native void ztrsv(int order, int Uplo, int TransA, int Diag, int N, double[] A, int lda, double[] X, int incX);

public static native void ztbsv(int order, int Uplo, int TransA, int Diag, int N, int K, double[] A, int lda, double[] X, int incX);

public static native void ztpsv(int order, int Uplo, int TransA, int Diag, int N, double[] Ap, double[] X, int incX);


/* 
 * Routines with S and D prefixes only
 */

public static native void ssymv(int order, int Uplo, int N, float alpha, float[] A, int lda, float[] X, int incX, float beta, float[] Y, int incY);

public static native void ssbmv(int order, int Uplo, int N, int K, float alpha, float[] A, int lda, float[] X, int incX, float beta, float[] Y, int incY);

public static native void sspmv(int order, int Uplo, int N, float alpha, float[] Ap, float[] X, int incX, float beta, float[] Y, int incY);

public static native void sger(int order, int M, int N, float alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void ssyr(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] A, int lda);

public static native void sspr(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] Ap);

public static native void ssyr2(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void sspr2(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] Y, int incY, float[] A);

public static native void dsymv(int order, int Uplo, int N, double alpha, double[] A, int lda, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dsbmv(int order, int Uplo, int N, int K, double alpha, double[] A, int lda, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dspmv(int order, int Uplo, int N, double alpha, double[] Ap, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dger(int order, int M, int N, double alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void dsyr(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] A, int lda);

public static native void dspr(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] Ap);

public static native void dsyr2(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void dspr2(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] Y, int incY, double[] A);


/* 
 * Routines with C and Z prefixes only
 */

public static native void chemv(int order, int Uplo, int N, float[] alpha, float[] A, int lda, float[] X, int incX,float[] beta, float[] Y, int incY);

public static native void chbmv(int order, int Uplo, int N, int K, float[] alpha, float[] A, int lda, float[] X, int incX,float[] beta, float[] Y, int incY);

public static native void chpmv(int order, int Uplo, int N, float[] alpha, float[] Ap, float[] X, int incX,float[] beta, float[] Y, int incY);

public static native void cgeru(int order, int M, int N, float[] alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void cgerc(int order, int M, int N, float[] alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void cher(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] A, int lda);

public static native void chpr(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] A);

public static native void cher2(int order, int Uplo, int N, float[] alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void chpr2(int order, int Uplo, int N, float[] alpha, float[] X, int incX, float[] Y, int incY, float[] Ap);

public static native void zhemv(int order, int Uplo, int N, double[] alpha, double[] A, int lda, double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void zhbmv(int order, int Uplo, int N, int K, double[] alpha, double[] A, int lda, double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void zhpmv(int order, int Uplo, int N, double[] alpha, double[] Ap, double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void zgeru(int order, int M, int N, double[] alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void zgerc(int order, int M, int N, double[] alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void zher(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] A, int lda);

public static native void zhpr(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] A);

public static native void zher2(int order, int Uplo, int N, double[] alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void zhpr2(int order, int Uplo, int N, double[] alpha, double[] X, int incX, double[] Y, int incY, double[] Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */

public static native void sgemm(int Order, int TransA, int TransB, int M, int N, int K, float alpha, float[] A, int lda, float[] B, int ldb, float beta, float[] C, int ldc);

public static native void ssymm(int Order, int Side, int Uplo, int M, int N, float alpha, float[] A, int lda, float[] B, int ldb, float beta, float[] C, int ldc);

public static native void ssyrk(int Order, int Uplo, int Trans, int N, int K, float alpha, float[] A, int lda, float beta, float[] C, int ldc);

public static native void ssyr2k(int Order, int Uplo, int Trans, int N, int K, float alpha, float[] A, int lda, float[] B, int ldb,  float beta, float[] C, int ldc);

public static native void strmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float alpha, float[] A, int lda, float[] B, int ldb);

public static native void strsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float alpha, float[] A, int lda, float[] B, int ldb);

public static native void dgemm(int Order, int TransA, int TransB, int M, int N, int K, double alpha, double[] A, int lda, double[] B, int ldb, double beta, double[] C, int ldc);

public static native void dsymm(int Order, int Side, int Uplo, int M, int N, double alpha, double[] A, int lda, double[] B, int ldb, double beta, double[] C, int ldc);

public static native void dsyrk(int Order, int Uplo, int Trans, int N, int K, double alpha, double[] A, int lda, double beta, double[] C, int ldc);

public static native void dsyr2k(int Order, int Uplo, int Trans, int N, int K, double alpha, double[] A, int lda, double[] B, int ldb, double beta, double[] C, int ldc);

public static native void dtrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double alpha, double[] A, int lda, double[] B, int ldb);

public static native void dtrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double alpha, double[] A, int lda, double[] B, int ldb);

public static native void cgemm(int Order, int TransA, int TransB, int M, int N, int K, float[] alpha, float[] A, int lda, float[] B, int ldb, float[] beta, float[] C, int ldc);

public static native void csymm(int Order, int Side, int Uplo, int M, int N, float[] alpha, float[] A, int lda, float[] B, int ldb, float[] beta, float[] C, int ldc);

public static native void csyrk(int Order, int Uplo, int Trans, int N, int K, float[] alpha, float[] A, int lda, float[] beta, float[] C, int ldc);

public static native void csyr2k(int Order, int Uplo, int Trans, int N, int K, float[] alpha, float[] A, int lda, float[] B, int ldb, float[] beta, float[] C, int ldc);

public static native void ctrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float[] alpha, float[] A, int lda, float[] B, int ldb);

public static native void ctrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float[] alpha, float[] A, int lda, float[] B, int ldb);

public static native void zgemm(int Order, int TransA, int TransB, int M, int N, int K, double[] alpha, double[] A, int lda, double[] B, int ldb, double[] beta, double[] C, int ldc);

public static native void zsymm(int Order, int Side, int Uplo, int M, int N, double[] alpha, double[] A, int lda, double[] B, int ldb, double[] beta, double[] C, int ldc);

public static native void zsyrk(int Order, int Uplo, int Trans, int N, int K, double[] alpha, double[] A, int lda, double[] beta, double[] C, int ldc);

public static native void zsyr2k(int Order, int Uplo, int Trans, int N, int K, double[] alpha, double[] A, int lda, double[] B, int ldb, double[] beta, double[] C, int ldc);

public static native void ztrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double[] alpha, double[] A, int lda, double[] B, int ldb);

public static native void ztrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double[] alpha, double[] A, int lda, double[] B, int ldb);


/* 
 * Routines with prefixes C and Z only
 */

public static native void chemm(int Order, int Side, int Uplo, int M, int N, float[] alpha, float[] A, int lda, float[] B, int ldb, float[] beta, float[] C, int ldc);

public static native void cherk(int Order, int Uplo, int Trans, int N, int K, float alpha, float[] A, int lda, float beta, float[] C, int ldc);

public static native void cher2k(int Order, int Uplo, int Trans, int N, int K, float[] alpha, float[] A, int lda, float[] B, int ldb, float beta, float[] C, int ldc);

public static native void zhemm(int Order, int Side, int Uplo, int M, int N, double[] alpha, double[] A, int lda, double[] B, int ldb, double[] beta, double[] C, int ldc);

public static native void zherk(int Order, int Uplo, int Trans, int N, int K, double alpha, double[] A, int lda, double beta, double[] C, int ldc);

public static native void zher2k(int Order, int Uplo, int Trans, int N, int K, double[] alpha, double[] A, int lda, double[] B, int ldb, double beta, double[] C, int ldc);

//public static native void xerbla(int p, char[] rout, char[] form, ...);

    
}

