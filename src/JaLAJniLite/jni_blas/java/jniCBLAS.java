package JaLAJniLite.jni_blas;

/**CBLAS.java*/

public class jniCBLAS {
 private jniCBLAS() {}
 static {
     
    /* load library (which will contain wrapper for cblas function.)*/
    System.loadLibrary("blas_lite");
 
 }
    public static class LAYOUT {
        private LAYOUT() {}
        /** row-major arrays */
        public final static int RowMajor= 101;
        /** column-major arrays */
        public final static int ColMajor= 102;
    }
    
    public static class TRANSPOSE {
        private TRANSPOSE() {}
        /** trans = 'N' */
        public final static int NoTrans = 111;
        /** trans = 'T' */
        public final static int Trans= 112;
        /** trans = 'C'*/
        public final static int ConjTrans= 113;
    }

    public static class UPLO {
        private UPLO() {}
        /** Upper triangular matrix */
        public final static int Upper = 121;
        /** Lower triangular matrix*/
        public final static int Lower = 122;
    }
    
    public static class DIAG {
        private DIAG() {}
        /** not assumed to be unit  */
        public final static int NonUnit = 131;
        /** assumed to be unit */
        public final static int Unit = 132;
    }
    
    public static class SIDE {
        private SIDE() {}
        /** B := alpha*op( A )*B */
        public final static int Left = 141;
        /** B := alpha*B*op( A ) */
        public final static int Right = 142;
    }
    
    /* Level 1: */
    public static native void dscal( int n, double alpha, double[] x, int incx);
    
    public static native void daxpy( int n, double alpha, double[] x, int incx,
                                    double[] y, int incy);
    
    public static native double ddot(int n, double[] x, int incx, double[] y,
                                     int incy);
    
    /* Level 2: */
    public static native void dgemv(int Layout, int Trans, int m, int n, double
                                    alpha, double[] A, double[] x, int incx,
                                    double beta, double[] y, int incy);
    
    public static native void dtrmv(int Layout, int Uplo, int Trans, int Diag,
                                    int n, double[] A, double[] x, int incx);
    
    public static native void dsymv(int Layout, int Uplo, int n, double alpha,
                                    double[] A, double[] x, int incx, double beta,
                                    double[] y, int incy);
    
    /* Level 3: */
    public static native void dgemm(int Layout, int TransA, int TransB, int m,
                                    int n, int k, double alpha, double[] A,
                                    double[] B, double beta, double[] C);
    
    public static native void dtrmm(int Layout, int Side, int Uplo, int TransA,
                                    int Diag, int m, int n, double alpha,
                                    double[] A, double[] B);
    
    public static native void dsymm(int Layout, int Side, int Uplo, int m, int n,
                                    double alpha, double[] A, double[] B,
                                    double beta, double[] C);
    
    /**inform java virtual machine that function is defined externally*/
    
}







