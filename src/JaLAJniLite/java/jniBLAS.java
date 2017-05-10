package JaLAJniLite;

/**CBLAS.java*/

public class jniBLAS {
 private jniBLAS() {}
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
    
    public static native void daxpy( int n, double alpha, double[] x, int incx, double[] y, int incy);
    
    public static native double ddot(int n, double[] x, int incx, double[] y, int incy);
    
    /* Level 2: */
    public static native void dgemv(int layout, int trans, int m, int n, double alpha, double[] a, double[] x, int incx, double beta, double[] y, int incy);
    
    public static native void dtrmv(int layout, int uplo, int trans, int diag, int n, double[] a, double[] x, int incx);
    
    public static native void dsymv(int layout, int uplo, int n, double alpha, double[] a, double[] x, int incx, double beta, double[] y, int incy);
    
    /* Level 3: */
    public static native void dgemm(int layout, int transA, int transB, int m, int n, int k, double alpha, double[] a, double[] b, double beta, double[] c);
    
    public static native void dtrmm(int layout, int side, int uplo, int transA, int diag, int m, int n, double alpha, double[] a, double[] b);
    
    public static native void dsymm(int layout, int side, int uplo, int m, int n, double alpha, double[] a, double[] b, double beta, double[] c);
    
    /**inform java virtual machine that function is defined externally*/
    
}







