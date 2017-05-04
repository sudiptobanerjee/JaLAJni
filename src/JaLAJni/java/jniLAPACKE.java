/**jniLAPACKE.java*/
package JaLAJni;

public final class jniLAPACKE {
 private jniLAPACKE() {}
 static {
     
    /* load library (which will contain wrapper for lapacke.)*/
    System.loadLibrary("jni_java_lapacke");
 
 }
    
    public final static class LAYOUT {
        private LAYOUT() {}
        public final static int RowMajor = 101;
        public final static int ColMajor = 102;
    }
    
    public final static class TRANSPOSE {
        private TRANSPOSE() {}
        public final static char NoTrans = 'N';         /** trans='N' */
        public final static char Trans= 'T';            /** trans='T' */
        public final static char ConjTrans= 'C';        /** trans='C' */
    }
    
    public final static class UPLO {
        private UPLO() {}
        public final static char Upper = 'U';           /** Upper triangular matrix */
        public final static char Lower= 'L';            /** Lower triangular matrix*/
    }
    
    public final static class JOBV {
        private JOBV() {}
        public final static char NoCompute = 'N';       /** eigenvectors are not computed */
        public final static char Compute= 'V';          /** eigenvectors are computed*/
    }
    
    public final static class JOB {
        private JOB() {}
        public final static char All = 'A';             /** all M columns of U are returned in array U */
        public final static char firstInU = 'S';        /** the first min(m,n) columns of U (the left singular
                                                         vectors) are returned in the array U;*/
        public final static char Overwritten = 'O';     /** the first min(m,n) columns of U (the left singular
                                                         vectors) are overwritten on the array A; */
        public final static char NoCompute = 'N';       /** no columns of U (no left singular vectors) are
                                                         computed.*/
    }
    
    public final static class ITYPE {
	private ITYPE(){}
	public final static int first = 1;
	public final static int second = 2;
	public final static int third = 3;
    }

    /* LU */
    public static native int dgetrf(int matrix_layout, int m, int n, double[] a, int lda, int[] ipiv);
    
    public static native int dgetrs(int matrix_layout, char trans, int n, int nrhs, double[] a, int lda, int[] ipiv, double[] b, int ldb);
    
    public static native int dgetri(int matrix_layout, int n, double[] a, int lda, int[] ipiv);
    
    /* Cholesky */
    public static native int dpotrf(int matrix_layout, char uplo, int n, double[] a, int lda);
    
    public static native int dpotri(int matrix_layout, char uplo, int n, double[] a, int lda);
    
    public static native int dpotrs(int matrix_layout, char uplo, int n, int nrhs, double[] a, int lda, double[] b, int ldb);
    
    /* QR */
    public static native int dgeqrf(int matrix_layout, int m, int n, double[] a, int lda, double[] tau);
    
    public static native int dorgqr(int matrix_layout, int m, int n, int k, double[] a, int lda, double[] tau);
    
    public static native int dgeqp3(int matrix_layout, int m, int n, double[] a, int lda, int[] jpvt, double[] tau);
    
    /* Eigenvector and SVD */
    public static native int dgeev(int matrix_layout, char jobvl, char jobvr, int n, double[] a, int lda, double[] wr, double[] wi, double[] vl, int ldvl, double[] vr, int ldvr);
    
    public static native int dgesvd(int matrix_layout, char jobu, char jobvt, int m, int n, double[] a, int lda, double[] s, double[] u, int ldu, double[] vt, int ldvt, double[] superb);
    
    public static native int dgesdd(int matrix_layout, char jobz, int m, int n, double[] a, int lda, double[] s, double[] u, int ldu, double[] vt, int ldvt);
    
    /**inform java virtual machine that function is defined externally*/

    public static native void dsyev(int matrix_layout, char jobz, char uplo, int n, double[] a, int lda, double[] w );
 

    public static native int dgesv( int matrix_layout, int n, int nrhs, double[] a, int lda, int[] ipiv, double[] b, int ldb );

    public static native int dggev( int matrix_layout, char jobvl, char jobvr, int n, double[] a, int lda, double[] b, int ldb, double[] alphar, double[] alphai, double[] beta, double[] vl, int ldvl, double[] vr, int ldvr );

    public static native int dsygv( int matrix_layout, int itype, char jobz, char uplo, int n, double[] a, int lda, double[] b, int ldb, double[] w );

   
}







