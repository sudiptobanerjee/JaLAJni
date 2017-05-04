//
//  main.c
//  jniLAPACKE
//
//  Created by Lu Zhang on 5/8/16.
//  Copyright © 2016 Lu Zhang. All rights reserved.
//
/* jniLAPACKE.c */
#include <jni.h>
#include <assert.h>
#include <lapacke.h>


/* LU */

JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dgetrf (JNIEnv *env, jclass klass, jint matrix_layout, jint m, jint n, jdoubleArray a, jint lda, jintArray ipiv){
    
    double *aElems;
    lapack_int *ipivElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    ipivElems = (*env)-> GetIntArrayElements (env, ipiv, NULL);
    
    assert(aElems && ipivElems);
    
    info = LAPACKE_dgetrf ((int) matrix_layout, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, ipivElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseIntArrayElements (env, ipiv, ipivElems, 0);
    
    return info;
}

JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dgetrs (JNIEnv *env, jclass klass, jint matrix_layout, jchar trans, jint n, jint nrhs, jdoubleArray a, jint lda, jintArray ipiv, jdoubleArray b, jint ldb){
    
    double *bElems;
    double *aElems;
    lapack_int *ipivElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    bElems = (*env)-> GetDoubleArrayElements (env, b, NULL);
    ipivElems = (*env)-> GetIntArrayElements (env, ipiv, NULL);
    
    assert(aElems && bElems && ipivElems);
    
    info = LAPACKE_dgetrs ((int) matrix_layout, (char) trans, (lapack_int) n, (lapack_int) nrhs, aElems, (lapack_int) lda, ipivElems, bElems, (lapack_int) ldb);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, b, bElems, 0);
    (*env)-> ReleaseIntArrayElements (env, ipiv, ipivElems, JNI_ABORT);
    
    return info;
}


JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dgetri (JNIEnv *env, jclass klass, jint matrix_layout, jint n, jdoubleArray a, jint lda, jintArray ipiv){
    
    double *aElems;
    lapack_int *ipivElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    ipivElems = (*env)-> GetIntArrayElements (env, ipiv, NULL);
    
    assert(aElems && ipivElems);
    
    info = LAPACKE_dgetri((int) matrix_layout, (lapack_int) n, aElems, (lapack_int) lda, ipivElems );
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseIntArrayElements (env, ipiv, ipivElems, JNI_ABORT);
    
    return info;
}

/* Cholesky */

JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dpotrf (JNIEnv *env, jclass klass, jint matrix_layout, jchar uplo, jint n, jdoubleArray a, jint lda){
    
    double *aElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);

    assert(aElems);
    
    info = LAPACKE_dpotrf((int) matrix_layout, (char) uplo, (lapack_int) n, aElems, (lapack_int) lda);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    
    return info;

}


JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dpotri(JNIEnv *env, jclass klass, jint matrix_layout, jchar uplo, jint n, jdoubleArray a, jint lda){
 
    double *aElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    
    assert(aElems);
    
    info = LAPACKE_dpotri((int) matrix_layout, (char) uplo, (lapack_int) n, aElems, (lapack_int) lda);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    
}

JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dpotrs (JNIEnv *env, jclass klass, jint matrix_layout, jchar uplo, jint n, jint nrhs, jdoubleArray a, jint lda, jdoubleArray b, jint ldb){
    
    double *aElems, *bElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    bElems = (*env)-> GetDoubleArrayElements (env, b, NULL);
    
    assert(aElems && bElems);
    
    info = LAPACKE_dpotrs ((int) matrix_layout, (char) uplo, (lapack_int) n, (lapack_int) nrhs, aElems, (lapack_int) lda, bElems, (lapack_int) ldb);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, b, bElems, 0);
    
    return info;
    
}

/* QR */

JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dgeqrf (JNIEnv *env, jclass klass, jint matrix_layout, jint m, jint n, jdoubleArray a, jint lda, jdoubleArray tau){
    
    double *aElems, *tauElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    tauElems = (*env)-> GetDoubleArrayElements (env, tau, NULL);
    
    assert(aElems && tauElems);
    
    info = LAPACKE_dgeqrf ((int) matrix_layout, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, tauElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, tau, tauElems, 0);
    
    return info;

}

JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dorgqr (JNIEnv *env, jclass klass, jint matrix_layout, jint m, jint n, jint k, jdoubleArray a, jint lda, jdoubleArray tau){
    
    double *aElems, *tauElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    tauElems = (*env)-> GetDoubleArrayElements (env, tau, NULL);
    
    assert(aElems && tauElems);
    
    info = LAPACKE_dorgqr((int) matrix_layout, (lapack_int) m, (lapack_int) n, (lapack_int) k, aElems, (lapack_int) lda, tauElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, tau, tauElems, 0);
    
    return info;
    
}


JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dgeqp3  (JNIEnv *env, jclass klass, jint matrix_layout, jint m, jint n, jdoubleArray a, jint lda, jintArray jpvt, jdoubleArray tau){
    
    double *aElems, *tauElems;
    lapack_int *jpvtElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    tauElems = (*env)-> GetDoubleArrayElements (env, tau, NULL);
    jpvtElems = (*env)-> GetIntArrayElements (env, jpvt, NULL);
    
    assert(aElems && tauElems && jpvtElems);
    
    info = LAPACKE_dgeqp3 ((int) matrix_layout, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, jpvtElems, tauElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, tau, tauElems, 0);
    (*env)-> ReleaseIntArrayElements (env, jpvt, jpvtElems, 0);
    
}

/*Eigenvector and SVD*/

JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dgeev (JNIEnv *env, jclass klass, jint matrix_layout, jchar jobvl, jchar jobvr, jint n, jdoubleArray a, jint lda, jdoubleArray wr, jdoubleArray wi, jdoubleArray vl, jint ldvl, jdoubleArray vr, jint ldvr){
    
    double *aElems, *wrElems, *wiElems, *vlElems, *vrElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    wrElems = (*env)-> GetDoubleArrayElements (env, wr, NULL);
    wiElems = (*env)-> GetDoubleArrayElements (env, wi, NULL);
    vlElems = (*env)-> GetDoubleArrayElements (env, vl, NULL);
    vrElems = (*env)-> GetDoubleArrayElements (env, vr, NULL);
    
    assert(aElems && wrElems && wiElems && vlElems && vrElems);
    
    info = LAPACKE_dgeev((int) matrix_layout, (char) jobvl, (char) jobvr, (lapack_int) n, aElems, (lapack_int) lda, wrElems, wiElems, vlElems, (lapack_int) ldvl, vrElems, (lapack_int) ldvr);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, vl, vlElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, vr, vrElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, wr, wrElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, wi, wiElems, 0);
    
    return info;

}


JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dgesvd (JNIEnv *env, jclass klass, jint matrix_layout, jchar jobu, jchar jobvt, jint m, jint n, jdoubleArray a, jint lda, jdoubleArray s, jdoubleArray u, jint ldu, jdoubleArray vt, jint ldvt, jdoubleArray superb){
    
    //superb contains the unconverged superdiagonal elements of an upper bidiagonal matrix B whose diagonal is in S (not necessarily sorted).
    
    double *aElems, *sElems, *uElems, *vtElems, *superbElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    sElems = (*env)-> GetDoubleArrayElements (env, s, NULL);
    uElems = (*env)-> GetDoubleArrayElements (env, u, NULL);
    vtElems = (*env)-> GetDoubleArrayElements (env, vt, NULL);
    superbElems = (*env)-> GetDoubleArrayElements (env, superb, NULL);
    
    assert(aElems && sElems && uElems && vtElems && superbElems);
    
    info = LAPACKE_dgesvd((int) matrix_layout, (char) jobu, (char) jobvt, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, sElems, uElems, (lapack_int) ldu, vtElems, (lapack_int) ldvt, superbElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, s, sElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, u, uElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, vt, vtElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, superb, superbElems, 0);
    
    return info;
}

JNIEXPORT jint Java_JaLAJni_jniLAPACKE_dgesdd (JNIEnv *env, jclass klass, jint matrix_layout, jchar jobz, jint m, jint n, jdoubleArray a, jint lda, jdoubleArray s, jdoubleArray u, jint ldu, jdoubleArray vt, jint ldvt){
    
    
    double *aElems, *sElems, *uElems, *vtElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    sElems = (*env)-> GetDoubleArrayElements (env, s, NULL);
    uElems = (*env)-> GetDoubleArrayElements (env, u, NULL);
    vtElems = (*env)-> GetDoubleArrayElements (env, vt, NULL);
    
    assert(aElems && sElems && uElems && vtElems);
    
    info = LAPACKE_dgesdd((int) matrix_layout, (char) jobz, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, sElems, uElems, (lapack_int) ldu, vtElems, (lapack_int) ldvt);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, s, sElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, u, uElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, vt, vtElems, 0);
    
    return info;
}


JNIEXPORT void JNICALL Java_JaLAJni_jniLAPACKE_dsyev
  (JNIEnv *env, jclass obj, jint layout, jchar jobz, jchar uplo, jint n, jdoubleArray ja, jint lda, jdoubleArray jw)
{
    jdouble *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    jdouble *w = (*env)-> GetDoubleArrayElements (env, jw, NULL);
    
    assert(a && w);
    
    LAPACKE_dsyev((int) layout, (char) jobz, (char) uplo, (lapack_int) n, a, (lapack_int) lda, w);
    
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jw, w, 0);
}



JNIEXPORT jint JNICALL Java_JaLAJni_jniLAPACKE_dgesv
  (JNIEnv *env, jclass obj, jint layout, jint n, jint nrhs, jdoubleArray ja, jint lda, jintArray jipiv, jdoubleArray jb, jint ldb)
{
    int result;    

    jdouble *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    jint *ipiv = (*env)-> GetIntArrayElements (env, jipiv, NULL);
    jdouble *b = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    
    assert(a && ipiv && b);
    
    result = LAPACKE_dgesv((int) layout, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, ipiv, b, (lapack_int) ldb);
    
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseIntArrayElements (env, jipiv, ipiv, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jb, b, 0);

    return result;
}



JNIEXPORT jint JNICALL Java_JaLAJni_jniLAPACKE_dggev
  (JNIEnv *env, jclass obj, jint layout, jchar jobvl, jchar jobvr, jint n, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdoubleArray jalphar, jdoubleArray jalphai, jdoubleArray jbeta, jdoubleArray jvl, jint ldvl, jdoubleArray jvr, jint ldvr)
{
    int result;

    jdouble *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    jdouble *b = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    jdouble *alphar = (*env)-> GetDoubleArrayElements (env, jalphar, NULL);
    jdouble *alphai = (*env)-> GetDoubleArrayElements (env, jalphai, NULL);
    jdouble *beta = (*env)-> GetDoubleArrayElements (env, jbeta, NULL);
    jdouble *vl = (*env)-> GetDoubleArrayElements (env, jvl, NULL);
    jdouble *vr = (*env)-> GetDoubleArrayElements (env, jvr, NULL);
    
    assert(a && b && alphar && alphai && beta && vl && vr);
    
    result = LAPACKE_dggev((int) layout, (char) jobvl, (char)jobvr, (lapack_int) n, a, (lapack_int)lda, b, (lapack_int)ldb, alphar, alphai, beta, vl ,(lapack_int)ldvl, vr, (lapack_int)ldvr);
    
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jb, b, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jalphar, alphar, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jalphai, alphai, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jbeta, beta, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jvl, vl, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jvr, vr, 0);

    return result;
}



JNIEXPORT jint JNICALL Java_JaLAJni_jniLAPACKE_dsygv
  (JNIEnv *env, jclass obj, jint layout, jint itype, jchar jobz, jchar uplo, jint n, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdoubleArray jw)
{
    int result;

    jdouble *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    jdouble *b = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    jdouble *w = (*env)-> GetDoubleArrayElements (env, jw, NULL);
    
    assert(a && b && w);
    
    result = LAPACKE_dsygv((int) layout, (lapack_int) itype, (char) jobz, (char) uplo, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, w);
    
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jb, b, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jw, w, 0);

    return result;
}






