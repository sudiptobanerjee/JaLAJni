import JaLAJni.*;

public final class cblastest {
    /** Incarnation prohibited. */
    /** No command-line options. */
    private cblastest(){}
    public static void main(String[] args) {
	
	
        //
        // Prepare the matrices and other parameters
        //
        int M=2, N=4, K=3;
        float alpha = 2, beta = -1;
	double alphad = 2, betad = -1;
	float[] arrayconst = new float[] {2,-1};	
	double[] arrayconstd = new double[] {2,-1};	
	
	int order = jniCBLAS.ORDER.ColMajor; 	
	int transA = jniCBLAS.TRANSPOSE.NoTrans;
	int transB = jniCBLAS.TRANSPOSE.NoTrans;
	int uplo;
	int diag;
	int side;

        float[] A = new float[] {1,4,2,5,3,6};//with dim M*K
        float[] cpA = new float[] {1,4,2,5,3,6};
	double[] dA = new double[] {1,4,2,5,3,6};
	double[] cpdA = new double[] {1,4,2,5,3,6};
	float[] Ac = new float[]{1,1,4,1,2,1,5,2,3,1,6,1};//dim M*K
	float[] cpAc = new float[]{1,1,4,1,2,1,5,2,3,1,6,1};//dim M*K
	double[] Az = new double[]{1,1,4,1,2,1,5,2,3,1,6,1};//dim M*K
	double[] cpAz = new double[]{1,1,4,1,2,1,5,2,3,1,6,1};//dim M*K
	
	
	float resultf;
	double resultd;
	float[] arrayresult1 = new float[] {0};
	float[] arrayresult2 = new float[] {0};
	double[] arrayresult1d = new double[] {0};
	double[] arrayresult2d = new double[] {0};
	
	float[] dotf = new float[] {1,0,1,-1};//with dim M
	float[] cpdotf = new float[] {1,3,1,0};//with dim M
	float[] dotuf = new float[]{0,0};
	double[] dotd = new double[] {1,0,1,-1};//with dim M
	double[] cpdotd = new double[] {1,3,1,0};//with dim M
	double[] dotud = new double[]{0,1};
	
	float[] dotrf = new float[] {1,0};
	float[] dotif = new float[] {1,-1};
	float[] cpdotrf = new float[] {1,0};
	float[] cpdotif = new float[] {1,-1};
	double[] dotr = new double[] {1,0};
	double[] doti = new double[] {1,-1};
	double[] cpdotr = new double[] {1,0};
	double[] cpdoti = new double[] {1,-1};

	float[] af = new float[]{8};
	float[] bf = new float[] {6};
	double[] ad = new double[]{8};
	double[] bd = new double[]{6};
	float[] tmpx = new float[]{1,2,3,4,5};
	float[] tmpy = new float[]{-1,0,-2,0,-3,0,-4,0,-5};
	double[] tmpxd = new double[]{1,2,3,4,5};
	double[] tmpyd = new double[]{-1,0,-2,0,-3,0,-4,0,-5};
	double[] bandfull = new double[]{1,4,0,6,2,5,0,7,3};//dim K*K
	double[] bandd = new double[] {0,1,6,4,2,7,5,3,0};//dim K*K
	float[] bandf = new float[] {0,1,6,4,2,7,5,3,0};//dim K*K
	
	float[] xf = new float[] {1,2,3};
	float[] yf = new float[] {0,1,-1};
	float[] cpxf = new float[] {1,2,3};
	float[] cpyf = new float[] {0,1,-1};
	
	double[] xd = new double[] {1,2,3};
	double[] yd = new double[] {0,1,-1};
	double[] cpxd = new double[] {1,2,3};
	double[] cpyd = new double[] {0,1,-1};

	float[] trif = new float[] {1,0,0,2,4,0,3,5,0};//dim K*K
	double[] trid = new double[] {1,0,0,2,4,0,3,5,0,};
	float[] tribandfull = new float[] {1,0,0,2,4,0,0,5,3};
	double[] tribandfulld = new double[] {1,0,0,2,4,0,0,5,3};
	float[] tribandf = new float[] {0,1,0,2,4,0,5,3,0};
	float[] packed = new float[] {1,2,4,0,5,3};
	double[] tribandd = new double[] {0,1,0,2,4,0,5,3,0};
	double[] packedd = new double[] {1,2,4,0,5,3};

	float[] ac = new float[]{1,1,1,-1,2,1,2,-1};
	float[] xc = new float[]{1,-1,1,1};
	float[] yc = new float[]{3,2,2,3};
	float[] cpxc = new float[]{1,-1,1,1};
	float[] cpyc = new float[]{3,2,2,3};
	float[] alphac = new float[]{0,1};
	float[] betac = new float[]{0,-1};
	
	double[] az = new double[]{1,1,1,-1,2,1,2,-1};
	double[] xz = new double[]{1,-1,1,1};
	double[] yz = new double[]{3,2,2,3};
	double[] cpxz = new double[]{1,-1,1,1};
	double[] cpyz = new double[]{3,2,2,3};
	double[] alphaz = new double[]{0,1};
	double[] betaz = new double[]{0,-1};

	float[] bandfullc = new float[]{1,1,2,0,0,0,1,-1,1,1,3,0,0,0,1,-1,1,1};//row-major
	float[] bandfullccol = new float[]{1,1,1,-1,0,0,2,0,1,1,1,-1,0,0,3,0,1,1};
	float[] bandccol = new float[]{0,0,1,1,1,-1,2,0,1,1,1,-1,3,0,1,1,0,0};
	float[] bandcrow = new float[]{0,0,1,1,2,0,1,-1,1,1,3,0,1,-1,1,1,0,0};
	// this is something need to be modified
	float[] bxc = new float[] {1,-1,1,-1,0,0};
	float[] byc = new float[] {2,-1,0,0,2,1};
	float[] cpbxc = new float[] {1,-1,1,-1,0,0};
	float[] cpbyc = new float[] {2,-1,0,0,2,1};
	float[] tric = new float[] {1,1,2,0,0,0,0,0,1,1,3,0,0,0,0,0,1,1};
	float[] tribandc = new float[] {0,0,1,1,0,0,2,0,1,1,0,0,3,0,1,1,0,0};
	float[] packedc = new float[] {1,1,2,0,1,1,0,0,3,0,1,1};
	
	double[] bandfullz = new double[]{1,1,2,0,0,0,1,-1,1,1,3,0,0,0,1,-1,1,1};
	double[] bandzcol = new double[]{0,0,1,1,1,-1,2,0,1,1,1,-1,3,0,1,1,0,0};
	double[] bandzrow = new double[]{0,0,1,1,2,0,1,-1,1,1,3,0,1,-1,1,1,0,0};
	double[] bxz = new double[] {1,-1,1,-1,0,0};
	double[] byz = new double[] {2,-1,0,0,2,1};
	double[] cpbxz = new double[] {1,-1,1,-1,0,0};
	double[] cpbyz = new double[] {2,-1,0,0,2,1};
	double[] triz = new double[] {1,1,2,0,0,0,0,0,1,1,3,0,0,0,0,0,1,1};
	double[] trizfull = new double[]{1,1,0,0,0,0,2,0,1,1,0,0,0,0,3,0,1,1};
	double[] tribandz = new double[] {0,0,1,1,0,0,2,0,1,1,0,0,3,0,1,1,0,0};
	double[] packedz = new double[] {1,1,2,0,1,1,0,0,3,0,1,1};

	float[] symfulls = new float[]{1,4,0,4,2,5,0,5,3};
	float[] cpsymfulls = new float[]{1,4,0,4,2,5,0,5,3};
	float[] symbands = new float[]{0,1,4,4,2,5,5,3,0};
	float[] packedsyms = new float[]{1,4,2,0,5,3};
	float[] cppackedsyms = new float[]{1,4,2,0,5,3};
	
	double[] symfulld = new double[]{1,4,0,4,2,5,0,5,3};
	double[] cpsymfulld = new double[]{1,4,0,4,2,5,0,5,3};
	double[] symbandd = new double[]{0,1,4,4,2,5,5,3,0};
	double[] packedsymd = new double[]{1,4,2,0,5,3};
	double[] cppackedsymd = new double[]{1,4,2,0,5,3};

	float[] hermc = new float[]{1,0,1,2,0,0,1,-2,1,0,3,1,0,0,3,-1,1,0};//dim K*K
	float[] cphermc = new float[]{1,0,1,2,0,0,1,-2,1,0,3,1,0,0,3,-1,1,0};
	float[] hermbandc = new float[]{0,0,1,0,1,-2,1,2,1,0,3,-1,3,1,1,0,0,0};
	float[] packedhermc = new float[]{1,0,1,2,1,0,0,0,3,1,1,0};
	float[] cppackedhermc = new float[]{1,0,1,2,1,0,0,0,3,1,1,0};	
	
	double[] hermz = new double[]{1,0,1,2,0,0,1,-2,1,0,3,1,0,0,3,-1,1,0};
	double[] cphermz = new double[]{1,0,1,2,0,0,1,-2,1,0,3,1,0,0,3,-1,1,0};
	double[] hermbandz = new double[]{0,0,1,0,1,-2,1,2,1,0,3,-1,3,1,1,0,0,0};
	double[] packedhermz = new double[]{1,0,1,2,1,0,0,0,3,1,1,0};
	double[] cppackedhermz = new double[]{1,0,1,2,1,0,0,0,3,1,1,0};
	
	float[] Bf = new float[] {7,11,15,8,12,16,9,13,17,10,14,18};//with dim K*N
        float[] Cf = new float[] {2,-3,-4,-7,3,0,7,18};//with dim M*N
        float[] cpCf = new float[] {2,-3,-4,-7,3,0,7,18};
	float[] myaf = new float[] {1,2,3,4,5,6};//dim M*K
	float[] mybf = new float[] {-1,0,1,1,0,-1};//dim M*K
	float[] mycf = new float[] {1,0,0,1};//dim M*M        
	float[] cpmyaf = new float[] {1,2,3,4,5,6};//dim M*K
	float[] cpmybf = new float[] {-1,0,1,1,0,-1};//dim M*K
	float[] cpmycf = new float[] {1,0,0,1};//dim M*M  
	float[] solvetrif = new float[] {1,0,0,2,4,0,3,5,1};      
	
	double[] Bd = new double[] {7,11,15,8,12,16,9,13,17,10,14,18};//with dim K*N
        double[] Cd = new double[] {2,-3,-4,-7,3,0,7,18};//with dim M*N
        double[] cpCd = new double[] {2,-3,-4,-7,3,0,7,18};
	double[] myad = new double[] {1,2,3,4,5,6};//dim M*K
	double[] mybd = new double[] {-1,0,1,1,0,-1};//dim M*K
	double[] mycd = new double[] {1,0,0,1};//dim M*M        
	double[] cpmyad = new double[] {1,2,3,4,5,6};//dim M*K
	double[] cpmybd = new double[] {-1,0,1,1,0,-1};//dim M*K
	double[] cpmycd = new double[] {1,0,0,1};//dim M*M  
	double[] solvetrid = new double[] {1,0,0,2,4,0,3,5,1};   

	float[] Bc = new float[] {7,1,11,1,15,1,8,0,12,1,16,2,9,0,13,1,17,1,10,0,14,1,18,0};//with dim K*N
        float[] Cc = new float[] {2,0,-3,1,-4,0,-7,1,3,0,0,1,7,0,18,0};//with dim M*N
        float[] cpCc = new float[] {2,0,-3,1,-4,0,-7,1,3,0,0,1,7,0,18,0};
	float[] myac = new float[] {1,0,2,-1,3,0,4,1,5,1,6,0};//dim M*K
	float[] mybc = new float[] {-1,1,0,1,1,0,1,0,0,0,-1,2};//dim M*K
	float[] mycc = new float[] {1,1,0,0,0,0,1,-1};//dim M*M        
	float[] cpmyac = new float[] {1,0,2,-1,3,0,4,1,5,1,6,0};//dim M*K
	float[] cpmybc = new float[] {-1,1,0,1,1,0,1,0,0,0,-1,2};//dim M*K
	float[] cpmycc = new float[] {1,1,0,0,0,0,1,-1};//dim M*M  
	float[] solvetric = new float[] {3,2,2,3,0,0,0,0,1,1,4,0,0,0,0,0,5,0}; 
	
	float[] symfullc = new float[]{1,1,1,-1,0,0,1,-1,1,1,0,0,0,0,0,0,1,1};
	float[] cpsymfullc = new float[]{1,1,1,-1,0,0,1,-1,1,1,0,0,0,0,0,0,1,1};
	     
	double[] Bz = new double[] {7,1,11,1,15,1,8,0,12,1,16,2,9,0,13,1,17,1,10,0,14,1,18,0};//with dim K*N
        double[] Cz = new double[] {2,0,-3,1,-4,0,-7,1,3,0,0,1,7,0,18,0};//with dim M*N
        double[] cpCz = new double[] {2,0,-3,1,-4,0,-7,1,3,0,0,1,7,0,18,0};
	double[] myaz = new double[] {1,0,2,-1,3,0,4,1,5,1,6,0};//dim M*K
	double[] mybz = new double[] {-1,1,0,1,1,0,1,0,0,0,-1,2};//dim M*K
	double[] mycz = new double[] {1,1,0,0,0,0,1,-1};//dim M*M        
	double[] cpmyaz = new double[] {1,0,2,-1,3,0,4,1,5,1,6,0};//dim M*K
	double[] cpmybz = new double[] {-1,1,0,1,1,0,1,0,0,0,-1,2};//dim M*K
	double[] cpmycz = new double[] {1,1,0,0,0,0,1,-1};//dim M*M  
	double[] solvetriz = new double[] {3,2,2,3,0,0,0,0,1,1,4,0,0,0,0,0,5,0}; 
	
	double[] symfullz = new double[]{1,1,1,-1,0,0,1,-1,1,1,0,0,0,0,0,0,1,1};
	double[] cpsymfullz = new double[]{1,1,1,-1,0,0,1,-1,1,1,0,0,0,0,0,0,1,1};
	   
	float[] cpBc = new float[] {7,1,11,1,15,1,8,0,12,1,16,2,9,0,13,1,17,1,10,0,14,1,18,0};//with dim K*N
        double[] cpBz = new double[] {7,1,11,1,15,1,8,0,12,1,16,2,9,0,13,1,17,1,10,0,14,1,18,0};//with dim K*N
        
	float[] paramf = new float[]{0,0,0,0,0};
	float[] d1 = new float[]{1};
	float[] d2 = new float[]{1};
	float[] x1 = new float[]{1};
	float[] rotxf = new float[]{0,1};
	float[] rotyf = new float[]{1,2};
	
	double[] paramd = new double[]{0,0,0,0,0};
	double[] d1d = new double[]{1};
	double[] d2d = new double[]{1};
	double[] x1d = new double[]{1};
	double[] rotxd = new double[]{0,1};
	double[] rotyd = new double[]{1,2};
	//
        // Print the parameters
        //
        
	
	// sdsdot
	System.out.println();
        System.out.println("sdsdot: sdsdot = alpha + inner product of vector x and y");
	System.out.println("alpha=" + string(alpha));
        CprintMatrix("vector x = ", A, M*K, 1);
        CprintMatrix("vector y = ", cpA, M*K, 1);

	//N = M*K
	//alpha = alpha
	//X = A
	//incX = 1
	//Y = cpA
	//incY = 1

        resultf = jniCBLAS.sdsdot(M*K, alpha, A, 1, cpA, 1);
        System.out.printf("result = %f\n", resultf);
	
	// dsdot
	System.out.println();
        System.out.println("dsdot: dsdot = inner product of vector x and y");
	CprintMatrix("vector x = ", A, M*K, 1);
        CprintMatrix("vector y = ", cpA, M*K, 1);

	//N = M*k
	//X = A
	//incX = 1
	//Y = cpA
	//incY = 1

        resultd = jniCBLAS.dsdot(M*K, A, 1, cpA, 1);
        System.out.printf("result = %f\n", resultd);
	
	// sdot
	System.out.println();
        System.out.println("sdot: dsdot = inner product of vector x and y");
	CprintMatrix("vector x = ", A, M*K, 1);
        CprintMatrix("vector y = ", cpA, M*K, 1);

	//N = M*k
	//X = A
	//incX = 1
	//Y = cpA
	//incY = 1
	
	resultf = jniCBLAS.sdot(M*K, A, 1, cpA, 1);
        System.out.printf("result = %f\n", resultf);
	
	// ddot
	System.out.println();
        System.out.println("ddot: ddot = inner product of vector x and y");
	CprintMatrix("vector x = ", dA, M*K, 1);
        CprintMatrix("vector y = ", cpdA, M*K, 1);

	//N = M*k
	//X = dA
	//incX = 1
	//Y = cpdA
	//incY = 1

        resultd = jniCBLAS.ddot(M*K, dA, 1, cpdA, 1);
        System.out.printf("result = %f\n", resultd);

 	//cdotu_sub
	System.out.println();
        System.out.println("cdotu_sub: cdotu_sub = inner product of vector x and y");
	RprintMatrix("vector x = ", dotf, M, 2);
        RprintMatrix("vector y = ", cpdotf, M, 2);

	//N = M
	//X = dotf
	//incX = 1
	//Y = cpdotf
	//incY = 1
	//dotu = dotuf
        
	jniCBLAS.cdotu_sub(M, dotf, 1, cpdotf, 1, dotuf);
        System.out.printf("result = %f + %fi\n", dotuf[0], dotuf[1]);
	
	//cdotc_sub
	System.out.println();
        System.out.println("cdotc_sub = X^H * Y");
	RprintMatrix("vector x = ", dotf, M, 2);
        RprintMatrix("vector y = ", cpdotf, M, 2);
	
	//N = M
	//X = dotf
	//incX = 1
	//Y = cpdotf
	//incY = 1
	//dotc = dotuf
        
        jniCBLAS.cdotc_sub(M, dotf, 1, cpdotf, 1, dotuf);
        System.out.printf("result = %f + %fi\n", dotuf[0], dotuf[1]);

	//zdotu_sub
	System.out.println();
        System.out.println("zdotu_sub = X * Y");
	RprintMatrix("vector x = ", dotd, M, 2);
        RprintMatrix("vector y = ", cpdotd, M, 2);

	//N = M
	//X = dotd
	//incX = 1
	//Y = cpdotd
	//incY = 1
	//dotu = dotud
        
        jniCBLAS.zdotu_sub(M, dotd, 1, cpdotd, 1, dotud);
        System.out.printf("result = %f + %fi\n", dotud[0], dotud[1]);

	//zdotc_sub
	System.out.println();
        System.out.println("zdotc_sub = X * Y");
	RprintMatrix("vector x = ", dotd, M, 2);
        RprintMatrix("vector y = ", cpdotd, M, 2);

	//N = M
	//X = dotd
	//incX = 1
	//Y = cpdotd
	//incY = 1
	//dotc = dotud
        
        jniCBLAS.zdotc_sub(M, dotd, 1, cpdotd, 1, dotud);
        System.out.printf("result = %f + %fi\n", dotud[0], dotud[1]);

	//snrm2
	System.out.println();
        System.out.println("snrm2 = sqrt(x'*x) ");
	CprintMatrix("vector x = ", A, M*K, 1);

	//N = M*K
	//X = A
	//incX = 1
        
        resultf = jniCBLAS.snrm2(M*K, A, 1);
        System.out.printf("result = %f\n", resultf);

	//sasum
	System.out.println();
        System.out.println("sasum = sum of the absolute values ");
	CprintMatrix("vector x = ", A, M*K, 1);

	//N = M*K
	//X = A
	//incX = 1
        
	resultf = jniCBLAS.sasum(M*K, A, 1);
        System.out.printf("result = %f\n", resultf);

	//dnrm2
	System.out.println();
        System.out.println("dnrm2 = sqrt(x'*x) ");
	CprintMatrix("vector x = ", dA, M*K, 1);

	//N = M*K
	//X = dA
	//incX = 1
        
	resultd = jniCBLAS.dnrm2(M*K, dA, 1);
        System.out.printf("result = %f\n", resultd);

	//dasum
	System.out.println();
        System.out.println("dasum = sum of the absolute values ");
	CprintMatrix("vector x = ", dA, M*K, 1);

	//N = M*K
	//X = dA
	//incX = 1
                
	resultd = jniCBLAS.dasum(M*K, dA, 1);
        System.out.printf("result = %f\n", resultd);

	//scnrm2
	System.out.println();
        System.out.println("scnrm2 = sqrt(x**H*x) ");
	RprintMatrix("vector x = ", dotf, M, 2);

	//N = M
	//X = dotf
	//incX = 1
                
	resultf = jniCBLAS.scnrm2(M, dotf, 1);
        System.out.printf("result = %f\n", resultf);

	//scasum
	System.out.println();
        System.out.println("scasum = sum of the (|Re(.)|+|Im(.)|)'s of a complex vector");
	RprintMatrix("vector x = ", dotf, M, 2);

	//N = M
	//X = dotf
	//incX = 1
                
	resultf = jniCBLAS.scasum(M, dotf, 1);
        System.out.printf("result = %f\n", resultf);

	//dznrm2
	System.out.println();
        System.out.println("dznrm2 = sqrt(x**H*x) ");
	RprintMatrix("vector x = ", dotd, M, 2);

	//N = M
	//X = dotd
	//incX = 1
	
        resultd = jniCBLAS.dznrm2(M, dotd, 1);
        System.out.printf("result = %f\n", resultd);

	//dzasum
	System.out.println();
        System.out.println("dzasum = sum of the (|Re(.)|+|Im(.)|)'s of a complex vector");
	CprintMatrix("vector x = ", dotd, M, 2);

	//N = M
	//X = dotd
	//incX = 1
        
        resultd = jniCBLAS.dzasum(M, dotd, 1);
        System.out.printf("result = %f\n", resultd);

	//sswap
	System.out.println();
        System.out.println("sswap: interchanges two vectors");
	CprintMatrix("x = ", dotrf, M, 1);
        CprintMatrix("y = ", dotif, M, 1);

	//N = M
	//X = dotrf
	//incX = 1
	//Y = dotif
	//incY = 1
	
        jniCBLAS.sswap(M, dotrf, 1, dotif, 1);
	System.out.println("after swap: ");	
        CprintMatrix("x = ", dotrf, M, 1);
        CprintMatrix("y = ", dotif, M, 1);
        
	//scopy
	System.out.println();
        System.out.println("scopy: copy a vector x to a vector y");
	CprintMatrix("x = ", dotrf, M, 1);
        CprintMatrix("y = ", dotif, M, 1);

	//N = M
	//X = dotrf
	//incX = 1
	//Y = dotif
	//incY = 1
	
        jniCBLAS.scopy(M, dotrf, 1, dotif, 1);
	System.out.println("after copy: ");	
        CprintMatrix("x = ", dotrf, M, 1);
        CprintMatrix("y = ", dotif, M, 1);
        
	//saxpy
	System.out.println();
        System.out.println("saxpy: y = a * x + y ");
	System.out.println("a = " + string(alpha));
	CprintMatrix("x = ", cpdotrf, M, 1);
        CprintMatrix("y = ", dotif, M, 1);

	//N = M
	//alpha = alpha
	//X = cpdotrf
	//incX = 1
	//Y = dotif
	//incY = 1
	
        jniCBLAS.saxpy(M, alpha, cpdotrf, 1, dotif, 1);
	CprintMatrix("Result = ", dotif, M, 1);
        
	//dswap
	System.out.println();
        System.out.println("dswap: interchanges two vectors");
	CprintMatrix("x = ", dotr, M, 1);
        CprintMatrix("y = ", doti, M, 1);

	//N = M
	//X = dotr
	//incX = 1
	//Y = doti
	//incY = 1
	
        jniCBLAS.dswap(M, dotr, 1, doti, 1);
	System.out.println("after swap: ");	
        CprintMatrix("x = ", dotr, M, 1);
        CprintMatrix("y = ", doti, M, 1);
        
	//dcopy
	System.out.println();
        System.out.println("dcopy: copy a vector x to a vector y");
	CprintMatrix("x = ", dotr, M, 1);
        CprintMatrix("y = ", doti, M, 1);

	//N = M
	//X = dotr
	//incX = 1
	//Y = doti
	//incY = 1
	
        jniCBLAS.dcopy(M, dotr, 1, doti, 1);
	System.out.println("after copy: ");	
        CprintMatrix("x = ", dotr, M, 1);
        CprintMatrix("y = ", doti, M, 1);
        
	//daxpy
	System.out.println();
        System.out.println("daxpy: y = a * x + y ");
	System.out.println("a = " + string(alphad));
	CprintMatrix("x = ", cpdotr, M, 1);
	doti = cpdoti.clone();
        CprintMatrix("y = ", doti, M, 1);

	//N = M
	//alpha = alphad
	//X = cpdotr
	//incX = 1
	//Y = doti
	//incY = 1
	
        jniCBLAS.daxpy(M, alphad, cpdotr, 1, doti, 1);
	CprintMatrix("Result = ", doti, M, 1);
        
	//cswap
	System.out.println();
        System.out.println("cswap: interchanges two vectors");
	RprintMatrix("x = ", dotf, M, 2);
        RprintMatrix("y = ", cpdotf, M, 2);

	//N = M
	//X = dotf
	//incX = 1
	//Y = cpdotf
	//incY = 1
	
        jniCBLAS.cswap(M, dotf, 1, cpdotf, 1);
	System.out.println("after swap: ");	
        RprintMatrix("x = ", dotf, M, 2);
        RprintMatrix("y = ", cpdotf, M, 2);
        
	//ccopy
	System.out.println();
        System.out.println("ccopy: copy a vector x to a vector y");
	dotf = new float[] {1,2,3,4};
	cpdotf = new float[] {5,6,7,8};
	RprintMatrix("x = ", dotf, M, 2);
        RprintMatrix("y = ", cpdotf, M, 2);

	//N = M
	//X = dotf
	//incX = 1
	//Y = cpdotf
	//incY = 1
	
        jniCBLAS.ccopy(M, dotf, 1, cpdotf, 1);
	System.out.println("after copy: ");	
        RprintMatrix("x = ", dotf, M, 2);
        RprintMatrix("y = ", cpdotf, M, 2);
         
	//caxpy
	System.out.println();
        System.out.println("caxpy: y = a * x + y ");
	System.out.printf("a = %f + %fi\n", alpha, beta);
	RprintMatrix("x = ", dotf, M, 2);
        RprintMatrix("y = ", cpdotf, M, 2);

	//N = M
	//alpha = arrayconst
	//X = dotrf
	//incX = 1
	//Y = dotif
	//incY = 1
	        
	jniCBLAS.caxpy(M, arrayconst, dotf, 1, cpdotf, 1);
	System.out.println("after calculation: ");	
        RprintMatrix("y = ", cpdotf, M, 2);
        
	//zswap
	System.out.println();
        System.out.println("zswap: interchanges two vectors");
	RprintMatrix("x = ", dotd, M, 2);
        RprintMatrix("y = ", cpdotd, M, 2);

	//N = M
	//X = dotd
	//incX = 1
	//Y = cpdotd
	//incY = 1
	
        jniCBLAS.zswap(M, dotd, 1, cpdotd, 1);
	System.out.println("after swap: ");	
        RprintMatrix("x = ", dotd, M, 2);
        RprintMatrix("y = ", cpdotd, M, 2);
        
	//zcopy
	System.out.println();
        System.out.println("zcopy: copy a vector x to a vector y");
	dotd = new double[] {1,2,3,4};
	cpdotd = new double[] {5,6,7,8};
	RprintMatrix("x = ", dotd, M, 2);
        RprintMatrix("y = ", cpdotd, M, 2);

	//N = M
	//X = dotd
	//incX = 1
	//Y = cpdotd
	//incY = 1
	        
	jniCBLAS.zcopy(M, dotd, 1, cpdotd, 1);
	System.out.println("after copy: ");	
        RprintMatrix("x = ", dotd, M, 2);
        RprintMatrix("y = ", cpdotd, M, 2);
   
	//zaxpy
	System.out.println();
        System.out.println("zaxpy: y = a * x + y ");
	System.out.printf("a = %f + %fi\n", alphad, betad);
	RprintMatrix("x = ", dotd, M, 2);
        RprintMatrix("y = ", cpdotd, M, 2);

	//N = M
	//alpha = arrayconstd
	//X = dotd
	//incX = 1
	//Y = cpdotd
	//incY = 1
	        
	jniCBLAS.zaxpy(M, arrayconstd, dotd, 1, cpdotd, 1);
	System.out.println("after calculation: ");	
        RprintMatrix("y = ", cpdotd, M, 2);
             
	//srotg
	System.out.println();
        System.out.println("srotg construct givens plane rotation");
	System.out.println("a = 8, b = 6");
	
	//a = af
	//b = bf
	//c = arrayresult1
	//s = arrayresult2
	
	jniCBLAS.srotg(af, bf, arrayresult1, arrayresult2);
	System.out.println("Output: ");	
        System.out.printf("a = %f\tb=%f\tc=%f\td=%f\n", af[0], bf[0], arrayresult1[0], arrayresult2[0]);	
        
	//srotmg
	System.out.println();
        System.out.println("srotmg construct the modified givens transformation matrix which zeros the second component of the 2-vector");
	System.out.println("d1 = 1, d2 = 1, x1 = 1, y1 = 1");

	//d1 = d1
	//d2 = d2
	//b1 = x1
	//b2 = 1
	//p = paramf

	jniCBLAS.srotmg(d1,d2,x1,1,paramf);
	System.out.println("Result: ");	
        RprintMatrix("param = ",paramf,5,1);
        
	//srot: applies a plane rotation
	System.out.println();
        System.out.println("srot: applies a plane rotation");
	System.out.println("n = 5, incx = 1, incy = 2, c = 0.5, s = sqrt(3)/2");
	RprintMatrix("x = ", tmpx, 1, 5);
        RprintMatrix("y = ", tmpy, 1, 5);

	//N = 5
	//X = tmpx
	//incX = 1
	//Y = tmpy
	//incY = 2
	//c = 0.5f
	//s = 0.8660254f

        jniCBLAS.srot(5,tmpx,1,tmpy,2,0.5f,0.8660254f);
	System.out.println("Output: ");	
        RprintMatrix("x = ", tmpx, 1, 5);
        RprintMatrix("y = ", tmpy, 1, 9);

	//srotm
	System.out.println();
        System.out.println("srotm: applies the modified givens transformation, H, to the 2 by n matrix");
	System.out.println("n = 2");
	RprintMatrix("x = ", rotxf, 2, 1);
        RprintMatrix("y = ", rotyf, 2, 1);
	RprintMatrix("param = ",paramf,1,5);

	//N = 2
	//X = rotxf
	//incX = 1
	//Y = rotyf
	//incY = 1
	//p = paramf

        jniCBLAS.srotm(2,rotxf,1,rotyf,1,paramf);
	System.out.println("Output: ");	
        RprintMatrix("x = ", rotxf, 2, 1);
        RprintMatrix("y = ", rotyf, 2, 1);

	//drotg
	System.out.println();
        System.out.println("drotg construct givens plane rotation");
	System.out.println("a = 8, b = 6");
	
	//a = ad
	//b = bd
	//c = arrayresult1d
	//s = arrayresult2d

	jniCBLAS.drotg(ad, bd, arrayresult1d, arrayresult2d);
	System.out.println("Output: ");	
        System.out.printf("a = %f\tb=%f\tc=%f\td=%f\n", ad[0], bd[0], arrayresult1d[0], arrayresult2d[0]);

	//drotmg
	System.out.println();
        System.out.println("drotmg construct the modified givens transformation matrix which zeros the second component of the 2-vector");
	System.out.println("d1 = 1, d2 = 1, x1 = 1, y1 = 1");

	//d1 = d1d
	//d2 = d2d
	//b1 = x1d
	//b2 = 1
	//p = paramd

	jniCBLAS.drotmg(d1d,d2d,x1d,1,paramd);
	System.out.println("Result: ");	
        RprintMatrix("param = ",paramd,5,1);
        
        
	//drot: applies a plane rotation
	System.out.println();
        System.out.println("drot: applies a plane rotation");
	System.out.println("n = 5, incx = 1, incy = 2, c = 0.5, s = sqrt(3)/2");
	RprintMatrix("x = ", tmpxd, 1, 5);
        RprintMatrix("y = ", tmpyd, 1, 5);

	//N = 5
	//X = tmpxd
	//incX = 1
	//Y = tmpyd
	//incY = 2
	//c = 0.5
	//s = 0.8660254

        jniCBLAS.drot(5,tmpxd,1,tmpyd,2,0.5,0.8660254);
	System.out.println("Output: ");	
        RprintMatrix("x = ", tmpxd, 1, 5);
        RprintMatrix("y = ", tmpyd, 1, 9);

	//drotm
	System.out.println();
        System.out.println("drotm: applies the modified givens transformation, H, to the 2 by n matrix");
	System.out.println("n = 2");
	RprintMatrix("x = ", rotxd, 2, 1);
        RprintMatrix("y = ", rotyd, 2, 1);
	RprintMatrix("param = ",paramd,1,5);

	//N = 2
	//X = rotxd
	//incX = 1
	//Y = rotyd
	//incY = 1
	//p = paramd

        jniCBLAS.drotm(2,rotxd,1,rotyd,1,paramd);
	System.out.println("Output: ");	
        RprintMatrix("x = ", rotxd, 2, 1);
        RprintMatrix("y = ", rotyd, 2, 1);

	//sscal
	System.out.println();
        System.out.println("sscal: scales a vector by a constant");
	CprintMatrix("vector x = ", A, M*K, 1);
	System.out.println("alpha=" + string(alpha));

	//N = M*K
	//alpha = alpha	
	//x = A
	//incx = 1

        jniCBLAS.sscal(M*K, alpha, A, 1);
        CprintMatrix("Afterscaling: vector x = ", A, M*K, 1);
	
	//dscal
	System.out.println();
        System.out.println("dscal: scales a vector by a constant");
	CprintMatrix("vector x = ", dA, M*K, 1);
	System.out.println("alpha=" + string(alpha));

	//N = M*K
	//alpha = alphad	
	//x = dA
	//incx = 1

        jniCBLAS.dscal(M*K, alphad, dA, 1);
        CprintMatrix("Afterscaling: vector x = ", dA, M*K, 1);
	
	//cscal
	System.out.println();
        System.out.println("cscal: scales a vector by a constant");
	RprintMatrix("vector x = ", dotf, M, 2);
	RprintMatrix("alpha = ", arrayconst, 1, 2);

	//N = M
	//alpha = arrayconst	
	//x = dotf
	//incx = 1

	jniCBLAS.cscal(M, arrayconst, dotf, 1);
        RprintMatrix("Afterscaling: vector x = ", dotf, M, 2);
	
	//zscal
	System.out.println();
        System.out.println("zscal: scales a vector by a constant");
	RprintMatrix("vector x = ", dotd, M, 2);
	RprintMatrix("alpha = ", arrayconstd, 1, 2);

	//N = M
	//alpha = arrayconstd	
	//x = dotd
	//incx = 1

	jniCBLAS.zscal(M, arrayconstd, dotd, 1);
        RprintMatrix("Afterscaling: vector x = ", dotd, M, 2);
	
	//csscal
	System.out.println();
        System.out.println("csscal: scales a complex vector by a real constant");
	RprintMatrix("vector x = ", dotf, M, 2);
	System.out.printf("alpha = %f\n", arrayconst[0]);

	//N = M
	//alpha = arrayconst[0]	
	//x = dotf
	//incx = 1

        jniCBLAS.csscal(M, arrayconst[0], dotf, 1);
        RprintMatrix("Afterscaling: vector x = ", dotf, M, 2);
	
	//zdscal
	System.out.println();
        System.out.println("zdscal: scales a complex vector by a real constant");
	RprintMatrix("vector x = ", dotd, M, 2);
	System.out.printf("alpha = %f\n", arrayconstd[0]);

	//N = M
	//alpha = arrayconstd[0]	
	//x = dotd
	//incx = 1

        jniCBLAS.zdscal(M, arrayconstd[0], dotd, 1);
        RprintMatrix("Afterscaling: vector x = ", dotd, M, 2);

	//sgemv
	System.out.println();
        System.out.println("sgemv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
        RprintMatrix("A = ", dotf, M, M);
	RprintMatrix("x = ", dotrf, M, 1);
	RprintMatrix("y = ", dotif, M, 1);
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.Trans;

	//order = order
	//TransA = transA
	//M = M
	//N = M
	//alpha = alpha	
	//a = dotf
	//lda = M
	//x = dotrf
	//incx = 1
	//beta = beta
	//y = dotif
	//incy = 1

	jniCBLAS.sgemv(order, transA, M, M, alpha, dotf, M, dotrf, 1, beta, dotif, 1);
        RprintMatrix("Result: y = ", dotif, M, 1);
	//note: if you input a matrix as row-major, then the meaning of lda is also changed

	//sgbmv
	System.out.println();
        System.out.println("sgbmv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y, and A is a band matrix");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
        RprintMatrix("A = ", bandfull, K, K);
	RprintMatrix("x = ", xf, K, 1);
	RprintMatrix("y = ", yf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//order = order
	//TransA = transA
	//M = K
	//N = K
	//KL = 1
	//KU = 1
	//alpha = alpha	
	//a = bandf
	//lda = K
	//x = xf
	//incx = 1
	//beta = beta
	//y = yf
	//incy = 1

	jniCBLAS.sgbmv(order, transA, K, K, 1, 1, alpha, bandf, K, xf, 1, beta, yf, 1);
        RprintMatrix("Result: y = ", yf, K, 1);
	
	//strmv
	System.out.println();
        System.out.println("strmv: x = A * x or x = A**T * x");
	CprintMatrix("A = ", trif, K, K);
	RprintMatrix("x = ", xf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//a = trif
	//lda = K
	//x = xf
	//incx = 1
	
	jniCBLAS.strmv(order, uplo, transA, diag, K, trif, K, xf, 1);
        RprintMatrix("Result: x = ", xf, K, 1);
	
	//stbmv
	System.out.println();
        System.out.println("stbmv: x = A * x or x = A**T * x, A is triangular matrix in banded storage");
	CprintMatrix("A = ", tribandfull, K, K);
	xf = cpxf.clone();
	RprintMatrix("x = ", xf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//K = 1	
	//a = tribandf
	//lda = K
	//x = xf
	//incx = 1

	jniCBLAS.stbmv(order, uplo, transA, diag, K, 1, tribandf, K, xf, 1);
        RprintMatrix("Result: x = ", xf, K, 1);
	
	//stpmv
	System.out.println();
        System.out.println("stpmv: x = A * x or x = A**T * x, A is triangular banded matrix in packed storage");
	CprintMatrix("A = ", tribandfull, K, K);
	xf = cpxf.clone();
	RprintMatrix("x = ", xf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//Ap = packed
	//x = xf
	//incx = 1

	jniCBLAS.stpmv(order, uplo, transA, diag, K, packed, xf, 1);
        RprintMatrix("Result: x = ", xf, K, 1);
	
	//strsv
	System.out.println();
        System.out.println("strsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular matrix");
	CprintMatrix("A = ", tribandfull, K, K);
	xf = cpxf.clone();
	RprintMatrix("b = ", xf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//A = tribandfull
	//lda = K
	//x = xf
	//incX = 1

	jniCBLAS.strsv(order, uplo, transA, diag, K, tribandfull, K, xf, 1);
        RprintMatrix("Result: x = ", xf, K, 1);
	
	//stbsv
	System.out.println();
        System.out.println("stbsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular banded matrix");
	CprintMatrix("A = ", tribandfull, K, K);
	xf = cpxf.clone();
	RprintMatrix("b = ", xf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//K = 1
	//A = tribandf
	//lda = K
	//x = xf
	//incX = 1

	jniCBLAS.stbsv(order, uplo, transA, diag, K, 1, tribandf, K, xf, 1);
        RprintMatrix("Result: x = ", xf, K, 1);
	
	//stpsv
	System.out.println();
        System.out.println("stpsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular matrix in packed form");
	CprintMatrix("A = ", tribandfull, K, K);
	xf = cpxf.clone();
	RprintMatrix("b = ", xf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//Ap = packed
	//x = xf
	//incX = 1
	
	jniCBLAS.stpsv(order, uplo, transA, diag, K, packed, xf, 1);
        RprintMatrix("Result: x = ", xf, K, 1);
	
	//dgemv
	System.out.println();
        System.out.println("dgemv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
	dotd = cpdotd.clone();
	doti = cpdoti.clone();
	dotr = cpdotr.clone();
        RprintMatrix("A = ", dotd, M, M);
	RprintMatrix("x = ", dotr, M, 1);
	RprintMatrix("y = ", doti, M, 1);
	order = jniCBLAS.ORDER.RowMajor; 	

	//order = order
	//TransA = transA
	//M = M
	//N = M
	//alpha = alpha
	//A = dotd
	//lda = M
	//x = dotr
	//incX = 1
	//beta = beta
	//y = doti
	//incY = 1

	transA = jniCBLAS.TRANSPOSE.Trans;
	jniCBLAS.dgemv(order, transA, M, M, alpha, dotd, M, dotr, 1, beta, doti, 1);
        RprintMatrix("Result: y = ", doti, M, 1);
	
	//dgbmv
	System.out.println();
        System.out.println("dgbmv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y, and A is a band matrix");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
        RprintMatrix("A = ", bandfull, K, K);
	RprintMatrix("x = ", xd, K, 1);
	RprintMatrix("y = ", yd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//order = order
	//TransA = transA
	//M = K
	//N = K
	//KL = 1
	//KU = 1
	//alpha = alpha
	//A = bandd
	//lda = K
	//x = xd
	//incX = 1
	//beta = beta
	//y = yd
	//incY = 1

	jniCBLAS.dgbmv(order, transA, K, K, 1, 1, alpha, bandd, K, xd, 1, beta, yd, 1);
        RprintMatrix("Result: y = ", yd, K, 1);

	//dtrmv
	System.out.println();
        System.out.println("dtrmv: x = A * x or x = A**T * x");
	CprintMatrix("A = ", trid, K, K);
	RprintMatrix("x = ", xd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//N = K	
	//A = trid
	//lda = K
	//x = xd
	//incX = 1
	
	jniCBLAS.dtrmv(order, uplo, transA, diag, K, trid, K, xd, 1);
        RprintMatrix("Result: x = ", xd, K, 1);
	
	//dtbmv
	System.out.println();
        System.out.println("dtbmv: x = A * x or x = A**T * x, A is triangular matrix in banded storage");
	CprintMatrix("A = ", tribandfull, K, K);
	xd = cpxd.clone();
	RprintMatrix("x = ", xd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//N = K	
	//K = 1
	//A = tribandd
	//lda = K
	//x = xd
	//incX = 1

	jniCBLAS.dtbmv(order, uplo, transA, diag, K, 1, tribandd, K, xd, 1);
        RprintMatrix("Result: x = ", xd, K, 1);
	
	//dtpmv
	System.out.println();
        System.out.println("dtpmv: x = A * x or x = A**T * x, A is triangular banded matrix in packed storage");
	CprintMatrix("A = ", tribandfull, K, K);
	xd = cpxd.clone();
	RprintMatrix("x = ", xd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//N = K	
	//Ap = packedd
	//x = xd
	//incX = 1

	jniCBLAS.dtpmv(order, uplo, transA, diag, K, packedd, xd, 1);
        RprintMatrix("Result: x = ", xd, K, 1);
	
	//dtrsv
	System.out.println();
        System.out.println("dtrsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular matrix");
	CprintMatrix("A = ", tribandfull, K, K);
	xd = cpxd.clone();
	RprintMatrix("b = ", xd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//N = K	
	//A = tribandfulld
	//lda = K
	//x = xd
	//incX = 1

	jniCBLAS.dtrsv(order, uplo, transA, diag, K, tribandfulld, K, xd, 1);
        RprintMatrix("Result: x = ", xd, K, 1);
	
	//dtbsv
	System.out.println();
        System.out.println("dtbsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular banded matrix");
	CprintMatrix("A = ", tribandfulld, K, K);
	xd = cpxd.clone();
	RprintMatrix("b = ", xd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//N = K	
	//K = 1
	//A = tribandd
	//lda = K
	//x = xd
	//incX = 1

	jniCBLAS.dtbsv(order, uplo, transA, diag, K, 1, tribandd, K, xd, 1);
        RprintMatrix("Result: x = ", xd, K, 1);
	
	//dtpsv
	System.out.println();
        System.out.println("dtpsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular matrix in packed form");
	CprintMatrix("A = ", tribandfulld, K, K);
	xd = cpxd.clone();
	RprintMatrix("b = ", xd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//N = K	
	//Ap = packedd
	//x = xd
	//incX = 1

	jniCBLAS.dtpsv(order, uplo, transA, diag, K, packedd, xd, 1);
        RprintMatrix("Result: x = ", xd, K, 1);

	//cgemv
	System.out.println();
        System.out.println("cgemv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y");
	System.out.printf("alpha = %f + %fi\n", alphac[0], alphac[1]);
        System.out.printf("beta = %f + %fi\n", betac[0], betac[1]);
        RprintMatrix("A = ", ac, M, 2*M);
	RprintMatrix("x = ", xc, M, 2);
	RprintMatrix("y = ", yc, M, 2);
	System.out.println("Case 1 : order = RowMajor");
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//order = order
	//TransA = transA
	//M = M
	//N = M
	//alpha = alphac
	//A = ac
	//lda = M
	//x = xc
	//incX = 1
	//beta = betac
	//y = yc
	//incY = 1

	jniCBLAS.cgemv(order, transA, M, M, alphac, ac, M, xc, 1, betac, yc, 1);
        RprintMatrix("Result: y = ", yc, M, 2);
	System.out.println("Case 2 : order = ColMajor");
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	ac = new float[]{1,1,2,1,1,-1,2,-1};
	yc = cpyc.clone();

	//order = order
	//TransA = transA
	//M = M
	//N = M
	//alpha = alphac
	//A = ac
	//lda = M
	//x = xc
	//incX = 1
	//beta = betac
	//y = yc
	//incY = 1
	
	jniCBLAS.cgemv(order, transA, M, M, alphac, ac, M, xc, 1, betac, yc, 1);
        RprintMatrix("Result: y = ", yc, M, 2);
	System.out.println();
	System.out.printf("alpha = %f + %fi\n", alphac[0], alphac[1]);
        System.out.printf("beta = %f + %fi\n", betac[0], betac[1]);
        RprintMatrix("A = ", bandfullc, K, 2*K);
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);
	System.out.println("Case 1 : order = ColMajor");
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//order = order
	//TransA = transA
	//M = K
	//N = K
	//alpha = alphac
	//A = bandfullccol
	//lda = K
	//x = bxc
	//incX = 1
	//beta = betac
	//y = byc
	//incY = 1

	jniCBLAS.cgemv(order, transA, K, K, alphac, bandfullccol, K, bxc, 1, betac, byc, 1);
        RprintMatrix("Result: y = ", byc, K, 2);
	System.out.println("Case 2 : order = RowMajor");
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	byc = cpbyc.clone();

	//order = order
	//TransA = transA
	//M = K
	//N = K
	//alpha = alphac
	//A = bandfullc
	//lda = K
	//x = bxc
	//incX = 1
	//beta = betac
	//y = byc
	//incY = 1

	jniCBLAS.cgemv(order, transA, K, K, alphac, bandfullc, K, bxc, 1, betac, byc, 1);
        RprintMatrix("Result: y = ", byc, K, 2);
	
	//cgbmv
	System.out.println();
        System.out.println("cgbmv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y, and A is a band matrix");
	System.out.printf("alpha = %f + %fi\n", alphac[0], alphac[1]);
        System.out.printf("beta = %f + %fi\n", betac[0], betac[1]);
        RprintMatrix("A = ", bandfullc, K, 2*K);
	System.out.println("Case 1 : order = ColMajor");
	RprintMatrix("x = ", bxc, K, 2);
	byc = cpbyc.clone();
	RprintMatrix("y = ", byc, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//order = order
	//TransA = transA
	//M = K
	//N = K
	//KL = 1
	//KU = 1
	//alpha = alphac
	//A = bandccol
	//lda = K
	//x = bxc
	//incX = 1
	//beta = betac
	//y = byc
	//incY = 1

	jniCBLAS.cgbmv(order, transA, K, K, 1, 1, alphac, bandccol, K, bxc, 1, betac, byc, 1);
        RprintMatrix("Result: y = ", byc, K, 2);
	//note: in case of row-major, band storage is different from that of col-majored ones
	System.out.println("Case 2 : order = RowMajor");
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	byc = cpbyc.clone();
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);

	//order = order
	//TransA = transA
	//M = K
	//N = K
	//KL = 1
	//KU = 1
	//alpha = alphac
	//A = bandcrow
	//lda = K
	//x = bxc
	//incX = 1
	//beta = betac
	//y = byc
	//incY = 1

	jniCBLAS.cgbmv(order, transA, K, K, 1, 1, alphac, bandcrow, K, bxc, 1, betac, byc, 1);
        RprintMatrix("Result: y = ", byc, K, 2);

	//ctrmv
	System.out.println();
        System.out.println("ctrmv: x = A * x or x = A**T * x");
	RprintMatrix("A = ", tric, K, 2*K);
	bxc = cpbxc.clone();
	RprintMatrix("x = ", bxc, K, 2);
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//A = tric
	//lda = K
	//x = bxc
	//incX = 1

	jniCBLAS.ctrmv(order, uplo, transA, diag, K, tric, K, bxc, 1);
        RprintMatrix("Result: x = ", bxc, K, 2);
	
	//ctbmv
	System.out.println();
        System.out.println("ctbmv: x = A * x or x = A**T * x, A is triangular matrix in banded storage");
	RprintMatrix("A = ", tric, K, 2*K);
	bxc = cpbxc.clone();
	RprintMatrix("x = ", bxc, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//K = 1
	//A = tribandc
	//lda = K
	//x = bxc
	//incX = 1

	jniCBLAS.ctbmv(order, uplo, transA, diag, K, 1, tribandc, K, bxc, 1);
        RprintMatrix("Result: x = ", bxc, K, 2);
	
	//ctpmv
	System.out.println();
        System.out.println("ctpmv: x = A * x or x = A**T * x, A is triangular banded matrix in packed storage");
	RprintMatrix("A = ", tric, K, 2*K);
	bxc = cpbxc.clone();
	RprintMatrix("x = ", bxc, K, 2);
	System.out.println("Case 1 : col-major");
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//Ap = packedc
	//x = bxc
	//incX = 1

	jniCBLAS.ctpmv(order, uplo, transA, diag, K, packedc, bxc, 1);
        RprintMatrix("Result: x = ", bxc, K, 2);
	System.out.println("Case 2 : row-major");
	bxc = cpbxc.clone();
	packedc = new float[]{1,1,2,0,1,-1,1,1,3,0,1,-1,1,1};
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//Ap = packedc
	//x = bxc
	//incX = 1

	jniCBLAS.ctpmv(order, uplo, transA, diag, K, packedc, bxc, 1);
        RprintMatrix("Result: x = ", bxc, K, 2);
	
	//ctrsv
	System.out.println();
        System.out.println("ctrsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular matrix");
	RprintMatrix("A = ", tric, K, 2*K);
	bxc = cpbxc.clone();
	RprintMatrix("b = ", bxc, K, 2);
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//A = tric
	//lda = K
	//x = bxc
	//incX = 1

	jniCBLAS.ctrsv(order, uplo, transA, diag, K, tric, K, bxc, 1);
        RprintMatrix("Result: x = ", bxc, K, 2);
	
	//ctbsv
	System.out.println();
        System.out.println("ctbsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular banded matrix");
	RprintMatrix("A = ", tric, K, 2*K);
	bxc = cpbxc.clone();
	RprintMatrix("b = ", bxc, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//K = 1
	//A = tribandc
	//lda = K
	//x = bxc
	//incX = 1

	jniCBLAS.ctbsv(order, uplo, transA, diag, K, 1, tribandc, K, bxc, 1);
        RprintMatrix("Result: x = ", bxc, K, 2);
	
	//ctpsv
	System.out.println();
        System.out.println("ctpsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular banded matrix");
	RprintMatrix("A = ", tric, K, 2*K);
	bxc = cpbxc.clone();
	RprintMatrix("b = ", bxc, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//Ap = packedc
	//x = bxc
	//incX = 1

	jniCBLAS.ctpsv(order, uplo, transA, diag, K, packedc, bxc, 1);
        RprintMatrix("Result: x = ", bxc, K, 2);
	
	//zgemv
	System.out.println();
        System.out.println("zgemv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y");
	System.out.printf("alpha = %f + %fi\n", alphac[0], alphac[1]);
        System.out.printf("beta = %f + %fi\n", betac[0], betac[1]);
        RprintMatrix("A = ", az, M, 2*M);
	RprintMatrix("x = ", xz, M, 2);
	RprintMatrix("y = ", yz, M, 2);
	System.out.println("Case 1 : order = RowMajor");
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//order = order
	//TransA = transA
	//M = M
	//N = M
	//alpha = alphaz
	//A = az
	//lda = M
	//x = xz
	//incX = 1
	//beta = betaz
	//y = yz
	//incY = 1

	jniCBLAS.zgemv(order, transA, M, M, alphaz, az, M, xz, 1, betaz, yz, 1);
        RprintMatrix("Result: y = ", yz, M, 2);
	System.out.println("Case 2 : order = ColMajor");
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	az = new double[]{1,1,2,1,1,-1,2,-1};
	yz = cpyz.clone();

	//order = order
	//TransA = transA
	//M = M
	//N = M
	//alpha = alphaz
	//A = az
	//lda = M
	//x = xz
	//incX = 1
	//beta = betaz
	//y = yz
	//incY = 1

	jniCBLAS.zgemv(order, transA, M, M, alphaz, az, M, xz, 1, betaz, yz, 1);
        RprintMatrix("Result: y = ", yz, M, 2);

	//zgbmv
	System.out.println();
        System.out.println("zgbmv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y, and A is a band matrix");
	System.out.printf("alpha = %f + %fi\n", alphaz[0], alphaz[1]);
        System.out.printf("beta = %f + %fi\n", betaz[0], betaz[1]);
        RprintMatrix("A = ", bandfullz, K, 2*K);
	RprintMatrix("x = ", bxz, K, 2);
	RprintMatrix("y = ", byz, K, 2);
	System.out.println("Case 1 : order = ColMajor");
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	byc = cpbyc.clone();

	//order = order
	//TransA = transA
	//M = K
	//N = K
	//KL = 1
	//KU = 1
	//alpha = alphaz
	//A = bandzcol
	//lda = K
	//x = bxz
	//incX = 1
	//beta = betaz
	//y = byz
	//incY = 1

	jniCBLAS.zgbmv(order, transA, K, K, 1, 1, alphaz, bandzcol, K, bxz, 1, betaz, byz, 1);
        RprintMatrix("Result: y = ", byz, K, 2);
	
	//ztrmv
	System.out.println();
        System.out.println("ztrmv: x = A * x or x = A**T * x");
	RprintMatrix("A = ", triz, K, 2*K);
	bxz = cpbxz.clone();
	RprintMatrix("x = ", bxz, K, 2);
	order = jniCBLAS.ORDER.RowMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//A = triz
	//lda = K
	//x = bxz
	//incX = 1

	jniCBLAS.ztrmv(order, uplo, transA, diag, K, triz, K, bxz, 1);
        RprintMatrix("Result: x = ", bxz, K, 2);
	
	//ztbmv
	//the case of row-major...
	System.out.println();
        System.out.println("ztbmv: x = A * x or x = A**T * x, A is triangular matrix in banded storage");
	RprintMatrix("A = ", triz, K, 2*K);
	bxz = cpbxz.clone();
	RprintMatrix("x = ", bxz, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//K = 1
	//A = tribandz
	//lda = K
	//x = bxz
	//incX = 1

	jniCBLAS.ztbmv(order, uplo, transA, diag, K, 1, tribandz, K, bxz, 1);
        RprintMatrix("Result: x = ", bxz, K, 2);

	//ztpmv
	System.out.println();
        System.out.println("ztpmv: x = A * x or x = A**T * x, A is triangular banded matrix in packed storage");
	RprintMatrix("A = ", triz, K, 2*K);
	bxz = cpbxz.clone();
	RprintMatrix("x = ", bxz, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//Ap = packedz
	//x = bxz
	//incX = 1

	jniCBLAS.ztpmv(order, uplo, transA, diag, K, packedz, bxz, 1);
        RprintMatrix("Result: x = ", bxz, K, 2);
	
	//ztbsv
	System.out.println();
        System.out.println("ztbsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular banded matrix");
	RprintMatrix("A = ", triz, K, 2*K);
	bxz = cpbxz.clone();
	RprintMatrix("b = ", bxz, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//K = 1
	//A = tribandz
	//lda = K
	//x = bxz
	//incX = 1

	jniCBLAS.ztbsv(order, uplo, transA, diag, K, 1, tribandz, K, bxz, 1);
        RprintMatrix("Result: x = ", bxz, K, 2);

	//ztpsv
	System.out.println();
        System.out.println("ztpsv: solves systems of equations: A * x = b, or A**T * x = b, where A is triangular banded matrix");
	RprintMatrix("A = ", triz, K, 2*K);
	bxz = cpbxz.clone();
	RprintMatrix("b = ", bxz, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	diag = jniCBLAS.DIAG.NonUnit;

	//order = order
	//Uplo = uplo
	//TransA = transA
	//Diag = diag	
	//N = K
	//Ap = packedz
	//x = bxz
	//incX = 1

	jniCBLAS.ztpsv(order, uplo, transA, diag, K, packedz, bxz, 1);
        RprintMatrix("Result: x = ", bxz, K, 2);
	
	//ssymv
	System.out.println();
        System.out.println("ssymv: y = alpha * A * x + beta * y, where A is symmetric matrix");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
	CprintMatrix("A = ", symfulls, K, K);
	xf = cpxf.clone();
	yf = cpyf.clone();
	RprintMatrix("x = ", xf, K, 1);
	RprintMatrix("y = ", yf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//order = order
	//Uplo = uplo
	//N = K
	//alpha = alpha
	//A = symfulls
	//lda = K
	//x = xf
	//incX = 1
	//beta = beta
	//y = yf
	//incY = 1

	jniCBLAS.ssymv(order, uplo, K, alpha, symfulls, K, xf, 1, beta, yf, 1);
        RprintMatrix("Result: x = ", yf, K, 1);
	
	//ssbmv
	System.out.println();
        System.out.println("ssbmv: y = alpha * A * x + beta * y, where A is symmetric matrix");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
	CprintMatrix("A = ", symfulls, K, K);
	xf = cpxf.clone();
	yf = cpyf.clone();
	RprintMatrix("x = ", xf, K, 1);
	RprintMatrix("y = ", yf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//order = order
	//Uplo = uplo
	//N = K
	//K = 1
	//alpha = alpha
	//A = symbands
	//lda = K
	//x = xf
	//incX = 1
	//beta = beta
	//y = yf
	//incY = 1

	jniCBLAS.ssbmv(order, uplo, K, 1, alpha, symbands, K, xf, 1, beta, yf, 1);
        RprintMatrix("Result: x = ", yf, K, 1);
	
	//sspmv
	System.out.println();
        System.out.println("sspmv: y = alpha * A * x + beta * y, where A is symmetric matrix");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
	CprintMatrix("A = ", symfulls, K, K);
	xf = cpxf.clone();
	yf = cpyf.clone();
	RprintMatrix("x = ", xf, K, 1);
	RprintMatrix("y = ", yf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//order = order
	//Uplo = uplo
	//N = K
	//alpha = alpha
	//Ap = packedsyms
	//x = xf
	//incX = 1
	//beta = beta
	//y = yf
	//incY = 1

	jniCBLAS.sspmv(order, uplo, K, alpha, packedsyms, xf, 1, beta, yf, 1);
        RprintMatrix("Result: x = ", yf, K, 1);
	
	//sger
	System.out.println();
        System.out.println("sger: A = alpha * x * y**T + A");
	System.out.printf("alpha = %f\n", alpha);
        CprintMatrix("A = ", symfulls, K, K);
	xf = cpxf.clone();
	yf = cpyf.clone();
	RprintMatrix("x = ", xf, K, 1);
	RprintMatrix("y = ", yf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//M = K 
	//N = K
	//alpha = alpha
	//x = xf
	//incx = 1
	//y = yf
	//incy = 1
	//a = symfulls
	//lda = K

	jniCBLAS.sger(order, K, K, alpha, xf, 1, yf, 1, symfulls, K);
        CprintMatrix("Result: A = ", symfulls, K, K);
	
	//ssyr
	System.out.println();
        System.out.println("ssyr: A = alpha * x * x**T + A");
	System.out.printf("alpha = %f\n", alpha);
	symfulls = cpsymfulls.clone();
        CprintMatrix("A = ", symfulls, K, K);
	xf = cpxf.clone();
	RprintMatrix("x = ", xf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alpha
	//x = xf
	//incx = 1
	//a = symfulls
	//lda = K
	
	jniCBLAS.ssyr(order, uplo, K, alpha, xf, 1, symfulls, K);
        CprintMatrix("Result: A = ", symfulls, K, K);
	
	//sspr
	System.out.println();
        System.out.println("sspr: A = alpha * x * x**T + A");
	System.out.printf("alpha = %f\n", alpha);
	symfulls = cpsymfulls.clone();
        CprintMatrix("A = ", symfulls, K, K);
	xf = cpxf.clone();
	RprintMatrix("x = ", xf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;
	packedsyms = cppackedsyms.clone();

	//uplo = uplo
	//N = K
	//alpha = alpha
	//x = xf
	//incx = 1
	//Ap = packedsyms

	jniCBLAS.sspr(order, uplo, K, alpha, xf, 1, packedsyms);
        CprintMatrix("Result: A = ", packedsyms, 1, K*(K+1)/2);
	
	//ssyr2
	System.out.println();
        System.out.println("ssyr2: A = alpha * x * y**T + alpha * y * x**T + A");
	System.out.printf("alpha = %f\n", alpha);
	symfulls = cpsymfulls.clone();
        CprintMatrix("A = ", symfulls, K, K);
	xf = cpxf.clone();
	yf = cpyf.clone();
	RprintMatrix("x = ", xf, K, 1);
	RprintMatrix("y = ", yf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alpha
	//x = xf	
	//incx = 1
	//y = yf
	//incy = 1
	//a = symfulls
	//lda = K
	
	jniCBLAS.ssyr2(order, uplo, K, alpha, xf, 1, yf, 1, symfulls, K);
        CprintMatrix("Result: A = ", symfulls, K, K);
	
	//sspr2
	System.out.println();
        System.out.println("sspr2: A = alpha * x * y**T + alpha * y * x**T + A");
	System.out.printf("alpha = %f\n", alpha);
	symfulls = cpsymfulls.clone();
        CprintMatrix("A = ", symfulls, K, K);
	xf = cpxf.clone();
	yf = cpyf.clone();
	RprintMatrix("x = ", xf, K, 1);
	RprintMatrix("y = ", yf, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;
	packedsyms = cppackedsyms.clone();

	//uplo = uplo
	//N = K
	//alpha = alpha
	//x = xf
	//incx = 1
	//y = yf
	//incy = 1
	//Ap = packedsyms

	jniCBLAS.sspr2(order, uplo, K, alpha, xf, 1, yf, 1, packedsyms);
        CprintMatrix("Result: A = ", packedsyms, 1, K*(K+1)/2);
	
	//dsymv
	System.out.println();
        System.out.println("dsymv: y = alpha * A * x + beta * y, where A is symmetric matrix");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
	CprintMatrix("A = ", symfulld, K, K);
	xd = cpxd.clone();
	yd = cpyd.clone();
	RprintMatrix("x = ", xd, K, 1);
	RprintMatrix("y = ", yd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alpha
	//a = symfulld
	//lda = K
	//x = xd
	//incx = 1
	//beta = beta
	//y = yd
	//incy = 1
	
	jniCBLAS.dsymv(order, uplo, K, alpha, symfulld, K, xd, 1, beta, yd, 1);
        RprintMatrix("Result: x = ", yd, K, 1);
	
	//dsbmv
	System.out.println();
        System.out.println("dsbmv: y = alpha * A * x + beta * y, where A is symmetric matrix");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
	CprintMatrix("A = ", symfulld, K, K);
	xd = cpxd.clone();
	yd = cpyd.clone();
	RprintMatrix("x = ", xd, K, 1);
	RprintMatrix("y = ", yd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//K = 1
	//alpha = alpha
	//a = symbandd
	//lda = K
	//x = xd
	//incx = 1
	//beta = beta
	//y = yd
	//incy = 1
	
	jniCBLAS.dsbmv(order, uplo, K, 1, alpha, symbandd, K, xd, 1, beta, yd, 1);
        RprintMatrix("Result: x = ", yd, K, 1);

	//dspmv
	System.out.println();
        System.out.println("dspmv: y = alpha * A * x + beta * y, where A is symmetric matrix");
	System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
	CprintMatrix("A = ", symfulld, K, K);
	xd = cpxd.clone();
	yd = cpyd.clone();
	RprintMatrix("x = ", xd, K, 1);
	RprintMatrix("y = ", yd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alpha
	//Ap = packedsymd
	//x = xd
	//incx = 1
	//beta = beta
	//y = yd
	//incy = 1

	jniCBLAS.dspmv(order, uplo, K, alpha, packedsymd, xd, 1, beta, yd, 1);
        RprintMatrix("Result: x = ", yd, K, 1);
	
	//dger
	System.out.println();
        System.out.println("dger: A = alpha * x * y**T + A");
	System.out.printf("alpha = %f\n", alpha);
        CprintMatrix("A = ", symfulld, K, K);
	xd = cpxd.clone();
	yd = cpyd.clone();
	RprintMatrix("x = ", xd, K, 1);
	RprintMatrix("y = ", yd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//M = K 
	//N = K
	//alpha = alpha
	//x = xd
	//incx = 1
	//y = yd
	//incy = 1
	//a = symfulld
	//lda = K

	jniCBLAS.dger(order, K, K, alpha, xd, 1, yd, 1, symfulld, K);
        CprintMatrix("Result: A = ", symfulld, K, K);
	
	//dsyr
	System.out.println();
        System.out.println("dsyr: A = alpha * x * x**T + A");
	System.out.printf("alpha = %f\n", alpha);
	symfulld = cpsymfulld.clone();
        CprintMatrix("A = ", symfulld, K, K);
	xd = cpxd.clone();
	RprintMatrix("x = ", xd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alpha
	//x = xd
	//incx = 1
	//a = symfulld
	//lda = K
	
	jniCBLAS.dsyr(order, uplo, K, alpha, xd, 1, symfulld, K);
        CprintMatrix("Result: A = ", symfulld, K, K);
	
	//dspr
	System.out.println();
        System.out.println("dspr: A = alpha * x * x**T + A");
	System.out.printf("alpha = %f\n", alpha);
	symfulld = cpsymfulld.clone();
        CprintMatrix("A = ", symfulld, K, K);
	xd = cpxd.clone();
	RprintMatrix("x = ", xd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;
	packedsymd = cppackedsymd.clone();

	//uplo = uplo
	//N = K
	//alpha = alpha
	//x = xd
	//incx = 1
	//Ap = packedsymd

	jniCBLAS.dspr(order, uplo, K, alpha, xd, 1, packedsymd);
        CprintMatrix("Result: A = ", packedsymd, 1, K*(K+1)/2);
	
	//dsyr2
	System.out.println();
        System.out.println("dsyr2: A = alpha * x * y**T + alpha * y * x**T + A");
	System.out.printf("alpha = %f\n", alpha);
	symfulld = cpsymfulld.clone();
        CprintMatrix("A = ", symfulld, K, K);
	xd = cpxd.clone();
	yd = cpyd.clone();
	RprintMatrix("x = ", xd, K, 1);
	RprintMatrix("y = ", yd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alpha
	//x = xd	
	//incx = 1
	//y = yd
	//incy = 1
	//a = symfulld
	//lda = K
	
	jniCBLAS.dsyr2(order, uplo, K, alpha, xd, 1, yd, 1, symfulld, K);
        CprintMatrix("Result: A = ", symfulld, K, K);
	
	//dspr2
	System.out.println();
        System.out.println("dspr2: A = alpha * x * y**T + alpha * y * x**T + A");
	System.out.printf("alpha = %f\n", alpha);
	symfulld = cpsymfulld.clone();
        CprintMatrix("A = ", symfulld, K, K);
	xd = cpxd.clone();
	yd = cpyd.clone();
	RprintMatrix("x = ", xd, K, 1);
	RprintMatrix("y = ", yd, K, 1);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;
	packedsymd = cppackedsymd.clone();

	//uplo = uplo
	//N = K
	//alpha = alphad
	//x = xd
	//incx = 1
	//y = yd
	//incy = 1
	//Ap = packedsymd

	jniCBLAS.dspr2(order, uplo, K, alphad, xd, 1, yd, 1, packedsymd);
        CprintMatrix("Result: A = ", packedsymd, 1, K*(K+1)/2);
	
	//chemv
	System.out.println();
        System.out.println("chemv: y = alpha * A * x + beta * y, where A is hermitian matrix");
	System.out.printf("alpha = %f + %fi\n", alphac[0],alphac[1]);
        System.out.printf("beta = %f + %fi\n", betac[0],betac[1]);
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	byc = cpbyc.clone();	
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);
	order = jniCBLAS.ORDER.RowMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alphac
	//a = hermc
	//lda = K
	//x = bxc
	//incx = 1
	//beta = beta
	//y = byc
	//incy = 1
	
	jniCBLAS.chemv(order, uplo, K, alphac, hermc, K, bxc, 1, betac, byc, 1);
        RprintMatrix("Result: y = ", byc, K, 2);
	
	//chbmv
	//case of Rowmajor....
	System.out.println();
        System.out.println("chbmv: y = alpha * A * x + beta * y, where A is hermitian band matrix");
	System.out.printf("alpha = %f + %fi\n", alphac[0],alphac[1]);
        System.out.printf("beta = %f + %fi\n", betac[0],betac[1]);
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	byc = cpbyc.clone();	
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//K = 1
	//alpha = alphac
	//a = hermbandc
	//lda = K
	//x = bxc
	//incx = 1
	//beta = betac
	//y = byc
	//incy = 1
	
	jniCBLAS.chbmv(order, uplo, K, 1, alphac, hermbandc, K, bxc, 1, betac, byc, 1);
        RprintMatrix("Result: y = ", byc, K, 2);
	
	//chpmv
	System.out.println();
        System.out.println("chpmv: y = alpha * A * x + beta * y, where A is hermitian matrix");
	System.out.printf("alpha = %f + %fi\n", alphac[0],alphac[1]);
        System.out.printf("beta = %f + %fi\n", betac[0],betac[1]);
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	byc = cpbyc.clone();	
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alphac
	//Ap = packedhermc
	//x = bxc
	//incx = 1
	//beta = betac
	//y = byc
	//incy = 1
	
	jniCBLAS.chpmv(order, uplo, K, alphac, packedhermc, bxc, 1, betac, byc, 1);
        RprintMatrix("Result: y = ", byc, K, 2);
	
	//cgeru
	System.out.println();
        System.out.println("cgeru: A = alpha * x * y**T + A");
	System.out.printf("alpha = %f + %fi\n", alphac[0],alphac[1]);
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	byc = cpbyc.clone();	
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);
	order = jniCBLAS.ORDER.RowMajor;

	//M = K
	//N = K
	//alpha = alphac
	//x = bxc
	//incx = 1
	//y = byc
	//incy = 1
	//a = hermc
	//lda = K
	 	
	jniCBLAS.cgeru(order, K, K, alphac, bxc, 1, byc, 1, hermc, K);
        RprintMatrix("Result: A = ", hermc, K, 2*K);
	
	//cgerc
	System.out.println();
        System.out.println("cgerc: A = alpha * x * y**H + A");
	System.out.printf("alpha = %f + %fi\n", alphac[0],alphac[1]);
	hermc = cphermc.clone();
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	byc = cpbyc.clone();	
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);
	order = jniCBLAS.ORDER.RowMajor;

	//M = K
	//N = K
	//alpha = alphac
	//x = bxc
	//incx = 1
	//y = byc
	//incy = 1
	//a = hermc
	//lda = K
 	
	jniCBLAS.cgerc(order, K, K, alphac, bxc, 1, byc, 1, hermc, K);
        RprintMatrix("Result: A = ", hermc, K, 2*K);
	
	//cher
	System.out.println();
        System.out.println("cher: A = alpha * x * x**H + A");
	System.out.printf("alpha = %f\n", alpha);
	hermc = cphermc.clone();
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	RprintMatrix("x = ", bxc, K, 2);
	order = jniCBLAS.ORDER.RowMajor;
	uplo = jniCBLAS.UPLO.Upper; 

	//uplo = uplo
	//N = K
	//alpha = alpha
	//x = bxc
	//a = hermc	
	//lda = K
	
	jniCBLAS.cher(order, uplo, K, alpha, bxc, 1, hermc, K);
        RprintMatrix("Result: A = ", hermc, K, 2*K);
	
	//chpr
	System.out.println();
        System.out.println("chpr: A = alpha * x * x**H + A");
	System.out.printf("alpha = %f\n", alpha);
	hermc = cphermc.clone();
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	RprintMatrix("x = ", bxc, K, 2);
	order = jniCBLAS.ORDER.ColMajor;
	uplo = jniCBLAS.UPLO.Upper; 	

	//chpr = uplo
	//N = K
	//alpha = alpha	
	//x = bxc
	//incx = 1
	//Ap = packedhermc

	jniCBLAS.chpr(order, uplo, K, alpha, bxc, 1, packedhermc);
        RprintMatrix("Result: A = ", packedhermc, 1, K*(K+1));
	
	//cher2
	System.out.println();
        System.out.println("cher2: A = alpha * x * y**H + conjg(alpha) * y * x**H + A");
	System.out.printf("alpha = %f + %fi\n", alphac[0],alphac[1]);
	hermc = cphermc.clone();
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	byc = cpbyc.clone();
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);
	order = jniCBLAS.ORDER.RowMajor;
	uplo = jniCBLAS.UPLO.Upper; 

	//uplo = uplo
	//N = K
	//alhpa = alphac
	//x = bxc
	//incx = 1
	//y = byc
	//incy = 1
	//a = hermc
	//lda = K
		
	jniCBLAS.cher2(order, uplo, K, alphac, bxc, 1, byc, 1, hermc, K);
        RprintMatrix("Result: A = ", hermc, K, 2*K);
	
	//chpr2
	System.out.println();
        System.out.println("chpr2: A = alpha * x * y**H + conjg(alpha) * y * x**H + A");
	System.out.printf("alpha = %f + %fi\n", alphac[0],alphac[1]);
	hermc = cphermc.clone();
        RprintMatrix("A = ", hermc, K, 2*K);
	bxc = cpbxc.clone();
	byc = cpbyc.clone();
	RprintMatrix("x = ", bxc, K, 2);
	RprintMatrix("y = ", byc, K, 2);
	order = jniCBLAS.ORDER.ColMajor;
	uplo = jniCBLAS.UPLO.Upper; 
	packedhermc = cppackedhermc.clone();

	//uplo = uplo
	//N = K
	//alpha = alphac
	//x = bxc
	//incx = 1
	//y = byc
	//incy = 1
	//Ap = packedhermc
	
	jniCBLAS.chpr2(order, uplo, K, alphac, bxc, 1, byc, 1, packedhermc);
        RprintMatrix("Result: A = ", packedhermc, 1, K*(K+1));
	
	//zhemv
	System.out.println();
        System.out.println("zhemv: y = alpha * A * x + beta * y, where A is hermitian matrix");
	System.out.printf("alpha = %f + %fi\n", alphaz[0],alphaz[1]);
        System.out.printf("beta = %f + %fi\n", betaz[0],betaz[1]);
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	byz = cpbyz.clone();	
	RprintMatrix("x = ", bxz, K, 2);
	RprintMatrix("y = ", byz, K, 2);
	order = jniCBLAS.ORDER.RowMajor; 	
	uplo = jniCBLAS.UPLO.Upper;	

	//uplo = uplo
	//N = K
	//alpha = alphaz
	//a = hermz
	//lda = K
	//x = bxz
	//incx = 1
	//beta = betaz
	//y = byz
	//incy = 1
	
	jniCBLAS.zhemv(order, uplo, K, alphaz, hermz, K, bxz, 1, betaz, byz, 1);
        RprintMatrix("Result: y = ", byz, K, 2);
	
	//zhbmv
	System.out.println();
        System.out.println("zhbmv: y = alpha * A * x + beta * y, where A is hermitian band matrix");
	System.out.printf("alpha = %f + %fi\n", alphaz[0],alphaz[1]);
        System.out.printf("beta = %f + %fi\n", betaz[0],betaz[1]);
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	byz = cpbyz.clone();	
	RprintMatrix("x = ", bxz, K, 2);
	RprintMatrix("y = ", byz, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;	

	//uplo = uplo
	//N = K
	//K = 1
	//alpha = alphaz
	//a = hermbandz
	//lda = K
	//x = bxz
	//incx = 1
	//beta = betaz
	//y = byz
	//incy = 1
	
	jniCBLAS.zhbmv(order, uplo, K, 1, alphaz, hermbandz, K, bxz, 1, betaz, byz, 1);
        RprintMatrix("Result: y = ", byc, K, 2);
	
	//zhpmv
	System.out.println();
        System.out.println("zhpmv: y = alpha * A * x + beta * y, where A is hermitian matrix");
	System.out.printf("alpha = %f + %fi\n", alphaz[0],alphaz[1]);
        System.out.printf("beta = %f + %fi\n", betaz[0],betaz[1]);
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	byz = cpbyz.clone();	
	RprintMatrix("x = ", bxz, K, 2);
	RprintMatrix("y = ", byz, K, 2);
	order = jniCBLAS.ORDER.ColMajor; 	
	uplo = jniCBLAS.UPLO.Upper;	

	//uplo = uplo
	//N = K
	//alpha = alphaz
	//Ap = packedhermz
	//x = bxz
	//incx = 1
	//beta = betaz
	//y = byz
	//incy = 1
	
	jniCBLAS.zhpmv(order, uplo, K, alphaz, packedhermz, bxz, 1, betaz, byz, 1);
        RprintMatrix("Result: y = ", byz, K, 2);
	
	//zgeru
	System.out.println();
        System.out.println("zgeru: A = alpha * x * y**T + A");
	System.out.printf("alpha = %f + %fi\n", alphaz[0],alphaz[1]);
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	byz = cpbyz.clone();	
	RprintMatrix("x = ", bxz, K, 2);
	RprintMatrix("y = ", byz, K, 2);
	order = jniCBLAS.ORDER.RowMajor; 

	//M = K
	//N = K
	//alpha = alphaz
	//x = bxz
	//incx = 1
	//y = byz
	//incy = 1
	//a = hermz
	//lda = K
		
	jniCBLAS.zgeru(order, K, K, alphaz, bxz, 1, byz, 1, hermz, K);
        RprintMatrix("Result: A = ", hermz, K, 2*K);
	
	//zgerc
	System.out.println();
        System.out.println("zgerc: A = alpha * x * y**H + A");
	System.out.printf("alpha = %f + %fi\n", alphaz[0],alphaz[1]);
	hermz = cphermz.clone();
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	byz = cpbyz.clone();	
	RprintMatrix("x = ", bxz, K, 2);
	RprintMatrix("y = ", byz, K, 2);
	order = jniCBLAS.ORDER.RowMajor; 

	//M = K
	//N = K
	//alpha = alphaz
	//x = bxz
	//incx = 1
	//y = byz
	//incy = 1
	//a = hermz
	//lda = K
 		
	jniCBLAS.zgerc(order, K, K, alphaz, bxz, 1, byz, 1, hermz, K);
        RprintMatrix("Result: A = ", hermz, K, 2*K);
	
	//zher
	System.out.println();
        System.out.println("zher: A = alpha * x * x**H + A");
	System.out.printf("alpha = %f\n", alphad);
	hermz = cphermz.clone();
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	RprintMatrix("x = ", bxz, K, 2);
	order = jniCBLAS.ORDER.RowMajor;
	uplo = jniCBLAS.UPLO.Upper;

	//uplo = uplo
	//N = K
	//alpha = alphad
	//x = bxz
	//a = hermz	
	//lda = K
	 	
	jniCBLAS.zher(order, uplo, K, alphad, bxz, 1, hermz, K);
        RprintMatrix("Result: A = ", hermz, K, 2*K);
	
	//zhpr
	System.out.println();
        System.out.println("zhpr: A = alpha * x * x**H + A");
	System.out.printf("alpha = %f\n", alphad);
	hermz = cphermz.clone();
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	RprintMatrix("x = ", bxz, K, 2);
	order = jniCBLAS.ORDER.ColMajor;
	uplo = jniCBLAS.UPLO.Upper; 	

	//chpr = uplo
	//N = K
	//alpha = alphad
	//x = bxz
	//incx = 1
	//Ap = packedhermz

	jniCBLAS.zhpr(order, uplo, K, alphad, bxz, 1, packedhermz);
        RprintMatrix("Result: A = ", packedhermz, 1, K*(K+1));
	
	//zher2
	System.out.println();
        System.out.println("zher2: A = alpha * x * y**H + conjg(alpha) * y * x**H + A");
	System.out.printf("alpha = %f + %fi\n", alphaz[0],alphaz[1]);
	hermz = cphermz.clone();
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	byz = cpbyz.clone();
	RprintMatrix("x = ", bxz, K, 2);
	RprintMatrix("y = ", byz, K, 2);
	order = jniCBLAS.ORDER.RowMajor;
	uplo = jniCBLAS.UPLO.Upper; 	

	//uplo = uplo
	//N = K
	//alhpa = alphaz
	//x = bxz
	//incx = 1
	//y = byz
	//incy = 1
	//a = hermz
	//lda = K
	
	jniCBLAS.zher2(order, uplo, K, alphaz, bxz, 1, byz, 1, hermz, K);
        RprintMatrix("Result: A = ", hermz, K, 2*K);
	
	//zhpr2
	System.out.println();
        System.out.println("zhpr2: A = alpha * x * y**H + conjg(alpha) * y * x**H + A");
	System.out.printf("alpha = %f + %fi\n", alphaz[0],alphaz[1]);
	hermz = cphermz.clone();
        RprintMatrix("A = ", hermz, K, 2*K);
	bxz = cpbxz.clone();
	byz = cpbyz.clone();
	RprintMatrix("x = ", bxz, K, 2);
	RprintMatrix("y = ", byz, K, 2);
	order = jniCBLAS.ORDER.ColMajor;
	uplo = jniCBLAS.UPLO.Upper; 
	packedhermz = cppackedhermz.clone();

	//uplo = uplo
	//N = K
	//alpha = alphaz
	//x = bxz
	//incx = 1
	//y = byz
	//incy = 1
	//Ap = packedhermz
		
	jniCBLAS.zhpr2(order, uplo, K, alphaz, bxz, 1, byz, 1, packedhermz);
        RprintMatrix("Result: A = ", packedhermz, 1, K*(K+1));
	
	//sgemm
	System.out.println();
        System.out.println("sgemm: C = alpha*op( A )*op( B ) + beta*C ");
        System.out.println("transA = N, transB = N, order = C");
	System.out.println("alpha = " + string(alpha));
	System.out.println("beta = " + string(beta));
	CprintMatrix("Matrix A", A, M, K);
        CprintMatrix("Matrix B", Bf, K, N);
        CprintMatrix("Matrix C", Cf, M, N);
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	transB = jniCBLAS.TRANSPOSE.NoTrans;
	order = jniCBLAS.ORDER.ColMajor;

	//TransA = transA
	//TransB = transB
	//M = M
	//N = N
	//K = K
	//alpha = alpha	
	//a = A
	//lda = M
	//b = Bf
	//ldb = K
	//beta = beta
	//c = Cf
	//ldc = M

        jniCBLAS.sgemm(order, transA, transB, M, N, K, alpha, A, M, Bf, K, beta, Cf, M);
        CprintMatrix("Resulting C", Cf, M, N);

	//ssymm
	System.out.println();
        System.out.println("ssymm: C = alpha * A * B + beta * C, or C = alpha * B * A + beta * C ");
        System.out.println("alpha = " + string(alpha));
	System.out.println("beta = " + string(beta));
	CprintMatrix("Matrix A", symfulls, K, K);
        CprintMatrix("Matrix B", cpA, K, M);
        CprintMatrix("Matrix C", A, K, M);
	order = jniCBLAS.ORDER.ColMajor;
	side = jniCBLAS.SIDE.Left;
	uplo = jniCBLAS.UPLO.Upper;

	//side = side
	//uplo = uplo
	//M = K
	//N = M
	//alpha = alpha
	//a = symfulls
	//lda = K
	//b = cpA
	//ldb = K
	//beta = beta
	//c = A
	//ldc = K
	
        jniCBLAS.ssymm(order, side, uplo, K, M, alpha, symfulls, K, cpA, K, beta, A, K);
        CprintMatrix("Resulting C", A, K, M);

	//ssyrk
	System.out.println();
        System.out.println("ssyrk: C = alpha * A * A**T + beta * C, or C = alpha * A**T * A + beta * C ");
        System.out.println("alpha = " + string(alpha));
	System.out.println("beta = " + string(beta));
	A = cpA.clone();
	CprintMatrix("Matrix A", A, K, M);
        CprintMatrix("Matrix C", symfulls, K, K);
	order = jniCBLAS.ORDER.ColMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;

	//Uplo = uplo
	//Trans = transA
	//N = K
	//K = M
	//alpha = alpha
	//a = A
	//lda = K
	//beta = beta
	//c = symfulls
	//ldc = K
	
        jniCBLAS.ssyrk(order, uplo, transA, K, M, alpha, A, K, beta, symfulls, K);
        CprintMatrix("Resulting C", symfulls, K, K);

        //ssyr2k
	System.out.println();
        System.out.println("ssyr2k: C = alpha * A * B**T + alpha * B * A**T + beta * C, or C = alpha * A**T * B + alpha * B**T * A + beta * C ");
        System.out.println("alpha = " + string(alpha));
	System.out.println("beta = " + string(beta));
	A = cpA.clone();
	CprintMatrix("Matrix A", myaf, M, K);
        CprintMatrix("Matrix B", mybf, M, K);
        CprintMatrix("Matrix C", mycf, M, M);
	order = jniCBLAS.ORDER.ColMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;

	//Uplo = uplo
	//Trans = transA
	//N = M
	//K = K
	//alpha = alpha
	//a = myaf
	//lda = M
	//b = mybf
	//ldb = M
	//beta = beta
	//c = mycf
	//ldc = M

        jniCBLAS.ssyr2k(order, uplo, transA, M, K, alpha, myaf, M, mybf, M, beta, mycf, M);
        CprintMatrix("Resulting C", mycf, M, M);

	//strmm
	System.out.println();
        System.out.println("strmm: B = alpha * op(A) * B, or B = alpha * B * op(A), where A is a triangular matrix");
        System.out.println("alpha = " + string(alpha));
	CprintMatrix("Matrix A", trif, K, K);
        CprintMatrix("Matrix B", myaf, K, M);
        order = jniCBLAS.ORDER.ColMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	side = jniCBLAS.SIDE.Left;
	diag = jniCBLAS.DIAG.NonUnit;

	//Side = side
	//Uplo = uplo
	//Trans = transA
	//Diag = diag
	//M = K
	//N = M
	//alpha = alpha
	//a = trif
	//lda = K
	//b = myaf
	//ldb = K

        jniCBLAS.strmm(order, side, uplo, transA, diag, K, M, alpha, trif, K, myaf, K);
        CprintMatrix("Resulting C", myaf, K, M);

	//strsm
	System.out.println();
        System.out.println("strsm: op(A) * X = alpha * B, or X * op(A) = alpha * B, where A is a triangular matrix");
        System.out.println("alpha = " + string(alpha));
	CprintMatrix("Matrix A", solvetrif, K, K);
	myaf = cpmyaf.clone();
        CprintMatrix("Matrix B", myaf, K, M);
        order = jniCBLAS.ORDER.ColMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	side = jniCBLAS.SIDE.Left;
	diag = jniCBLAS.DIAG.NonUnit;

	//Side = side
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//M = K
	//N = M
	//alpha = alpha
	//a = solvetrif
	//lda = K
	//b = myaf
	//ldb = K

        jniCBLAS.strsm(order, side, uplo, transA, diag, K, M, alpha, solvetrif, K, myaf, K);
        CprintMatrix("Solution X = ", myaf, K, M);

	//dgemm
	System.out.println();
        System.out.println("dgemm: C = alpha*op( A )*op( B ) + beta*C ");
        System.out.println("transA = N, transB = N, order = C");
	System.out.println("alpha = " + string(alpha));
	System.out.println("beta = " + string(beta));
	CprintMatrix("Matrix A", dA, M, K);
        CprintMatrix("Matrix B", Bd, K, N);
        CprintMatrix("Matrix C", Cd, M, N);
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	transB = jniCBLAS.TRANSPOSE.NoTrans;
	order = jniCBLAS.ORDER.ColMajor;

	//TransA = transA
	//TransB = transB
	//M = M
	//N = N
	//K = K
	//alpha = alphad	
	//a = dA
	//lda = M
	//b = Bd
	//ldb = K
	//beta = betad
	//c = Cd
	//ldc = M

        jniCBLAS.dgemm(order, transA, transB, M, N, K, alphad, dA, M, Bd, K, betad, Cd, M);
        CprintMatrix("Resulting C", Cd, M, N);

	//dsymm
	System.out.println();
        System.out.println("dsymm: C = alpha * A * B + beta * C, or C = alpha * B * A + beta * C ");
        System.out.println("alpha = " + string(alpha));
	System.out.println("beta = " + string(beta));
	CprintMatrix("Matrix A", symfulld, K, K);
        CprintMatrix("Matrix B", cpdA, K, M);
        CprintMatrix("Matrix C", dA, K, M);
	order = jniCBLAS.ORDER.ColMajor;
	side = jniCBLAS.SIDE.Left;
	uplo = jniCBLAS.UPLO.Upper;

	//side = side
	//uplo = uplo
	//M = K
	//N = M
	//alpha = alphad
	//a = symfulld
	//lda = K
	//b = cpdA
	//ldb = K
	//beta = betad
	//c = dA
	//ldc = K
	
        jniCBLAS.dsymm(order, side, uplo, K, M, alphad, symfulld, K, cpdA, K, betad, dA, K);
        CprintMatrix("Resulting C", dA, K, M);

	//dsyrk
	System.out.println();
        System.out.println("dsyrk: C = alpha * A * A**T + beta * C, or C = alpha * A**T * A + beta * C ");
        System.out.println("alpha = " + string(alpha));
	System.out.println("beta = " + string(beta));
	dA = cpdA.clone();
	CprintMatrix("Matrix A", dA, K, M);
        CprintMatrix("Matrix C", symfulld, K, K);
	order = jniCBLAS.ORDER.ColMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;

	//Uplo = uplo
	//Trans = transA
	//N = K
	//K = M
	//alpha = alphad
	//a = dA
	//lda = K
	//beta = betad
	//c = symfulld
	//ldc = K
	
        jniCBLAS.dsyrk(order, uplo, transA, K, M, alphad, dA, K, betad, symfulld, K);
        CprintMatrix("Resulting C", symfulld, K, K);

	//dsyr2k
	System.out.println();
        System.out.println("dsyr2k: C = alpha * A * B**T + alpha * B * A**T + beta * C, or C = alpha * A**T * B + alpha * B**T * A + beta * C ");
        System.out.println("alpha = " + string(alpha));
	System.out.println("beta = " + string(beta));
	CprintMatrix("Matrix A", myad, M, K);
        CprintMatrix("Matrix B", mybd, M, K);
        CprintMatrix("Matrix C", mycd, M, M);
	order = jniCBLAS.ORDER.ColMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;

	//Uplo = uplo
	//Trans = transA
	//N = M
	//K = K
	//alpha = alphad
	//a = myad
	//lda = M
	//b = mybd
	//ldb = M
	//beta = betad
	//c = mycd
	//ldc = M

        jniCBLAS.dsyr2k(order, uplo, transA, M, K, alphad, myad, M, mybd, M, betad, mycd, M);
        CprintMatrix("Resulting C", mycd, M, M);

	//dtrmm
	System.out.println();
        System.out.println("dtrmm: B = alpha * op(A) * B, or B = alpha * B * op(A), where A is a triangular matrix");
        System.out.println("alpha = " + string(alpha));
	CprintMatrix("Matrix A", trid, K, K);
        CprintMatrix("Matrix B", myad, K, M);
        order = jniCBLAS.ORDER.ColMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	side = jniCBLAS.SIDE.Left;
	diag = jniCBLAS.DIAG.NonUnit;

	//Side = side
	//Uplo = uplo
	//Trans = transA
	//Diag = diag
	//M = K
	//N = M
	//alpha = alphad
	//a = trid
	//lda = K
	//b = myad
	//ldb = K

        jniCBLAS.dtrmm(order, side, uplo, transA, diag, K, M, alphad, trid, K, myad, K);
        CprintMatrix("Resulting C", myad, K, M);

	//dtrsm
	System.out.println();
        System.out.println("dtrsm: op(A) * X = alpha * B, or X * op(A) = alpha * B, where A is a triangular matrix");
        System.out.println("alpha = " + string(alpha));
	CprintMatrix("Matrix A", solvetrid, K, K);
	myad = cpmyad.clone();
        CprintMatrix("Matrix B", myad, K, M);
        order = jniCBLAS.ORDER.ColMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	side = jniCBLAS.SIDE.Left;
	diag = jniCBLAS.DIAG.NonUnit;

	//Side = side
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//M = K
	//N = M
	//alpha = alphad
	//a = solvetrid
	//lda = K
	//b = myad
	//ldb = K

        jniCBLAS.dtrsm(order, side, uplo, transA, diag, K, M, alphad, solvetrid, K, myad, K);
        CprintMatrix("Solution X = ", myad, K, M);

	//cgemm
	System.out.println();
        System.out.println("cgemm: C = alpha*op( A )*op( B ) + beta*C ");
        System.out.println("transA = N, transB = N, order = C");
	System.out.printf("alpha = %f + %fi\n",alphac[0],alphac[1]);
	System.out.printf("beta = %f + %fi\n",betac[0],betac[1]);
	RprintMatrix("Matrix A", Ac, M, 2*K);
        RprintMatrix("Matrix B", Bc, K, 2*N);
        RprintMatrix("Matrix C", Cc, M, 2*N);
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	transB = jniCBLAS.TRANSPOSE.NoTrans;
	order = jniCBLAS.ORDER.RowMajor;

	//TransA = transA
	//TransB = transB
	//M = M
	//N = N
	//K = K
	//alpha = alphac
	//a = Ac
	//lda = K
	//b = Bc
	//ldb = N
	//beta = betac
	//c = Cc
	//ldc = N

        jniCBLAS.cgemm(order, transA, transB, M, N, K, alphac, Ac, K, Bc, N, betac, Cc, N);
        RprintMatrix("Resulting C", Cc, M, 2*N);

	//csymm
	System.out.println();
        System.out.println("csymm: C = alpha * A * B + beta * C, or C = alpha * B * A + beta * C ");
        System.out.printf("alpha = %f + %fi\n",alphac[0],alphac[1]);
	System.out.printf("beta = %f + %fi\n",betac[0],betac[1]);
	RprintMatrix("Matrix A", symfullc, K, 2*K);
        RprintMatrix("Matrix B", cpAc, K, 2*M);
        RprintMatrix("Matrix C", Ac, K, 2*M);
	order = jniCBLAS.ORDER.RowMajor;
	side = jniCBLAS.SIDE.Left;
	uplo = jniCBLAS.UPLO.Upper;

	//side = side
	//uplo = uplo
	//M = K
	//N = M
	//alpha = alphac
	//a = symfullc
	//lda = K
	//b = cpAc
	//ldb = M
	//beta = betac
	//c = Ac
	//ldc = M
	
        jniCBLAS.csymm(order, side, uplo, K, M, alphac, symfullc, K, cpAc, M, betac, Ac, M);
        RprintMatrix("Resulting C", Ac, K, 2*M);

	//csyrk
	System.out.println();
        System.out.println("csyrk: C = alpha * A * A**T + beta * C, or C = alpha * A**T * A + beta * C ");
        System.out.printf("alpha = %f + %fi\n",alphac[0],alphac[1]);
	System.out.printf("beta = %f + %fi\n",betac[0],betac[1]);
	Ac = cpAc.clone();
	RprintMatrix("Matrix A", Ac, K, 2*M);
        RprintMatrix("Matrix C", symfullc, K, 2*K);
	order = jniCBLAS.ORDER.RowMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;

	//Uplo = uplo
	//Trans = transA
	//N = K
	//K = M
	//alpha = alphac
	//a = Ac
	//lda = M
	//beta = betac
	//c = symfullc
	//ldc = K
	
        jniCBLAS.csyrk(order, uplo, transA, K, M, alphac, Ac, M, betac, symfullc, K);
        RprintMatrix("Resulting C", symfullc, K, 2*K);

	//csyr2k
	System.out.println();
        System.out.println("csyr2k: C = alpha * A * B**T + alpha * B * A**T + beta * C, or C = alpha * A**T * B + alpha * B**T * A + beta * C ");
        System.out.printf("alpha = %f + %fi\n",alphac[0],alphac[1]);
	System.out.printf("beta = %f + %fi\n",betac[0],betac[1]);
	Ac = cpAc.clone();
	RprintMatrix("Matrix A", myac, M, 2*K);
        RprintMatrix("Matrix B", mybc, M, 2*K);
        RprintMatrix("Matrix C", mycc, M, 2*M);
	order = jniCBLAS.ORDER.RowMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;

	//Uplo = uplo
	//Trans = transA
	//N = M
	//K = K
	//alpha = alphac
	//a = myac
	//lda = K
	//b = mybc
	//ldb = K
	//beta = betac
	//c = mycc
	//ldc = M

        jniCBLAS.csyr2k(order, uplo, transA, M, K, alphac, myac, K, mybc, K, betac, mycc, M);
        RprintMatrix("Resulting C", mycc, M, 2*M);

	//ctrmm
	System.out.println();
        System.out.println("ctrmm: B = alpha * op(A) * B, or B = alpha * B * op(A), where A is a triangular matrix");
        System.out.printf("alpha = %f + %fi\n",alphac[0],alphac[1]);
	RprintMatrix("Matrix A", tric, K, 2*K);
        RprintMatrix("Matrix B", myac, K, 2*M);
        order = jniCBLAS.ORDER.RowMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	side = jniCBLAS.SIDE.Left;
	diag = jniCBLAS.DIAG.NonUnit;

	//Side = side
	//Uplo = uplo
	//Trans = transA
	//Diag = diag
	//M = K
	//N = M
	//alpha = alphac
	//a = tric
	//lda = K
	//b = myac
	//ldb = M

        jniCBLAS.ctrmm(order, side, uplo, transA, diag, K, M, alphac, tric, K, myac, M);
        RprintMatrix("Resulting B", myac, K, 2*M);

	//ctrsm
	System.out.println();
        System.out.println("ctrsm: op(A) * X = alpha * B, or X * op(A) = alpha * B, where A is a triangular matrix");
        System.out.printf("alpha = %f + %fi\n",alphac[0],alphac[1]);
	RprintMatrix("Matrix A", solvetric, K, 2*K);
	myac = cpmyac.clone();
        RprintMatrix("Matrix B", myac, K, 2*M);
        order = jniCBLAS.ORDER.RowMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	side = jniCBLAS.SIDE.Left;
	diag = jniCBLAS.DIAG.NonUnit;

	//Side = side
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//M = K
	//N = M
	//alpha = alphac
	//a = solvetric
	//lda = K
	//b = myac
	//ldb = M

        jniCBLAS.ctrsm(order, side, uplo, transA, diag, K, M, alphac, solvetric, K, myac, M);
        RprintMatrix("Solution X = ", myac, K, 2*M);

	//zgemm
	System.out.println();
        System.out.println("zgemm: C = alpha*op( A )*op( B ) + beta*C ");
        System.out.println("transA = N, transB = N, order = C");
	System.out.printf("alpha = %f + %fi\n",alphaz[0],alphaz[1]);
	System.out.printf("beta = %f + %fi\n",betaz[0],betaz[1]);
	RprintMatrix("Matrix A", Az, M, 2*K);
        RprintMatrix("Matrix B", Bz, K, 2*N);
        RprintMatrix("Matrix C", Cz, M, 2*N);
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	transB = jniCBLAS.TRANSPOSE.NoTrans;
	order = jniCBLAS.ORDER.RowMajor;

	//TransA = transA
	//TransB = transB
	//M = M
	//N = N
	//K = K
	//alpha = alphaz	
	//a = Az
	//lda = K
	//b = Bz
	//ldb = N
	//beta = betaz
	//c = Cz
	//ldc = N

        jniCBLAS.zgemm(order, transA, transB, M, N, K, alphaz, Az, K, Bz, N, betaz, Cz, N);
        RprintMatrix("Resulting C", Cz, M, 2*N);

	//zsymm
	System.out.println();
        System.out.println("zsymm: C = alpha * A * B + beta * C, or C = alpha * B * A + beta * C ");
        System.out.printf("alpha = %f + %fi\n",alphaz[0],alphaz[1]);
	System.out.printf("beta = %f + %fi\n",betaz[0],betaz[1]);
	RprintMatrix("Matrix A", symfullz, K, 2*K);
        RprintMatrix("Matrix B", cpAz, K, 2*M);
        RprintMatrix("Matrix C", Az, K, 2*M);
	order = jniCBLAS.ORDER.RowMajor;
	side = jniCBLAS.SIDE.Left;
	uplo = jniCBLAS.UPLO.Upper;

	//side = side
	//uplo = uplo
	//M = K
	//N = M
	//alpha = alphaz
	//a = symfullz
	//lda = K
	//b = cpAz
	//ldb = M
	//beta = betaz
	//c = Az
	//ldc = M
	
        jniCBLAS.zsymm(order, side, uplo, K, M, alphaz, symfullz, K, cpAz, M, betaz, Az, M);
        RprintMatrix("Resulting C", Az, K, 2*M);

	//zsyrk
	System.out.println();
        System.out.println("zsyrk: C = alpha * A * A**T + beta * C, or C = alpha * A**T * A + beta * C ");
        System.out.printf("alpha = %f + %fi\n",alphaz[0],alphaz[1]);
	System.out.printf("beta = %f + %fi\n",betaz[0],betaz[1]);
	Az = cpAz.clone();
	RprintMatrix("Matrix A", Az, K, 2*M);
        RprintMatrix("Matrix C", symfullz, K, 2*K);
	order = jniCBLAS.ORDER.RowMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;

	//Uplo = uplo
	//Trans = transA
	//N = K
	//K = M
	//alpha = alphaz
	//a = Az
	//lda = M
	//beta = betaz
	//c = symfullz
	//ldc = K
	
        jniCBLAS.zsyrk(order, uplo, transA, K, M, alphaz, Az, M, betaz, symfullz, K);
        RprintMatrix("Resulting C", symfullz, K, 2*K);

	//zsyr2k
	System.out.println();
        System.out.println("zsyr2k: C = alpha * A * B**T + alpha * B * A**T + beta * C, or C = alpha * A**T * B + alpha * B**T * A + beta * C ");
        System.out.printf("alpha = %f + %fi\n",alphaz[0],alphaz[1]);
	System.out.printf("beta = %f + %fi\n",betaz[0],betaz[1]);
	Az = cpAz.clone();
	RprintMatrix("Matrix A", myaz, M, 2*K);
        RprintMatrix("Matrix B", mybz, M, 2*K);
        RprintMatrix("Matrix C", mycz, M, 2*M);
	order = jniCBLAS.ORDER.RowMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;

	//Uplo = uplo
	//Trans = transA
	//N = M
	//K = K
	//alpha = alphaz
	//a = myaz
	//lda = K
	//b = mybz
	//ldb = K
	//beta = betaz
	//c = mycz
	//ldc = M

        jniCBLAS.zsyr2k(order, uplo, transA, M, K, alphaz, myaz, K, mybz, K, betaz, mycz, M);
        RprintMatrix("Resulting C", mycz, M, 2*M);

	//ztrmm
	System.out.println();
        System.out.println("ztrmm: B = alpha * op(A) * B, or B = alpha * B * op(A), where A is a triangular matrix");
        System.out.printf("alpha = %f + %fi\n",alphaz[0],alphaz[1]);
	RprintMatrix("Matrix A", triz, K, 2*K);
        RprintMatrix("Matrix B", myaz, K, 2*M);
        order = jniCBLAS.ORDER.RowMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	side = jniCBLAS.SIDE.Left;
	diag = jniCBLAS.DIAG.NonUnit;

	//Side = side
	//Uplo = uplo
	//Trans = transA
	//Diag = diag
	//M = K
	//N = M
	//alpha = alphaz
	//a = triz
	//lda = K
	//b = myaz
	//ldb = K

        jniCBLAS.ztrmm(order, side, uplo, transA, diag, K, M, alphaz, triz, K, myaz, M);
        RprintMatrix("Resulting B", myaz, K, 2*M);

	//ztrsm
	System.out.println();
        System.out.println("ztrsm: op(A) * X = alpha * B, or X * op(A) = alpha * B, where A is a triangular matrix");
        System.out.printf("alpha = %f + %fi\n",alphaz[0],alphaz[1]);
	RprintMatrix("Matrix A", solvetriz, K, 2*K);
	myaz = cpmyaz.clone();
        RprintMatrix("Matrix B", myaz, K, 2*M);
        order = jniCBLAS.ORDER.RowMajor;
	transA = jniCBLAS.TRANSPOSE.NoTrans;
	uplo = jniCBLAS.UPLO.Upper;
	side = jniCBLAS.SIDE.Left;
	diag = jniCBLAS.DIAG.NonUnit;

	//Side = side
	//Uplo = uplo
	//TransA = transA
	//Diag = diag
	//M = K
	//N = M
	//alpha = alphaz
	//a = solvetriz
	//lda = K
	//b = myaz
	//ldb = M

        jniCBLAS.ztrsm(order, side, uplo, transA, diag, K, M, alphaz, solvetriz, K, myaz, M);
        RprintMatrix("Solution X = ", myaz, K, 2*M);

	//chemm
	System.out.println();
        System.out.println("chemm: C = alpha*A*B + beta*C, or C = alpha*B*A + beta*C, where A is hermitian matrix");
        System.out.printf("alpha = %f + %fi\n",alphac[0],alphac[1]);
	System.out.printf("beta = %f + %fi\n",betac[0],betac[1]);
	hermc = cphermc.clone();
	RprintMatrix("Matrix A", hermc, K, 2*K);
        RprintMatrix("Matrix B", cpBc, K, 2*N);
        RprintMatrix("Matrix C", Bc, K, 2*N);
	order = jniCBLAS.ORDER.RowMajor;
	side = jniCBLAS.SIDE.Left;
	uplo = jniCBLAS.UPLO.Upper;

	//Side = side
	//Uplo = uplo
	//M = K
	//N = N
	//alpha = alphac
	//a = hermc
	//lda = K
	//b = cpBc
	//ldb = N
	//beta = betac
	//c = Bc
	//ldc = N

        jniCBLAS.chemm(order, side, uplo, K, N, alphac, hermc, K, cpBc, N, betac, Bc, N);
        RprintMatrix("Resulting C", Bc, K, 2*N);

	//cherk
	System.out.println();
        System.out.println("cherk: C = alpha*A*A**H + beta*C, or C = alpha*A**H*A + beta*C, where C is hermitian matrix");
        System.out.printf("alpha = %f\n",alpha);
	System.out.printf("beta = %f\n",beta);
	hermc = cphermc.clone();
	Bc = cpBc.clone();
	RprintMatrix("Matrix A", Bc, K, 2*N);
        RprintMatrix("Matrix C", hermc, K, 2*K);
	order = jniCBLAS.ORDER.RowMajor;
	uplo = jniCBLAS.UPLO.Upper;
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//Uplo = uplo
	//Trans = transA
	//N = K
	//K = N
	//alpha = alpha
	//a = Bc
	//lda = N
	//beta = beta
	//c = hermc
	//ldc = K
	
        jniCBLAS.cherk(order, uplo, transA, K, N, alpha, Bc, N, beta, hermc, K);
        RprintMatrix("Resulting C", hermc, K, 2*K);

	//cher2k
	System.out.println();
        System.out.println("cher2k: C = alpha*A*B**H + conjg(alpha)*B*A**H + beta*C, or C = alpha*A**H*B + conjg(alpha)*B**H*A + beta*C, where C is hermitian matrix");
        System.out.printf("alpha = %f + %fi\n",alphac[0],alphac[1]);
	System.out.printf("beta = %f\n",beta);
	hermc = cphermc.clone();
	Bc = cpBc.clone();
	RprintMatrix("Matrix A", Bc, K, 2*N);
        RprintMatrix("Matrix B", cpBc, K, 2*N);
        RprintMatrix("Matrix C", hermc, K, 2*K);
	order = jniCBLAS.ORDER.RowMajor;
	uplo = jniCBLAS.UPLO.Upper;
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//Uplo = uplo
	//Trans = transA
	//N = K
	//K = N
	//alpha = alphac
	//a = Bc
	//lda = N
	//b = cpBc
	//ldb = N
	//beta = beta
	//c = hermc
	//ldc = K

        jniCBLAS.cher2k(order, uplo, transA, K, N, alphac, Bc, N, cpBc, N, beta, hermc, K);
        RprintMatrix("Resulting C", hermc, K, 2*K);

	//zhemm
	System.out.println();
        System.out.println("zhemm: C = alpha*A*B + beta*C, or C = alpha*B*A + beta*C, where A is hermitian matrix");
        System.out.printf("alpha = %f + %fi\n",alphaz[0],alphaz[1]);
	System.out.printf("beta = %f + %fi\n",betaz[0],betaz[1]);
	hermz = cphermz.clone();
	RprintMatrix("Matrix A", hermz, K, 2*K);
        RprintMatrix("Matrix B", cpBz, K, 2*N);
        RprintMatrix("Matrix C", Bz, K, 2*N);
	order = jniCBLAS.ORDER.RowMajor;
	side = jniCBLAS.SIDE.Left;
	uplo = jniCBLAS.UPLO.Upper;

	//Side = side
	//Uplo = uplo
	//M = K
	//N = N
	//alpha = alphaz
	//a = hermz
	//lda = K
	//b = cpBz
	//ldb = N
	//beta = betaz	
	//c = Bz
	//ldc = N

        jniCBLAS.zhemm(order, side, uplo, K, N, alphaz, hermz, K, cpBz, N, betaz, Bz, N);
        RprintMatrix("Resulting C", Bz, K, 2*N);

	//zherk
	System.out.println();
        System.out.println("zherk: C = alpha*A*A**H + beta*C, or C = alpha*A**H*A + beta*C, where C is hermitian matrix");
        System.out.printf("alpha = %f\n",alpha);
	System.out.printf("beta = %f\n",beta);
	hermz = cphermz.clone();
	Bz = cpBz.clone();
	RprintMatrix("Matrix A", Bz, K, 2*N);
        RprintMatrix("Matrix C", hermz, K, 2*K);
	order = jniCBLAS.ORDER.RowMajor;
	uplo = jniCBLAS.UPLO.Upper;
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//Uplo = uplo
	//Trans = transA
	//N = K
	//K = N
	//alpha = alphad
	//a = Bz
	//lda = N
	//beta =  betad
	//c = hermz	
	//ldc = K

        jniCBLAS.zherk(order, uplo, transA, K, N, alphad, Bz, N, betad, hermz, K);
        RprintMatrix("Resulting C", hermz, K, 2*K);

	//zher2k
	System.out.println();
        System.out.println("zher2k: C = alpha*A*B**H + conjg(alpha)*B*A**H + beta*C, or C = alpha*A**H*B + conjg(alpha)*B**H*A + beta*C, where C is hermitian matrix");
        System.out.printf("alpha = %f + %fi\n",alphaz[0],alphaz[1]);
	System.out.printf("beta = %f\n",betad);
	hermz = cphermz.clone();
	Bz = cpBz.clone();
	RprintMatrix("Matrix A", Bz, K, 2*N);
        RprintMatrix("Matrix B", cpBz, K, 2*N);
        RprintMatrix("Matrix C", hermz, K, 2*K);
	order = jniCBLAS.ORDER.RowMajor;
	uplo = jniCBLAS.UPLO.Upper;
	transA = jniCBLAS.TRANSPOSE.NoTrans;

	//Uplo = uplo
	//Trans = transA
	//N = K
	//K = N
	//alpha = alphaz
	//a = Bz
	//lda = N
	//b = cpBz
	//ldb = N
	//beta = betad
	//c = hermz
	//ldc = K

        jniCBLAS.zher2k(order, uplo, transA, K, N, alphaz, Bz, N, cpBz, N, betad, hermz, K);
        RprintMatrix("Resulting C", hermz, K, 2*K);
 	

    }
    
    //print vector into matrix, assume vector is col-majored
    private static void CprintMatrix(String prompt, double[] X, int I, int J) {
        System.out.println(prompt);
        double[][] mat = new double[I][J]; 
        int count = 0;
        
        for (int j=0; j<J; j++){
          for (int i=0; i<I; i++){
            mat[i][j] = X[count];
            count = count+1;
          }
        }
          
        for (int i=0; i<I; i++) {
            for (int j=0; j<J; j++)
                //System.out.print("\t" + string(mat[i][j]));
                System.out.printf("%.4f\t", mat[i][j]);
            System.out.println();
        }
    }
    
    private static void CprintMatrix(String prompt, float[] X, int I, int J) {
        System.out.println(prompt);
        float[][] mat = new float[I][J]; 
        int count = 0;
        
        for (int j=0; j<J; j++){
          for (int i=0; i<I; i++){
            mat[i][j] = X[count];
            count = count+1;
          }
        }
          
        for (int i=0; i<I; i++) {
            for (int j=0; j<J; j++)
                //System.out.print("\t" + string(mat[i][j]));
                System.out.printf("%.4f\t", mat[i][j]);
            System.out.println();
        }
    }
    
    //print vector into matrix, assume vector is row-majored
    private static void RprintMatrix(String prompt, double[] X, int I, int J) {
        System.out.println(prompt);
        for (int i=0; i<I; i++) {
            for (int j=0; j<J; j++)
                System.out.printf("%.4f\t", X[i*J+j]);
            System.out.println();
        }
    }
    
    private static void RprintMatrix(String prompt, float[] X, int I, int J) {
        System.out.println(prompt);
        for (int i=0; i<I; i++) {
            for (int j=0; j<J; j++)
                System.out.printf("%.4f\t", X[i*J+j]);
            System.out.println();
        }
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



