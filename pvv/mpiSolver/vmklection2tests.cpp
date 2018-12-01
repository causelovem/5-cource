#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>        // openmp
#include <emmintrin.h>  // simd 

int main(int, char **){
 
#define NN 330000123 //123 - to check correctnes, so that NN%2, NN%4... !=0 

#if(1) // ADD explicit unroll
    {
        int N=NN; 
        double a = 1.0;
        double result, t;
        printf("Testing loop unrolling \n"); 
        t = omp_get_wtime(); // starting of wall clock time measurememt
        { 
            result = 0;
            for(int i=0; i<N; ++i){
                result += a;
            }
        }
        t = omp_get_wtime() - t;  // finishing timing
        printf("ADD X1 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));
 
        result=0;
        t = omp_get_wtime(); // starting of wall clock time measurememt
        {
            double result_tmp[2] = {0,0};
            for(int i=0; i<N/2; ++i){
                result_tmp[0] += a;
                result_tmp[1] += a;
            }
            for(int i=0; i<N%2; ++i) result += a; // finishing remaining iterations
            result += result_tmp[0] + result_tmp[1]; // reduction of result
        }
        t = omp_get_wtime() - t; // finishing timing
        printf("ADD X2 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=0;
        t = omp_get_wtime();
        {
            double result_tmp[3] = {0,0,0};
            for(int i=0; i<N/3; ++i){
                result_tmp[0] += a;
                result_tmp[1] += a;
                result_tmp[2] += a;
            }
            for(int i=0; i<N%3; ++i) result += a;
            result += result_tmp[0] + result_tmp[1] + result_tmp[1];
        }
        t = omp_get_wtime() - t;
        printf("ADD X3 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=0;
        t = omp_get_wtime();
        {
            double result_tmp[4] = {0,0,0,0};
            for(int i=0; i<N/4; ++i){
                result_tmp[0] += a;
                result_tmp[1] += a;
                result_tmp[2] += a;
                result_tmp[3] += a;
            }
            for(int i=0; i<N%4; ++i) result += a;
            result += result_tmp[0] + result_tmp[1] + result_tmp[2] + result_tmp[3];
        }
        t = omp_get_wtime() - t;
        printf("ADD X4 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=0;
        t = omp_get_wtime();
        {
            double result_tmp[5] = {0,0,0,0,0};
            for(int i=0; i<N/5; ++i){
                result_tmp[0] += a;
                result_tmp[1] += a;
                result_tmp[2] += a;
                result_tmp[3] += a;
                result_tmp[4] += a;
            }
            for(int i=0; i<N%5; ++i) result += a;
            result += result_tmp[0] + result_tmp[1] + result_tmp[2] + result_tmp[3] + result_tmp[4];
        }
        t = omp_get_wtime() - t;
        printf("ADD X5 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));
    }
#endif

#if(1) // MUL explicit unroll
    {
        int N=NN; 
        double a = 1.0000001;
        double result, t;
        printf("\n"); 
        t = omp_get_wtime();
        {
            result = 1.0;
            for(int i=0; i<N; ++i){
                result *= a;
            }
        }
        t = omp_get_wtime() - t;
        // note roundoff error may accumulate in last digits
        printf("MUL X1 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));
        
        result=1.0;
        t = omp_get_wtime();
        {
            double result_tmp[2]={1,1};
            for(int i=0; i<N/2; ++i){
                // using constant loop assuming compiler will unroll it 
                for(int j=0; j<2; ++j) result_tmp[j] *= a; 
            }
            for(int i=0; i<N%2; ++i) result *= a;
            for(int j=0; j<2;   ++j) result *= result_tmp[j];
        }
        t = omp_get_wtime() - t;
        printf("MUL X2 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=1.0;
        t = omp_get_wtime();
        {
            double result_tmp[3]={1,1,1};
            for(int i=0; i<N/3; ++i){
                for(int j=0; j<3; ++j) result_tmp[j] *= a;
            }
            for(int i=0; i<N%3; ++i) result *= a;
            for(int j=0; j<3;   ++j) result *= result_tmp[j];
        }
        t = omp_get_wtime() - t;
        printf("MUL X3 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=1.0;
        t = omp_get_wtime();
        {
            double result_tmp[4]={1,1,1,1};
            for(int i=0; i<N/4; ++i){
                for(int j=0; j<4; ++j) result_tmp[j] *= a;
            }
            for(int i=0; i<N%4; ++i) result *= a;
            for(int j=0; j<4;   ++j) result *= result_tmp[j];
        }
        t = omp_get_wtime() - t;
        printf("MUL X4 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=1.0;
        t = omp_get_wtime();
        {
            double result_tmp[5]={1,1,1,1,1};
            for(int i=0; i<N/5; ++i){
                for(int j=0; j<5; ++j) result_tmp[j] *= a;
            }
            for(int i=0; i<N%5; ++i) result *= a;
            for(int j=0; j<5;   ++j) result *= result_tmp[j];
        }
        t = omp_get_wtime() - t;
        printf("MUL X5 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=1.0;
        t = omp_get_wtime();       
        {
            double result_tmp[6]={1,1,1,1,1,1};
            for(int i=0; i<N/6; ++i){
                for(int j=0; j<6; ++j) result_tmp[j] *= a;
            }
            for(int i=0; i<N%6; ++i) result *= a;
            for(int j=0; j<6;   ++j) result *= result_tmp[j];
        }
        t = omp_get_wtime() - t;
        printf("MUL X6 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));
    }
#endif

#if(1) // VEC ADD explicit unroll (SIMD extension must be enabled in compiler settings)
    {
        printf("\n");
        printf("Testing simd add\n"); 
        int N=NN; 
        double a = 1.0;
        double result, t;

        result=0;
        {
            __m128d result_tmp = _mm_set_pd(0.0,0.0); // 128 bit vector time = 2 FP64 values 
            __m128d aa = _mm_set_pd(a,a);

            t = omp_get_wtime();
            for(int i=0; i<N/2; ++i){
                result_tmp = _mm_add_pd(result_tmp, aa); // vector addition intrinsic 
            }
            for(int i=0; i<N%2; ++i) result += a;   
            result += ((double*)(&result_tmp))[0]+((double*)(&result_tmp))[1]; // reduction 
            t = omp_get_wtime() - t;
        }
        printf("VECADD X1 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=0;
        {
            __m128d result_tmp[2] = {_mm_set_pd(0.0,0.0),_mm_set_pd(0.0,0.0)};
            __m128d aa = _mm_set_pd(a,a);

            t = omp_get_wtime();
            for(int i=0; i<N/4; ++i){
                result_tmp[0] = _mm_add_pd(result_tmp[0], aa);
                result_tmp[1] = _mm_add_pd(result_tmp[1], aa);
            }
            for(int i=0; i<N%4; ++i) result += a;   
            for(int j=0; j<2;   ++j) result += ((double*)(&result_tmp[j]))[0]+((double*)(&result_tmp[j]))[1];
            t = omp_get_wtime() - t;

        }
        printf("VECADD X2 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        result=0;
        {
            __m128d result_tmp[4];
            __m128d aa = _mm_set_pd(a,a);

            for(int i=0; i<4; ++i) result_tmp[i]=_mm_set_pd(0.0,0.0);

            t = omp_get_wtime();
            for(int i=0; i<N/8; ++i){
                result_tmp[0] = _mm_add_pd(result_tmp[0], aa);
                result_tmp[1] = _mm_add_pd(result_tmp[1], aa);
                result_tmp[2] = _mm_add_pd(result_tmp[2], aa);
                result_tmp[3] = _mm_add_pd(result_tmp[3], aa);
            }
            for(int i=0; i<N%8; ++i) result += a;   
            for(int j=0; j<4; ++j)   result += ((double*)(&result_tmp[j]))[0]+((double*)(&result_tmp[j]))[1];
            t = omp_get_wtime() - t;
        }
        printf("VECADD X4 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));
    }

#endif
 
#if(1)
    { // Elementwise summation of vectors X=X+Y
#define MM 66000123

        int N=MM; 
        printf("\n");
        double *x=new double[N];
        double *y=new double[N];

        double result, t, t1, t4;
        printf("Vector size %d bytes\n", N*(int)sizeof(double));

        for(int i=0; i<N; ++i){x[i]=i%2; y[i]=i%3; } // initialization with some crap
        result = 0.0;
        t = omp_get_wtime();
        {
            for(int i=0; i<N; ++i){
                x[i] = x[i]+y[i];
            }
        }
        t = omp_get_wtime() - t;
        for(int i=0; i<N; ++i) result += x[i];
        printf("XPY STRD1 T1 X1 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));
        t1 = t; 
 
        for(int i=0; i<N; ++i){x[i]=i%2; y[i]=i%3;}
        result = 0.0;
        t = omp_get_wtime();
        {
            const int N4 = (N/4)*4; 
            for(int i=0; i<N4; i+=4){
                x[i  ] = x[i]  +y[i];
                x[i+1] = x[i+1]+y[i+1];
                x[i+2] = x[i+2]+y[i+2];
                x[i+3] = x[i+3]+y[i+3];
            }
            for(int i=N4; i<N; ++i) x[i] = x[i] + y[i];
        }
        t = omp_get_wtime() - t;
        for(int i=0; i<N; ++i) result += x[i];
        printf("XPY STRD1 T1 X4 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        for(int i=0; i<N; ++i){x[i]=i%2; y[i]=i%3;}
        result = 0.0;
        t = omp_get_wtime();
        {
            #pragma omp parallel for num_threads(2)
            for(int i=0; i<N; ++i){
                x[i] = x[i]+y[i];
            }
        }
        t = omp_get_wtime() - t;
        for(int i=0; i<N; ++i) result += x[i];
        printf("XPY STRD1 T2 X1 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));

        for(int i=0; i<N; ++i){x[i]=i%2; y[i]=i%3;}
        result = 0.0;
        t = omp_get_wtime();
        {
            #pragma omp parallel for num_threads(4)
            for(int i=0; i<N; ++i){
                x[i] = x[i]+y[i];
            }
        }
        t = omp_get_wtime() - t;
        for(int i=0; i<N; ++i) result += x[i];
        printf("XPY STRD1 T4 X1 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));
        t4 = t; 

        for(int i=0; i<N; ++i){x[i]=i%2; y[i]=i%3;}
        result = 0.0;
        t = omp_get_wtime();
        {
            const int N4 = (N/4)*4; 
            #pragma omp parallel for num_threads(4)
            for(int i=0; i<N4; i+=4){
                x[i  ] = x[i]  +y[i];
                x[i+1] = x[i+1]+y[i+1];
                x[i+2] = x[i+2]+y[i+2];
                x[i+3] = x[i+3]+y[i+3];
            }
            for(int i=N4; i<N; ++i) x[i] = x[i] + y[i];
        }
        t = omp_get_wtime() - t;
        for(int i=0; i<N; ++i) result += x[i];
        printf("XPY STRD1 T4 X4 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", result, t, N/(t*1E9));
        printf("OpenMP speedup: %g\n\n", t1/t4); 
 
        {
            double tsin1, tsin4;
            printf("Checking OpenMP speedup with compute intensive payload\n");
            for(int i=0; i<N; ++i){x[i]=i%2; y[i]=i%3;}
            result = 0.0;
            tsin1 = omp_get_wtime();
            for(int i=0; i<N; ++i){
                x[i] += sin(y[i]);
            }
            tsin1 = omp_get_wtime() - tsin1;
            for(int i=0; i<N; ++i) result += x[i];
            printf("SIN T1 %15.12E time: %5.3f s SINPS: %5.3f\n", result, t, N/(tsin1*1E9));
            
            for(int i=0; i<N; ++i){x[i]=i%2; y[i]=i%3;}
            result = 0.0;
            tsin4 = omp_get_wtime();
            #pragma omp parallel for num_threads(4)
            for(int i=0; i<N; ++i){
                x[i] += sin(y[i]);
            }
            tsin4 = omp_get_wtime() - tsin4;
            for(int i=0; i<N; ++i) result += x[i];
            printf("SIN T4 %15.12E time: %5.3f s SINPS: %5.3f\n", result, t, N/(tsin4*1E9));
            printf("OpenMP speedup: %g\n\n", tsin1/tsin4); 
        }
 
        printf("Testing access with non-unit stride\n", t1/t); 
        for(int k=2; k<=64; k*=2){
            for(int i=0; i<N; ++i){x[i]=i%2; y[i]=i%3;}
            result = 0.0;
            t = omp_get_wtime();
            for(int j=0; j<k; ++j){
                for(int i=j; i<N; i+=k) x[i] = x[i]+y[i];
            }
            t = omp_get_wtime() - t;
            for(int i=0; i<N; ++i) result += x[i];
            printf("XPY STRD %2d T1 res: %15.12E time: %5.3f s GFLOPS: %5.3f\n", k, result, t, N/(t*1E9));
        } 
        printf("slowdown: %g times\n", t/t1);

        delete[] x;
        delete[] y;
    }
#endif
}
