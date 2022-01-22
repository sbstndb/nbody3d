#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

int n = 1024 ; 
int repetition = 1000; 

typedef struct{
	float x , y  ; 
} vector ;


unsigned long long rdtsc(void)
{
  unsigned long long a, d;
  
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  
  return (d << 32) | a;
}



int randxy(int x, int y)
{
  return (rand() % (y - x + 1)) + x; 
}

//
double randreal()
{
  int s = (randxy(0, 1)) ? 1 : -1;
  int a = randxy(1, RAND_MAX), b = randxy(1, RAND_MAX);

  return s * ((double)a / (double)b); 
}



vector *v_aos; 
float *v_x, *v_y, *v_x_al, *v_y_al ; 
float *result ; 



int main(){
	printf("This is the add - mul benchmark\n");
	v_aos = malloc(n * sizeof(vector));
	v_x = malloc(n * sizeof(float)) ; 
	v_y = malloc(n * sizeof(float)) ; 	
	v_x_al = aligned_alloc(n, n * sizeof(float)) ; 
	v_y_al = aligned_alloc(n, n * sizeof(float)) ; 
	result = aligned_alloc(64, 64*sizeof(float));	



	
	// initialisation 
  	
	for (int i = 0 ; i < n ; i++){
		float value_x = randreal(); 
		float value_y = randreal();	
		v_aos[i].x = value_x ; 
		v_aos[i].y = value_y ; 	
		v_x[i] = value_x ; 
		v_y[i] = value_y ; 
		v_x_al[i] = value_x ; 
		v_y_al[i] = value_y ; 			
	}
	for (int i = 0 ; i < 64 ; i++){
		result[i] = 0.0 ; 
	}
	//AOS version
	double sum = 0.0; 
	double before = (double)rdtsc();
	for (int j = 0 ; j < repetition ; j++){
	
		for (int i = 0 ; i < n ; i++){
			sum += v_aos[i].x * v_aos[i].y ; 
		}
	}
	double after = (double)rdtsc();
	printf("AOS time : %lf\n", ((after - before)/repetition));
	printf("AOS sum : %f\n", sum) ; 
	
	// SOA version
	sum = 0.0 ; 		
	before = (double)rdtsc();
	for (int j = 0 ; j < repetition ; j++){

		for (int i = 0 ; i < n ; i++){
			sum += v_x[i] * v_y[i] ; 
		}
	}
	after = (double)rdtsc();
	printf("SOA time : %lf\n", ((after - before)/repetition));
	printf("SOA sum : %f\n", sum) ; 
	
	// bad intrinsic version 
	__m256 rx, ry,rp, rs, rx2, ry2, rp2, rs2, rx3, ry3, rp3, rs3, rx4, ry4, rp4, rs4 ;
	sum = 0 ;                
	before = (double)rdtsc();
	for (int j = 0 ; j < repetition ; j++){
		for ( int i = 0 ; i < n ; i+=8){
			rx = _mm256_load_ps(&v_x_al[i]);
			ry = _mm256_load_ps(&v_y_al[i] );
			rp = _mm256_mul_ps(rx, ry) ; 
			_mm256_store_ps(result, rp) ; 
			for (int k = 0 ; k < 8 ; k++){
				//printf("%f\n", result[i]);
				sum += result[k] ; 		
			}	
		} 
	}
	after = (double)rdtsc();
	printf("bad intrinsic time : %lf\n", ((after - before)/repetition));
	printf("bad intrinsic sum : %f\n", sum) ; 	
	
	
	// fma intrinsic version 
	sum = 0 ;                             
	before = (double)rdtsc();
	for (int j = 0 ; j < repetition ; j++){	
	rs = _mm256_setzero_ps() ;	
		for ( int i = 0 ; i < n ; i+=8){
			rx = _mm256_load_ps(&v_x_al[i]);
			ry = _mm256_load_ps(&v_y_al[i] );
			rp = _mm256_mul_ps(rx, ry) ; 
			rs = _mm256_add_ps(rp, rs);
		} 
		_mm256_store_ps(result, rs) ; 		
		for (int k = 0 ; k < 8 ; k++){
			//printf("%f\n", result[k]);
			sum += result[k] ; 		
		}		
		
	}
	after = (double)rdtsc();
	printf("fma intrinsic time : %lf\n", ((after - before)/repetition));
	printf("fma intrinsic sum : %f\n", sum) ; 	


	// unrolled fma intrinsic version 
	sum = 0 ;                             
	before = (double)rdtsc();
	for (int j = 0 ; j < repetition ; j++){	
	rs = _mm256_setzero_ps() ;	
	rs2 = _mm256_setzero_ps() ;
		for ( int i = 0 ; i < n ; i+=16){
			rx = _mm256_load_ps(&v_x_al[i]);
			ry = _mm256_load_ps(&v_y_al[i] );
			rp = _mm256_mul_ps(rx, ry) ; 
			rs = _mm256_add_ps(rp, rs);
			
			rx2 = _mm256_load_ps(&v_x_al[i+8]);
			ry2 = _mm256_load_ps(&v_y_al[i+8] );
			rp2 = _mm256_mul_ps(rx2, ry2) ; 
			rs2 = _mm256_add_ps(rp2, rs2);
		} 
		_mm256_store_ps(&result[0], rs) ;
		_mm256_store_ps(&result[8], rs2) ;		 		
		for (int k = 0 ; k < 16 ; k++){
			//printf("%f\n", result[k]);
			sum += result[k] ; 		
		}	
	}
	after = (double)rdtsc();
	printf("unrolled fma intrinsic time : %lf\n", ((after - before)/repetition));
	printf("unrolled fma intrinsic sum : %f\n", sum) ; 
	
	
	// double unrolled fma intrinsic version 
	sum = 0 ;                             
	before = (double)rdtsc();
	for (int j = 0 ; j < repetition ; j++){	
	rs = _mm256_setzero_ps() ;	
	rs2 = _mm256_setzero_ps() ;
	rs3 = _mm256_setzero_ps() ;	
	rs4 = _mm256_setzero_ps() ;	
	
		for ( int i = 0 ; i < n ; i+=32){
			rx = _mm256_load_ps(&v_x_al[i]);
			ry = _mm256_load_ps(&v_y_al[i] );
			rp = _mm256_mul_ps(rx, ry) ; 
			rs = _mm256_add_ps(rp, rs);
			
			rx2 = _mm256_load_ps(&v_x_al[i+8]);
			ry2 = _mm256_load_ps(&v_y_al[i+8] );
			rp2 = _mm256_mul_ps(rx2, ry2) ; 
			rs2 = _mm256_add_ps(rp2, rs2);
			
			rx3 = _mm256_load_ps(&v_x_al[i+16]);
			ry3 = _mm256_load_ps(&v_y_al[i+16] );
			rp3 = _mm256_mul_ps(rx3, ry3) ; 
			rs3 = _mm256_add_ps(rp3, rs3);
			
			rx4 = _mm256_load_ps(&v_x_al[i+24]);
			ry4 = _mm256_load_ps(&v_y_al[i+24] );
			rp4 = _mm256_mul_ps(rx4, ry4) ; 
			rs4 = _mm256_add_ps(rp4, rs4);			
			
		} 
		_mm256_store_ps(&result[0], rs) ; 	
		_mm256_store_ps(&result[8], rs2) ; 
		_mm256_store_ps(&result[16], rs3) ; 
		_mm256_store_ps(&result[24], rs4) ; 	
		for (int k = 0 ; k < 64 ; k++){
			//printf("%f\n", result[k]);
			sum += result[k] ;	
		}	

		
	}
	after = (double)rdtsc();
	printf("double unrolled fma intrinsic time : %lf\n", ((after - before)/repetition));
	printf("double unrolled fma intrinsic sum : %f\n", sum) ; 	


}

