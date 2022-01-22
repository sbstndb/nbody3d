#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

int n = 10240 ; 
int repetition = 100000; 

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
double *v_x, *v_y, *v_x_al, *v_y_al, *asx, *asy, *asxy, *asx2, *asy2; 
double *result ; 



int main(){
	printf("-- This is the Pearson benchmark --\n");
	v_aos = malloc(n * sizeof(vector));
	v_x = malloc(n * sizeof(double)) ; 
	v_y = malloc(n * sizeof(double)) ; 	
	v_x_al = aligned_alloc(n, n * sizeof(double)) ; 
	v_y_al = aligned_alloc(n, n * sizeof(double)) ; 
	asx = aligned_alloc(64, 64*sizeof(double));	
	asy = aligned_alloc(64, 64*sizeof(double));	
	asxy = aligned_alloc(64, 64*sizeof(double));	
	asx2 = aligned_alloc(64, 64*sizeof(double));	
	asy2 = aligned_alloc(64, 64*sizeof(double));	
	printf("allocated\n");

	double before, after ;

	
	// initialisation 
  	
	for (int i = 0 ; i < n ; i++){
		double value_x = randreal(); 
		double value_y = randreal() ;	
		v_aos[i].x = value_x ; 
		v_aos[i].y = value_y ; 	
		v_x[i] = value_x ; 
		v_y[i] = value_y ; 
		v_x_al[i] = value_x ; 
		v_y_al[i] = value_y ; 			
	}
	for (int i = 0 ; i < 64 ; i++){
		asx[i] = 0.0 ; 
		asy[i] = 0.0 ; 
		asxy[i] = 0.0 ; 
		asx2[i] = 0.0 ; 
		asy2[i] = 0.0 ; 
	}
	printf("initialized\n");	
	//bad AOS version
	
	double rho = 0.0 ; 
	double sum_x = 0.0 ; 
	double sum_y = 0.0 ; 
	double sum_x_2 = 0.0 ; 
	double sum_y_2 = 0.0 ;
	double sum_xy = 0.0 ; 
	double num = 0.0 ;
	double den = 0.0 ; 
	
	before = (double)rdtsc();	
	for (int j = 0 ; j < repetition ; j++){
		rho = 0.0 ; 
		sum_x = 0.0 ; 
		sum_y = 0.0 ; 
		sum_x_2 = 0.0 ; 
		sum_y_2 = 0.0 ;
		sum_xy = 0.0 ; 
		num = 0.0 ;
		den = 0.0 ; 	
		for (int i = 0 ; i < n ; i++){
			sum_xy += v_aos[i].x * v_aos[i].y ;
			sum_x += v_aos[i].x ;
			sum_x_2 += v_aos[i].x * v_aos[i].x ;		
			sum_y += v_aos[i].y ;
			sum_y_2 += v_aos[i].y * v_aos[i].y ;
			//printf(": %f\n", v_aos[i].x) ; 			
		}
		num = n * sum_xy - sum_x * sum_y ; 
		den = sqrt(n * sum_x_2 - sum_x * sum_x) * sqrt(n * sum_y_2 - sum_y * sum_y) ;
		rho = num / den  ;		
	}
	after = (double)rdtsc();
	
	printf(" \n") ; 			
	printf("AOS time : %lf\n", ((after - before)/repetition));
	printf("AOS num : %lf\n", num) ; 
	printf("AOS den : %lf\n", den) ; 
	printf("AOS rho : %lf\n", rho) ; 	
	
	
	
	//SOA version
	before = (double)rdtsc();	
	for (int j = 0 ; j < repetition ; j++){	
		rho = 0.0 ; 
		sum_x = 0.0 ; 
		sum_y = 0.0 ; 
		sum_x_2 = 0.0 ; 
		sum_y_2 = 0.0 ;
		sum_xy = 0.0 ; 
		num = 0.0 ;
		den = 0.0 ;
		for (int i = 0 ; i < n ; i++){
			sum_xy += v_x[i] * v_y[i] ;
			sum_x += v_x[i] ;
			sum_x_2 += v_x[i] * v_x[i] ;		
			sum_y += v_y[i] ;
			sum_y_2 += v_y[i] * v_y[i] ;
			//printf(": %f\n", v_aos[i].x) ; 			
		}
		num = n * sum_xy - sum_x * sum_y ; 
		den = sqrt(n * sum_x_2 - sum_x * sum_x) * sqrt(n * sum_y_2 - sum_y * sum_y) ;
		rho = num / den  ;		
	}
	after = (double)rdtsc();
	
	printf(" \n") ; 		
	printf("SOA time : %lf\n", ((after - before)/repetition));
	printf("SOA num : %lf\n", num) ; 
	printf("SOA den : %lf\n", den) ; 
	printf("SOA rho : %lf\n", rho) ; 	
	
	//Aligned SOA version	
	before = (double)rdtsc();	
	for (int j = 0 ; j < repetition ; j++){	
		rho = 0.0 ; 
		sum_x = 0.0 ; 
		sum_y = 0.0 ; 
		sum_x_2 = 0.0 ; 
		sum_y_2 = 0.0 ;
		sum_xy = 0.0 ; 
		num = 0.0 ;
		den = 0.0 ;
		for (int i = 0 ; i < n ; i++){
			sum_xy += v_x_al[i] * v_y_al[i] ;
			sum_x += v_x_al[i] ;
			sum_x_2 += v_x_al[i] * v_x_al[i] ;		
			sum_y += v_y_al[i] ;
			sum_y_2 += v_y_al[i] * v_y_al[i] ;
			//printf(": %f\n", v_aos[i].x) ; 			
		}
		num = n * sum_xy - sum_x * sum_y ; 
		den = sqrt(n * sum_x_2 - sum_x * sum_x) * sqrt(n * sum_y_2 - sum_y * sum_y) ;
		rho = num / den  ;		
	}
	after = (double)rdtsc();
	
	printf(" \n") ; 		
	printf("Aligned SOA time : %lf\n", ((after - before)/repetition));
	printf("Aligned SOA num : %lf\n", num) ; 
	printf("Aligned SOA den : %lf\n", den) ; 
	printf("Aligned SOA rho : %lf\n", rho) ; 	



	// intrinsic version 	
	__m256d rx, ry, rxy, rx2, ry2, sx, sy, sxy, sx2, sy2; 
 
	before = (double)rdtsc();	
	for (int j = 0 ; j < repetition ; j++){
		rho = 0.0 ; 
		sum_x = 0.0 ; 
		sum_y = 0.0 ; 
		sum_x_2 = 0.0 ; 
		sum_y_2 = 0.0 ;
		sum_xy = 0.0 ; 
		num = 0.0 ;
		den = 0.0 ;
		
		sx = _mm256_setzero_pd() ;
		sy = _mm256_setzero_pd() ;
		sxy = _mm256_setzero_pd() ;
		sx2 = _mm256_setzero_pd() ;	
		sy2 = _mm256_setzero_pd() ;	
		for (int i = 0 ; i < n ; i+=4){
			rx = _mm256_load_pd(&v_x_al[i]);			
			ry = _mm256_load_pd(&v_y_al[i]);			
			rx2 = _mm256_mul_pd(rx, rx) ; 			
			ry2 = _mm256_mul_pd(ry, ry) ;			
			rxy = _mm256_mul_pd(rx, ry) ;	
			sx = _mm256_add_pd(sx, rx);
			sy = _mm256_add_pd(sy, ry);
			sx2 = _mm256_add_pd(sx2, rx2);
			sy2 = _mm256_add_pd(sy2, ry2);
			sxy = _mm256_add_pd(sxy, rxy);						
		}
		_mm256_store_pd(&asx[0], sx) ;
		_mm256_store_pd(&asy[0], sy) ;	
		_mm256_store_pd(&asxy[0], sxy) ;	
		_mm256_store_pd(&asx2[0], sx2) ;	
		_mm256_store_pd(&asy2[0], sy2) ;
		for (int i = 0 ; i < 4 ; i++){
			sum_x += asx[i] ; 
			sum_y += asy[i] ;
			sum_xy += asxy[i] ; 
			sum_x_2 += asx2[i] ;
			sum_y_2 += asy2[i] ;
			
		}	
		num = n * sum_xy - sum_x * sum_y ; 
		den = sqrt(n * sum_x_2 - sum_x * sum_x) * sqrt(n * sum_y_2 - sum_y * sum_y) ;
		rho = num / den  ;							

		rho = num / den  ;		
	}
	after = (double)rdtsc();
	printf("\n") ; 			
	printf("intrinsic SOA time : %lf\n", ((after - before)/repetition));
	printf("intrinsic SOA num : %lf\n", num) ; 
	printf("intrinsic SOA den : %lf\n", den) ; 
	printf("intrinsic SOA rho : %lf\n", rho) ; 	

	// unrolled intrinsic version 	
	__m256d rxu, ryu, rxyu, rx2u, ry2u, sxu, syu, sxyu, sx2u, sy2u; 
 
	before = (double)rdtsc();	
	for (int j = 0 ; j < repetition ; j++){
		rho = 0.0 ; 
		sum_x = 0.0 ; 
		sum_y = 0.0 ; 
		sum_x_2 = 0.0 ; 
		sum_y_2 = 0.0 ;
		sum_xy = 0.0 ; 
		num = 0.0 ;
		den = 0.0 ;
		
		sx = _mm256_setzero_pd() ;
		sy = _mm256_setzero_pd() ;
		sxy = _mm256_setzero_pd() ;
		sx2 = _mm256_setzero_pd() ;	
		sy2 = _mm256_setzero_pd() ;	
		sxu = _mm256_setzero_pd() ;
		syu = _mm256_setzero_pd() ;
		sxyu = _mm256_setzero_pd() ;
		sx2u = _mm256_setzero_pd() ;	
		sy2u = _mm256_setzero_pd() ;
		for (int i = 0 ; i < n ; i+=8){
			rx = _mm256_load_pd(&v_x_al[i]);			
			ry = _mm256_load_pd(&v_y_al[i]);			
			rx2 = _mm256_mul_pd(rx, rx) ; 			
			ry2 = _mm256_mul_pd(ry, ry) ;			
			rxy = _mm256_mul_pd(rx, ry) ;	
			sx = _mm256_add_pd(sx, rx);
			sy = _mm256_add_pd(sy, ry);
			sx2 = _mm256_add_pd(sx2, rx2);
			sy2 = _mm256_add_pd(sy2, ry2);
			sxy = _mm256_add_pd(sxy, rxy);	
			
			rxu = _mm256_load_pd(&v_x_al[i+4]);			
			ryu = _mm256_load_pd(&v_y_al[i+4]);			
			rx2u = _mm256_mul_pd(rxu, rxu) ; 			
			ry2u = _mm256_mul_pd(ryu, ryu) ;			
			rxyu = _mm256_mul_pd(rxu, ryu) ;	
			sxu = _mm256_add_pd(sxu, rxu);
			syu = _mm256_add_pd(syu, ryu);
			sx2u = _mm256_add_pd(sx2u, rx2u);
			sy2u = _mm256_add_pd(sy2u, ry2u);
			sxyu = _mm256_add_pd(sxyu, rxyu);				
					
		}
		_mm256_store_pd(&asx[0], sx) ;
		_mm256_store_pd(&asy[0], sy) ;	
		_mm256_store_pd(&asxy[0], sxy) ;	
		_mm256_store_pd(&asx2[0], sx2) ;	
		_mm256_store_pd(&asy2[0], sy2) ;
		
		_mm256_store_pd(&asx[4], sxu) ;
		_mm256_store_pd(&asy[4], syu) ;	
		_mm256_store_pd(&asxy[4], sxyu) ;	
		_mm256_store_pd(&asx2[4], sx2u) ;	
		_mm256_store_pd(&asy2[4], sy2u) ;
		for (int i = 0 ; i < 8 ; i++){
			sum_x += asx[i] ; 
			sum_y += asy[i] ;
			sum_xy += asxy[i] ; 
			sum_x_2 += asx2[i] ;
			sum_y_2 += asy2[i] ;
			
		}	
		num = n * sum_xy - sum_x * sum_y ; 
		den = sqrt(n * sum_x_2 - sum_x * sum_x) * sqrt(n * sum_y_2 - sum_y * sum_y) ;
		rho = num / den  ;							

		rho = num / den  ;		
	}
	after = (double)rdtsc();
	printf("\n") ; 			
	printf("unrolled intrinsic SOA time : %lf\n", ((after - before)/repetition));
	printf("unrolled intrinsic SOA num : %lf\n", num) ; 
	printf("unrolled intrinsic SOA den : %lf\n", den) ; 
	printf("unrolled intrinsic SOA rho : %lf\n", rho) ;



	// double unrolled intrinsic version 	
	__m256d rxu3, ryu3, rxyu3, rx2u3, ry2u3, sxu3, syu3, sxyu3, sx2u3, sy2u3, rxu4, ryu4, rxyu4, rx2u4, ry2u4, sxu4, syu4, sxyu4, sx2u4, sy2u4; 
 
	before = (double)rdtsc();	
	for (int j = 0 ; j < repetition ; j++){
		rho = 0.0 ; 
		sum_x = 0.0 ; 
		sum_y = 0.0 ; 
		sum_x_2 = 0.0 ; 
		sum_y_2 = 0.0 ;
		sum_xy = 0.0 ; 
		num = 0.0 ;
		den = 0.0 ;
		
		sx = _mm256_setzero_pd() ;
		sy = _mm256_setzero_pd() ;
		sxy = _mm256_setzero_pd() ;
		sx2 = _mm256_setzero_pd() ;	
		sy2 = _mm256_setzero_pd() ;	
		sxu = _mm256_setzero_pd() ;
		syu = _mm256_setzero_pd() ;
		sxyu = _mm256_setzero_pd() ;
		sx2u = _mm256_setzero_pd() ;	
		sy2u = _mm256_setzero_pd() ;
		sxu3 = _mm256_setzero_pd() ;
		syu3 = _mm256_setzero_pd() ;
		sxyu3 = _mm256_setzero_pd() ;
		sx2u3 = _mm256_setzero_pd() ;	
		sy2u3 = _mm256_setzero_pd() ;
		sxu4 = _mm256_setzero_pd() ;
		syu4 = _mm256_setzero_pd() ;
		sxyu4 = _mm256_setzero_pd() ;
		sx2u4 = _mm256_setzero_pd() ;	
		sy2u4 = _mm256_setzero_pd() ;
		for (int i = 0 ; i < n ; i+=16){
			rx = _mm256_load_pd(&v_x_al[i]);			
			ry = _mm256_load_pd(&v_y_al[i]);			
			rx2 = _mm256_mul_pd(rx, rx) ; 			
			ry2 = _mm256_mul_pd(ry, ry) ;			
			rxy = _mm256_mul_pd(rx, ry) ;	
			sx = _mm256_add_pd(sx, rx);
			sy = _mm256_add_pd(sy, ry);
			sx2 = _mm256_add_pd(sx2, rx2);
			sy2 = _mm256_add_pd(sy2, ry2);
			sxy = _mm256_add_pd(sxy, rxy);	
			
			rxu = _mm256_load_pd(&v_x_al[i+4]);			
			ryu = _mm256_load_pd(&v_y_al[i+4]);			
			rx2u = _mm256_mul_pd(rxu, rxu) ; 			
			ry2u = _mm256_mul_pd(ryu, ryu) ;			
			rxyu = _mm256_mul_pd(rxu, ryu) ;	
			sxu = _mm256_add_pd(sxu, rxu);
			syu = _mm256_add_pd(syu, ryu);
			sx2u = _mm256_add_pd(sx2u, rx2u);
			sy2u = _mm256_add_pd(sy2u, ry2u);
			sxyu = _mm256_add_pd(sxyu, rxyu);				

			rxu3 = _mm256_load_pd(&v_x_al[i+8]);			
			ryu3 = _mm256_load_pd(&v_y_al[i+8]);			
			rx2u3 = _mm256_mul_pd(rxu3, rxu3) ; 			
			ry2u3 = _mm256_mul_pd(ryu3, ryu3) ;			
			rxyu3 = _mm256_mul_pd(rxu3, ryu3) ;	
			sxu3 = _mm256_add_pd(sxu3, rxu3);
			syu3 = _mm256_add_pd(syu3, ryu3);
			sx2u3 = _mm256_add_pd(sx2u3, rx2u3);
			sy2u3 = _mm256_add_pd(sy2u3, ry2u3);
			sxyu3 = _mm256_add_pd(sxyu3, rxyu3);

			rxu4 = _mm256_load_pd(&v_x_al[i+12]);			
			ryu4 = _mm256_load_pd(&v_y_al[i+12]);			
			rx2u4 = _mm256_mul_pd(rxu4, rxu4) ; 			
			ry2u4 = _mm256_mul_pd(ryu4, ryu4) ;			
			rxyu4 = _mm256_mul_pd(rxu4, ryu4) ;	
			sxu4 = _mm256_add_pd(sxu4, rxu4);
			syu4 = _mm256_add_pd(syu4, ryu4);
			sx2u4 = _mm256_add_pd(sx2u4, rx2u4);
			sy2u4 = _mm256_add_pd(sy2u4, ry2u4);
			sxyu4 = _mm256_add_pd(sxyu4, rxyu4);
					
		}
		_mm256_store_pd(&asx[0], sx) ;
		_mm256_store_pd(&asy[0], sy) ;	
		_mm256_store_pd(&asxy[0], sxy) ;	
		_mm256_store_pd(&asx2[0], sx2) ;	
		_mm256_store_pd(&asy2[0], sy2) ;
		
		_mm256_store_pd(&asx[4], sxu) ;
		_mm256_store_pd(&asy[4], syu) ;	
		_mm256_store_pd(&asxy[4], sxyu) ;	
		_mm256_store_pd(&asx2[4], sx2u) ;	
		_mm256_store_pd(&asy2[4], sy2u) ;
		
		_mm256_store_pd(&asx[8], sxu3) ;
		_mm256_store_pd(&asy[8], syu3) ;	
		_mm256_store_pd(&asxy[8], sxyu3) ;	
		_mm256_store_pd(&asx2[8], sx2u3) ;	
		_mm256_store_pd(&asy2[8], sy2u3) ;
		
		_mm256_store_pd(&asx[12], sxu4) ;
		_mm256_store_pd(&asy[12], syu4) ;	
		_mm256_store_pd(&asxy[12], sxyu4) ;	
		_mm256_store_pd(&asx2[12], sx2u4) ;	
		_mm256_store_pd(&asy2[12], sy2u4) ;
		for (int i = 0 ; i < 16 ; i++){
			sum_x += asx[i] ; 
			sum_y += asy[i] ;
			sum_xy += asxy[i] ; 
			sum_x_2 += asx2[i] ;
			sum_y_2 += asy2[i] ;
			
		}	
		num = n * sum_xy - sum_x * sum_y ; 
		den = sqrt(n * sum_x_2 - sum_x * sum_x) * sqrt(n * sum_y_2 - sum_y * sum_y) ;
		rho = num / den  ;							

		rho = num / den  ;		
	}
	after = (double)rdtsc();
	printf("\n") ; 			
	printf("double unrolled intrinsic SOA time : %lf\n", ((after - before)/repetition));
	printf("double unrolled intrinsic SOA num : %lf\n", num) ; 
	printf("double unrolled intrinsic SOA den : %lf\n", den) ; 
	printf("double unrolled intrinsic SOA rho : %lf\n", rho) ;


	
	
}

