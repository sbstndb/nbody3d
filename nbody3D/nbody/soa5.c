//
// structure : soa
// optimization : avx512 intrinsic functions
//
//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <unistd.h>

#define SAVE

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
//typedef struct particle_s {
//
//  f32 x, y, z;
//  f32 vx, vy, vz;
//  
//} particle_t;

//
void init(float *x, float *y, float *z, float *vx, float *vy, float *vz, float *softening, float softening_value, float *dt, float dt_value,  u64 n)
{
	for (u64 i = 0 ; i < 16 ; i++){
		x[i] = 0.0 ; 
		y[i] = 0.0 ; 
		z[i] = 0.0 ;
		vx[i] = 0.0 ; 
		vy[i] = 0.0 ; 
		vz[i] = 0.0 ; 
	}
	for (u64 i = 8; i < n + 16; i++)
	{
		//
		u64 r1 = (u64)rand();
		u64 r2 = (u64)rand();
		f32 sign = (r1 > r2) ? 1 : -1;

		//
		x[i] = sign * (f32)rand() / (f32)RAND_MAX;
		y[i] = (f32)rand() / (f32)RAND_MAX;
		z[i] = sign * (f32)rand() / (f32)RAND_MAX;

		//
		vx[i] = 0.0;
		vy[i] = 0.0;
		vz[i] = 0.0;
	}
	for (u64 i  = 0  ; i < 128 ; i++){
		softening[i] = softening_value;
		dt[i] = dt_value ; 
	}
	for (u64 i = n + 16 ; i < n + 32 ; i++){
		x[i] = 0.0 ; 
		y[i] = 0.0 ; 
		z[i] = 0.0 ;
		vx[i] = 0.0 ; 
		vy[i] = 0.0 ; 
		vz[i] = 0.0 ; 
	}	
	for (int i = 0 ; i < n+32 ; i++){
	//printf("%f\n",vx[i]);
	//sleep(1);
	}
	
}

//
void move_particles(float *x, float *y, float *z, float *vx, float *vy, float *vz, float *softening, float *dt, u64 n){

	float test ; 
	
	//printf("valeur test : %f\n", x[n+8]) ; 
	//
	__m512 rxi, ryi, rzi, rxj, ryj, rzj , rfx, rfy, rfz;
	__m512 rdx  , rdy , rdz , rdxyz , rsoft, rt;
	
	rxi = _mm512_setzero_ps() ;
	ryi = _mm512_setzero_ps() ;	
	rzi = _mm512_setzero_ps() ;	

	
	rxj = _mm512_setzero_ps() ;
	ryj = _mm512_setzero_ps() ;	
	rzj = _mm512_setzero_ps() ;	

	
	rdx = _mm512_setzero_ps() ;
	rdy = _mm512_setzero_ps() ;	
	rzj = _mm512_setzero_ps() ;	
	rdz = _mm512_setzero_ps() ;	
	rdxyz = _mm512_setzero_ps() ;	
	rsoft = _mm512_load_ps(&softening[0]) ;
	rt = _mm512_load_ps(&dt[0]) ;	
	for (u64 i = 1; i < n + 16; i++)
	{


		rfx = _mm512_setzero_ps() ;	
		rfy = _mm512_setzero_ps() ;	
		rfz = _mm512_setzero_ps() ; 
		// load i values
		rxi = _mm512_loadu_ps(&x[i]);
		ryi = _mm512_loadu_ps(&y[i]);
		rzi = _mm512_loadu_ps(&z[i]);			
		for (u64 j = 16; j < n + 1; j+=16)
		{

						
			rxj = _mm512_load_ps(&x[j]);	
			ryj = _mm512_load_ps(&y[j]);
			rzj = _mm512_load_ps(&z[j]);
			


			rdx = _mm512_sub_ps(rxj, rxi) ; 
			rdy = _mm512_sub_ps(ryj, ryi) ; 
			rdz = _mm512_sub_ps(rzj, rzi) ;


			 
			rxj = _mm512_mul_ps(rdx, rdx) ;
			rzj = _mm512_mul_ps(rdz, rdz) ;
			

						
			
			rdxyz = _mm512_add_ps(rxj, rzj) ;
			ryj =   _mm512_mul_ps(rdy, rdy) ;
			rdxyz = _mm512_add_ps(rdxyz, ryj) ;				
			rdxyz = _mm512_add_ps(rsoft, rdxyz) ;	

			rxj = _mm512_rsqrt14_ps(rdxyz);
		
			
			rdxyz = _mm512_mul_ps(rxj, rxj) ;	
			rdxyz = _mm512_mul_ps(rdxyz, rxj) ;
						
			rdx = _mm512_mul_ps(rdx, rdxyz) ; 			
			rdy = _mm512_mul_ps(rdy, rdxyz) ; 
			rdz = _mm512_mul_ps(rdz, rdxyz) ; 
	
			rfx = _mm512_add_ps(rfx, rdx) ; 	
			rfy = _mm512_add_ps(rfy, rdy) ; 	
			rfz = _mm512_add_ps(rfz, rdz) ; 
			
							
		}
		
			//sleep(1) ;	
		rfx = _mm512_mul_ps(rfx, rt) ; 	
		rfy = _mm512_mul_ps(rfy, rt) ; 	
		rfz = _mm512_mul_ps(rfz, rt) ; 
 				
		
		rxi = _mm512_loadu_ps(&vx[i]);
		ryi = _mm512_loadu_ps(&vy[i]);
		rzi = _mm512_loadu_ps(&vz[i]);
		
	
		
		rxi = _mm512_add_ps(rxi, rfx) ; 	
		ryi = _mm512_add_ps(ryi, rfy) ; 
		rzi = _mm512_add_ps(rzi, rfz) ; 
		
		
		
		_mm512_storeu_ps(&vx[i], rxi) ;				
		_mm512_storeu_ps(&vy[i], ryi) ;	
		_mm512_storeu_ps(&vz[i], rzi) ;
		
		
		rfx = _mm512_setzero_ps() ;	
		rfy = _mm512_setzero_ps() ;	
		rfz = _mm512_setzero_ps() ;
		


		
	}
	
	
	
	//3 floating-point operations

	for (u64 i = 0 ; i < 16 ; i++){
		x[i] = 0.0 ; 
		y[i] = 0.0 ; 
		z[i] = 0.0 ;
		vx[i] = 0.0 ; 
		vy[i] = 0.0 ; 
		vz[i] = 0.0 ; 
	}
	for (u64 i = n + 16 ; i < n + 32 ; i++){
		x[i] = 0.0 ; 
		y[i] = 0.0 ; 
		z[i] = 0.0 ;
		vx[i] = 0.0 ; 
		vy[i] = 0.0 ; 
		vz[i] = 0.0 ; 
	}	
			
	for (u64 i = 16; i < n + 16; i++)
	{
	x[i] += dt[0] * vx[i];
	y[i] += dt[0] * vy[i];
	z[i] += dt[0] * vz[i];
	}
}

//
int main(int argc, char **argv)
{
	//
  const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
  const u64 warmup = (argc > 2) ? atoll(argv[2]) : 3;   
  const u64 steps = (argc > 3) ? atoll(argv[3]) : 10;  
	const f32 dt_value = 0.01;
	const f32 softening_value = 1e-20;
	
	
	// declaration of file for saving 1st coordinates during time
	#ifdef SAVE
		FILE *xfilePtr = NULL ; 
		xfilePtr = fopen("nbodyx5.txt", "w");
		FILE *vfilePtr = NULL ; 
		vfilePtr = fopen("nbodyv5.txt", "w");		
		if (xfilePtr == NULL || vfilePtr == NULL){
		printf("Issue in writing in file\n") ; 
		}
		char buf[100] ;   
	#endif	
	
	//
	f64 rate = 0.0, drate = 0.0;


	//

	float *x = aligned_alloc(64, sizeof(float) * (n+128)) ; 
	float *y = aligned_alloc(64, sizeof(float) * (n+128)) ; 	
	float *z = aligned_alloc(64, sizeof(float) * (n+128)) ; 	
	float *vx = aligned_alloc(64, sizeof(float) * (n+128)) ; 	
	float *vy = aligned_alloc(64, sizeof(float) * (n+128)) ; 	
	float *vz = aligned_alloc(64, sizeof(float) * (n+128)) ; 	
	float *softening = aligned_alloc(64, sizeof(float) * 128) ; 
	float *dt = aligned_alloc(64, sizeof(float) * 128) ; 	
	//
	init(x, y,z, vx, vy, vz, softening, softening_value, dt, dt_value,  n);

	const u64 sx = 6 *  sizeof(double) * n;

	printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", sx, sx >> 10, sx >> 20);

	//
	printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);

	//
	for (u64 i = 0; i < steps; i++)
	{
	

		
	
	//Measure
	const f64 start = omp_get_wtime();

	move_particles(x, y, z, vx, vy, vz, softening, dt, n);

	const f64 end = omp_get_wtime();

	//Number of interactions/iterations
	const f32 h1 = (f32)(n) * (f32)(n - 1);

	//GFLOPS
	const f32 h2 = (23.0 * h1 + 3.0 * (f32)n) * 1e-9;

	if (i >= warmup)
	{
		rate += h2 / (end - start);
		drate += (h2 * h2) / ((end - start) * (end - start));
	}

	//
	printf("%5llu %10.3e %10.3e %8.1f %s\n",
	i,
	(end - start),
	h1 / (end - start),
	h2 / (end - start),
	(i < warmup) ? "*" : "");

	fflush(stdout);
	}
	
	
      #ifdef SAVE
	      for (unsigned long long b = 0 ; b < n ; b++){
		fputs(gcvt(x[b+16], 16, buf), xfilePtr)  ; 
		fputs(" ", xfilePtr)  ; 
		fputs(gcvt(y[b+16], 16, buf), xfilePtr)  ; 
		fputs(" ", xfilePtr)  ;       
		fputs(gcvt(z[b+16], 16, buf), xfilePtr)  ; 
		fputs(" \n", xfilePtr)  ;     
		
		fputs(gcvt(vx[b+16], 16, buf), vfilePtr)  ; 
		fputs(" ", vfilePtr)  ; 
		fputs(gcvt(vy[b+16], 16, buf), vfilePtr)  ; 
		fputs(" ", vfilePtr)  ;       
		fputs(gcvt(vz[b+16], 16, buf), vfilePtr)  ; 
		fputs(" \n", vfilePtr)  ;  	
		}	         
      #endif 		
	

	//
	rate /= (f64)(steps - warmup);
	drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

	printf("-----------------------------------------------------\n");
	printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
	"Average performance:", "", rate, drate);
	printf("-----------------------------------------------------\n");

	//
	free(x);
	free(y);
	free(z);
	free(vx);
	free(vy);
	free(vz);
	free(softening);
	free(dt);
	
	#ifdef SAVE
		fclose(xfilePtr); 
		fclose(vfilePtr) ; 
	#endif	

	//
	return 0;
}
