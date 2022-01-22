//
// structure : soa
// optimization of sqrt etc, vectorization and unrolling loops
//
//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

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
void init(float *x, float *y, float *z, float *vx, float *vy, float *vz, u64 n)
{
	for (u64 i = 0; i < n; i++)
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
}

//
void move_particles(float *x, float *y, float *z, float *vx, float *vy, float *vz, const f32 dt, u64 n)
{
	//
	const f32 softening = 1e-20;
	float tmp = 0.0 ; 
	float tmp2 = 0.0 ; 
	float tmp3 = 0.0 ; 
	float tmp4 = 0.0 ; 
	

	//
	

	for (u64 i = 0; i < n; i++)
	{
	//
		f32 fx = 0.0;
		f32 fy = 0.0;
		f32 fz = 0.0;

		f32 fx2 = 0.0;
		f32 fy2 = 0.0;
		f32 fz2 = 0.0;

		//23 floating-point operations
		for (u64 j = 0; j < n; j+=2)
		{
		  //Newton's law
		  const f32 dx = x[j] - x[i]; //1
		  const f32 dx2 = x[j+1] - x[i]; //1
		  const f32 dy = y[j] - y[i]; //2
		  const f32 dy2 = y[j+1] - y[i]; //2
		  const f32 dz = z[j] - z[i]; //3
		  const f32 dz2 = z[j+1] - z[i]; //3
		  const f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; //9
		  const f32 d_22 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2) + softening; //9
		  tmp = sqrt(d_2) ; 
		  tmp3 = sqrt(d_22) ; 
		  tmp2 = tmp * tmp * tmp ; 
		  tmp4 = tmp3 * tmp3 * tmp3 ; 
		  const f32 d_3_over_2 = 1./tmp2 ; 
		  const f32 d_3_over_22 = 1./tmp4 ; 		  

		  //Net force
		  fx += dx / d_3_over_2; //13
		  fy += dy / d_3_over_2; //15
		  fz += dz / d_3_over_2; //17

		  //Net force
		  fx2 += dx2 / d_3_over_22; //13
		  fy2 += dy2 / d_3_over_22; //15
		  fz2 += dz2 / d_3_over_22; //17		  
		  
		}

		//
		vx[i] += dt * fx; //19
		vy[i] += dt * fy; //21
		vz[i] += dt * fz; //23
	
		vx[i] += dt * fx2; //19
		vy[i] += dt * fy2; //21
		vz[i] += dt * fz2; //23
	}

	//3 floating-point operations
	for (u64 i = 0; i < n; i++)
	{
	x[i] += dt * vx[i];
	y[i] += dt * vy[i];
	z[i] += dt * vz[i];
	}
}

//
int main(int argc, char **argv)
{
	//
	const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
	const u64 steps= 10;
	const f32 dt = 0.01;

	// declaration of file for saving 1st coordinates during time
	#ifdef SAVE
	  FILE *xfilePtr = NULL ; 
	  xfilePtr = fopen("soa3x.txt", "w");
	  FILE *vfilePtr = NULL ; 
	  vfilePtr = fopen("soa3v.txt", "w");	  
	  if (xfilePtr == NULL || vfilePtr == NULL){
	  	printf("Issue in writing in file\n") ; 
	  }
	  char buf[100] ;   
	#endif

	//
	f64 rate = 0.0, drate = 0.0;

	//Steps to skip for warm up
	const u64 warmup = 3;

	//
	//particle_t *p = malloc(sizeof(particle_t) * n);
	float *x = aligned_alloc(64, sizeof(float) * n) ; 
	float *y = aligned_alloc(64, sizeof(float) * n) ; 	
	float *z = aligned_alloc(64, sizeof(float) * n) ; 	
	float *vx = aligned_alloc(64, sizeof(float) * n) ; 	
	float *vy = aligned_alloc(64, sizeof(float) * n) ; 	
	float *vz = aligned_alloc(64, sizeof(float) * n) ; 	
	
	//
	init(x, y,z, vx, vy, vz, n);

	const u64 sx = 6 *  sizeof(double) * n;

	printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", sx, sx >> 10, sx >> 20);

	//
	printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);

	//
	for (u64 i = 0; i < steps; i++)
	{
	
	#ifdef SAVE
		// write 1st trajectory particle for comparison 
		fputs(gcvt(x[0], 16, buf), xfilePtr)  ; 
		fputs(" ", xfilePtr)  ; 
		fputs(gcvt(y[0], 16, buf), xfilePtr)  ; 
		fputs(" ", xfilePtr)  ;       
		fputs(gcvt(z[0], 16, buf), xfilePtr)  ; 
		fputs(" \n", xfilePtr)  ;     
		
		fputs(gcvt(vx[0], 16, buf), vfilePtr)  ; 
		fputs(" ", vfilePtr)  ; 
		fputs(gcvt(vy[0], 16, buf), vfilePtr)  ; 
		fputs(" ", vfilePtr)  ;       
		fputs(gcvt(vz[0], 16, buf), vfilePtr)  ; 
		fputs(" \n", vfilePtr)  ;  		 
	#endif	
	
	//Measure
	const f64 start = omp_get_wtime();

	move_particles(x, y, z, vx, vy, vz, dt, n);

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
	
	#ifdef SAVE
		fclose(xfilePtr);
		fclose(vfilePtr) ;  
	#endif	

	//
	return 0;
}
