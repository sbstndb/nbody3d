/*
- this example shows differnts methods to measure time or tics
- time is measured on all processes
*/

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>





typedef unsigned long long uint64 ; 


// these rdtsc are from @dssgabriel

static inline uint64 fenced_rdtscp()
{
    uint64 tsc;
    asm volatile(
        "rdtscp                  \n\t"
        "lfence                  \n\t"
        "shl     $0x20, %%rdx    \n\t"
        "or      %%rdx, %%rax    \n\t"
        : "=a" (tsc)
        :
        : "rdx", "rcx");
    return tsc;
}

uint64 rdtsc(void)
{
  uint64 a, d;
 
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  
  return (d << 32) | a;
}



void measure_time(uint64 nn , uint64 *time_rdtsc, uint64 *time_fenced_rdtscp ){

		double inc = 0.0 ;		
		uint64 start ; 
		uint64 end ; 
		uint64 boost = 1000000 ;
		
		// first burn loop + show out of order effects
		for (uint64 i = 0 ; i < boost ; i++){
			inc = inc + 1.; 
		}
		inc = 0.0 ;
		
		
		// rdtsc basic measure
		start = rdtsc() ; 
		for (uint64 i = 0 ; i < nn ; i++){
			inc = inc + 1.; 
		}
		end = rdtsc() ; 
		(*time_rdtsc) = end - start ; 		
		printf("     rdtsc timer : %llu %f\n", (*time_rdtsc), inc);
		// new burn loop
		inc = 0.0 ;
		for (uint64 i = 0 ; i < boost ; i++){
			inc = inc + 1.; 
		}
		inc = 0.0 ;	
		
			
		// fenced_rdtsc basic measure
		start = fenced_rdtscp() ; 
		for (uint64 i = 0 ; i < nn ; i++){
			inc = inc + 1.; 
		}
		end = fenced_rdtscp() ; 
		(*time_fenced_rdtscp) = end - start ; 		
		printf("   fenced  rdtscp timer : %llu %f\n", (*time_fenced_rdtscp), inc);
					
		// loop to test out of order effects
		for (uint64 i = 0 ; i < boost ; i++){
			inc = inc + 1.; 
		}		



}

void treatment(uint64 *rdtsc_mean, uint64 *rdtsc_min, uint64 *rdtsc_max, uint64 *tmp, int world_size){
	(*rdtsc_mean) += (*tmp)/(world_size - 1) ; 
	if ((*tmp) < (*rdtsc_min)){
		(*rdtsc_min) = (*tmp) ; 
	}
	if ((*tmp) > (*rdtsc_max)){
		(*rdtsc_max) = (*tmp) ; 
	}

}


int main(int argc, char** argv){

	int world_rank;
	int world_size;
	int world_res ; 
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	double inc = 0.0 ;
	uint64 max_iter = 400; 
	uint64 nn ;  	

	//if (argc==2){
	//	char *tmpstring ; 
	//	nn = strtoul(argv[1], &tmpstring, 10);
	//}	
	
	uint64 start ; 
	uint64 end ; 
	
	
	uint64 time_rdtsc ; 
	uint64 time_fenced_rdtscp ;	

	
	if (world_rank !=0 && world_rank < world_size)
	{	
		for (uint64 i = 1 ; i < max_iter ; i++){
			nn = i ; 
			measure_time(nn, &time_rdtsc, &time_fenced_rdtscp)  ;
			MPI_Send(&time_rdtsc, 1, MPI_UNSIGNED_LONG_LONG , 0 ,  4*i , MPI_COMM_WORLD) ;  // for rdtsc
			MPI_Send(&time_fenced_rdtscp, 1, MPI_UNSIGNED_LONG_LONG , 0 , 4*i + 3 , MPI_COMM_WORLD) ; // for fenced_rdtscp
		}		
			
	}
	
	if (world_rank == 0){
		FILE *fp;
		fp = fopen("time.txt" , "w") ;
		
		
		for (uint64 i = 1 ; i < max_iter ; i++){
			nn = i ; 
			printf("\nvaleur de i : %llu\n\n", i);
		  
			unsigned long long tmp = 0 ; 
		
			unsigned long long rdtsc_mean = 0 ; 	
			unsigned long long fenced_rdtscp_mean = 0 ;  


			unsigned long long rdtsc_min = 1000000000000000 ; 
			unsigned long long fenced_rdtscp_min = 1000000000000000 ; 	
		
			unsigned long long rdtsc_max = 0 ;  
			unsigned long long fenced_rdtscp_max = 0 ; 


			for (int j = 1 ; j < world_size ; j++){
			
				MPI_Recv(&tmp , 1 , MPI_UNSIGNED_LONG_LONG , j , 4 * i ,  MPI_COMM_WORLD, MPI_STATUS_IGNORE) ; // for rdtsc
				printf("     rdtsc timer : %llu with value %f\n", tmp, 0.0);
				treatment(&rdtsc_mean, &rdtsc_min, &rdtsc_max, &tmp, world_size);			
				
				MPI_Recv(&tmp , 1 , MPI_UNSIGNED_LONG_LONG , j , 4 * i + 3 ,  MPI_COMM_WORLD, MPI_STATUS_IGNORE) ; // for fenced_rdtscp
				printf("fenced rdtscp timer : %llu with value %f\n", tmp, 0.0);
				treatment(&fenced_rdtscp_mean, &fenced_rdtscp_min, &fenced_rdtscp_max, &tmp, world_size);			
				
				
								
			}
			
			printf("     rdtsc timer : %llu\n", rdtsc_mean);
			printf("fenced rdtscp timer : %llu\n", fenced_rdtscp_mean);	
			

			fprintf(fp, "%llu %llu %llu %llu %llu %llu %llu\n" , nn, rdtsc_mean , rdtsc_min , rdtsc_max , fenced_rdtscp_mean , fenced_rdtscp_min , fenced_rdtscp_max );
		}

	fclose(fp);
			
	}



	MPI_Finalize();

}
