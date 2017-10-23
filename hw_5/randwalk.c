#include <errno.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <sys/time.h>





const double DOUBLE_ERROR = 1e-7;




long l;
long a, b;
long n;
long N;
double p_left;
double p_right;
double p_up;
double p_down;




void print_help(void) {
	printf("\x1B[4mUsage\x1B[0m: randwalk l a b n N p_left p_right p_up p_down\n");
	printf("Monte-Carlo random walk modelling MPI-based program\n");
	
	printf("\n\n");
	
	printf("\t\x1B[1ml\x1B[0m\n");
	printf("\tSquare side length\n");
	printf("\n");
	
	printf("\t\x1B[1ma\x1B[0m\n");
	printf("\tNumber of squares at X axis (horizontal)\n");
	printf("\n");
	
	printf("\t\x1B[1mb\x1B[0m\n");
	printf("\tNumber of squares at Y axis (vertical)\n");
	printf("\n");
	
	printf("\t\x1B[1mn\x1B[0m\n");
	printf("\tNumber of steps (for each particle)\n");
	printf("\n");
	
	printf("\t\x1B[1mN\x1B[0m\n");
	printf("\tNumber of particles\n");
	printf("\n");
	
	printf("\t\x1B[1mp_left\x1B[0m\n");
	printf("\tProbability of going to the left\n");
	printf("\n");
	
	printf("\t\x1B[1mp_right\x1B[0m\n");
	printf("\tProbability of going to the right\n");
	printf("\n");
	
	printf("\t\x1B[1mp_up\x1B[0m\n");
	printf("\tProbability of going up\n");
	printf("\n");
	
	printf("\t\x1B[1mp_down\x1B[0m\n");
	printf("\tProbability of going down\n");
	printf("\n");
	
	printf("\n");
}




int main(int argc, char** argv) {
	// Read arguments //
	{
		// Check argc
		if (argc != 10) {
			print_help();
			return -1;
		}
		
		// Convert arguments
		l = strtol(argv[1], NULL, 10);
		a = strtol(argv[2], NULL, 10);
		b = strtol(argv[3], NULL, 10);
		n = strtol(argv[4], NULL, 10);
		N = strtol(argv[5], NULL, 10);
		p_left = strtod(argv[6], NULL);
		p_right = strtod(argv[7], NULL);
		p_up = strtod(argv[8], NULL);
		p_down = strtod(argv[9], NULL);
		
		// Check conversions
		if (errno != 0) {
			print_help();
			return -1;
		}
	}
	
	
	
	
	
	return 0;
}




