#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <errno.h>




void print_help(void) {
	printf("\x1B[4mUsage\x1B[0m: mccalc a b x N p P\n");
	printf("Multi-threaded Monte-Carlo probability calculator\n");
	
	printf("\n\n");
	
	printf("\t\x1B[1ma\x1B[0m\n");
	printf("\tLeft border\n");
	printf("\n");
	
	printf("\t\x1B[1mb\x1B[0m\n");
	printf("\tRight border\n");
	printf("\n");
	
	printf("\t\x1B[1mx\x1B[0m\n");
	printf("\tStart position\n");
	printf("\n");
	
	printf("\t\x1B[1mN\x1B[0m\n");
	printf("\tNumber of tests\n");
	printf("\n");
	
	printf("\t\x1B[1mp\x1B[0m\n");
	printf("\tTrue probability\n");
	printf("\n");
	
	printf("\t\x1B[1mP\x1B[0m\n");
	printf("\tNumber of threads to run\n");
	printf("\n");
	
	printf("\n");
}




int main(int argc, char** argv) {
	const double DOUBLE_ERROR = 1e-8;
	
	long a, b;
	long x_start;
	long N;
	double p_right;
	int threads_total;
	
	unsigned long *wandering_time;
	char *wandering_result;
	
	int i;
	
	struct timeval t_start;
	struct timeval t_end;
	
	
	// Prepare data //
	{
		// Check argc
		if (argc != 7) {
			print_help();
			return -1;
		}
		
		// Read arguments
		a = strtol(argv[1], NULL, 10);
		b = strtol(argv[2], NULL, 10);
		x_start = strtol(argv[3], NULL, 10);
		N = strtol(argv[4], NULL, 10);
		p_right = strtod(argv[5], NULL);
		threads_total = (int)strtol(argv[6], NULL, 10);
		
		// Check conversions
		if (errno == ERANGE) {
			print_help();
			return -1;
		}
		
		// Allocate memory
		wandering_time = malloc(sizeof(unsigned long) * N);
		wandering_result = malloc((size_t)N);
		if ((wandering_time == NULL) ||
			(wandering_result == NULL)) {
			free(wandering_result);
			free(wandering_time);
			fprintf(stderr, "Malloc error\n");
			return -1;
		}
	}
	
	// Calculations //
	{
		long current_position;
		unsigned long current_wandering_time;
		
		struct drand48_data rand_seed;
		double current_rand;
		
		
		gettimeofday(&t_start, NULL);
		
		for (i = 0; i < N; i++) {
			current_position = x_start;
			current_wandering_time = 0;
			
			srand48_r(t_start.tv_usec / 1000 * i + t_start.tv_sec, &rand_seed);
			
			while (1) {
				// Check break conditions
				if (current_position <= a) {
					wandering_result[i] = -1;
					wandering_time[i] = current_wandering_time;
					break;
				}
				if (current_position >= b) {
					wandering_result[i] = 1;
					wandering_time[i] = current_wandering_time;
					break;
				}
				
				// Change position
				drand48_r(&rand_seed, &current_rand);
				if (current_rand - p_right < DOUBLE_ERROR) {
					current_position += 1;
				}
				else {
					current_position -= 1;
				}
				
				current_wandering_time += 1;
			}
		}
		
		gettimeofday(&t_end, NULL);
	}
	
	// Output //
	{
		long double average_wandering_time = 0.0;
		double p_experimental = 0.0;
		double elapsed_time = ((double)t_end.tv_sec +
							   ((double)t_end.tv_usec / 1000000.0) -
							   (double)t_start.tv_sec -
							   ((double)t_start.tv_usec / 1000000.0));
		
		for (i = 0; i < N; i++) {
			average_wandering_time += (long double)wandering_time[i];
			if (wandering_result[i] > 0) {
				p_experimental += 1;
			}
		}
		average_wandering_time /= (long double)N;
		p_experimental /= (double)N;
		
		FILE* o_stream = fopen("stats.txt", "a");
		fprintf(o_stream, "%.4f %.1Lf %.3fs %s %s %s %s %s %s\n",
				p_experimental, average_wandering_time, elapsed_time,
				argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
		fclose(o_stream);
		
		free(wandering_result);
		free(wandering_time);
	}
	
	return 0;
}



