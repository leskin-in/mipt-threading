#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <errno.h>




int n;  // Array length
int *data;




void print_help(void) {
	printf("\x1B[4mUsage\x1B[0m: mergesort n m P\n");
	printf("Qsort\n");
	
	printf("\n\n");
	
	printf("\t\x1B[1mn\x1B[0m\n");
	printf("\tNumber of elements in data.txt\n");
	printf("\n");
	
	printf("\n");
}




int qsort_compare(const void *a, const void *b) {
	int element_a = *((int*)a);
	int element_b = *((int*)b);
	
	if (element_a < element_b) {
		return -1;
	}
	else if (element_a > element_b) {
		return 1;
	}
	else {
		return 0;
	}
}




int main(int argc, char** argv) {
	
	// Launch program and read input //
	
	{
		// Check argc
		if (argc != 2) {
			print_help();
			return -1;
		}
		
		// Read arguments
		n = (int) strtol(argv[1], NULL, 10);
		if (n < 1) {
			fprintf(stderr, "Incorrect parameters.\n");
			print_help();
			return -1;
		}
		
		// Check conversions
		if (errno == ERANGE) {
			print_help();
			return -1;
		}
		
		// Allocate memory
		data = malloc(sizeof(int) * n);
		if (data == NULL) {
			fprintf(stderr, "Malloc error\n");
			return -1;
		}
		
		// Read data
		FILE *i_stream = fopen("data.txt", "r");
		if (i_stream == NULL) {
			free(data);
			fprintf(stderr, "Cannot open data file\n");
			return -1;
		}
		for (int i = 0; i < n; i++) {
			if (fscanf(i_stream, "%d", &(data[i])) == EOF) {
				free(data);
				fprintf(stderr, "File is corrupted or in wrong format\n");
				return -1;
			}
		}
		fclose(i_stream);
	}
	
	struct timeval t_start;
	gettimeofday(&t_start, NULL);
	
	qsort(data, (size_t)n, sizeof(int), qsort_compare);
	
	struct timeval t_end;
	gettimeofday(&t_end, NULL);
	
	// Output //
	
	{
		double elapsed_time = ((double) t_end.tv_sec +
							   ((double) t_end.tv_usec / 1000000.0) -
							   (double) t_start.tv_sec -
							   ((double) t_start.tv_usec / 1000000.0));
		
		FILE *o_stream;
		
		o_stream = fopen("data.txt", "a");
		if (o_stream == NULL) {
			fprintf(stderr, "Cannot open data file\n");
		}
		else {
			for (int i = 0; i < n; i++) {
				fprintf(o_stream, "%d ", data[i]);
			}
			fprintf(o_stream, "\n");
			fclose(o_stream);
		}
		
		o_stream = fopen("stats.txt", "a");
		if (o_stream == NULL) {
			fprintf(stderr, "Cannot open stats file\n");
		}
		else {
			fprintf(o_stream, "%.5fs %s\n",
					elapsed_time, argv[1]);
			fclose(o_stream);
		}
	}
	
	free(data);
	
	return 0;
}



