#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <sys/time.h>
#include <errno.h>




void print_help(void) {
	printf("\x1B[4mUsage\x1B[0m: mergesort n m P\n");
	printf("Multi-threaded qsort-based mergesort\n");
	
	printf("\n\n");
	
	printf("\t\x1B[1mn\x1B[0m\n");
	printf("\tNumber of elements in data.txt\n");
	printf("\n");
	
	printf("\t\x1B[1mm\x1B[0m\n");
	printf("\tMaximum chunk size\n");
	printf("\n");
	
	printf("\t\x1B[1mP\x1B[0m\n");
	printf("\tNumber of threads\n");
	printf("\n");
	
	printf("\n");
}




int int_comparator(const void *a, const void *b) {
	if (a > b) {
		return -1;
	}
	else if (a < b) {
		return 1;
	}
	else {
		return 0;
	}
}




int main(int argc, char** argv) {
	int n;  // Array length
	
	int m;  // Max chunk size
	
	int P;  // Number of OpenMP threads
	
	int *data;
	int *data_sorted;
	
	int i;
	struct timeval t_start;
	struct timeval t_end;
	
	// Prepare data //
	{
		// Check argc
		if (argc != 4) {
			print_help();
			return -1;
		}
		
		// Read arguments
		n = (int)strtol(argv[1], NULL, 10);
		m = (int)strtol(argv[2], NULL, 10);
		P = (int)strtol(argv[3], NULL, 10);
		
		// Check conversions
		if (errno == ERANGE) {
			print_help();
			return -1;
		}
		
		// Allocate memory
		data = malloc(sizeof(int) * n);
		data_sorted = malloc(sizeof(int) * n);
		if ((data == NULL) ||
				(data_sorted == NULL)) {
			free(data);
			free(data_sorted);
			fprintf(stderr, "Malloc error\n");
			return -1;
		}
		
		// Read data
		FILE* i_stream = fopen("data.txt", "r");
		if (i_stream == NULL) {
			free(data);
			free(data_sorted);
			fprintf(stderr, "Cannot open data file\n");
			return -1;
		}
		for (i = 0; i < n; i++) {
			fscanf(i_stream, "%d", &(data[i]));
		}
		fclose(i_stream);
		
		// Prepare OpenMP
		omp_set_num_threads(P);
	}
	
	// Sorting //
	{
		int m_actual[P];
		int start_pos[P];
		
		// Define start positions and boundaries
		start_pos[0] = 0;
		m_actual[0] = n / P;
		if (m_actual[0] >= m) {
			m_actual[0] = m;
			for (i = 1; i < P; i++) {
				m_actual[i] = m_actual[0];
				start_pos[i] = m_actual[0] * i;
			}
		}
		else {
			int sum = m_actual[0];
			for (i = 1; i < P; i++) {
				m_actual[i] = m_actual[0];
				start_pos[i] = m_actual[0] * i;
				sum += m_actual[0];
			}
			if (sum < n) {
				m_actual[P - 1] += (n - sum);
			}
		}
		
		gettimeofday(&t_start, NULL);

		#pragma omp parallel
		{
			#pragma omp single
			{
				for (i = 0; i < P; i++) {
					#pragma omp task shared(data, m_actual, n, P)
					{
						int cpos = start_pos[omp_get_thread_num()];
						int epos = cpos + m_actual[omp_get_thread_num()];
						do {
							if (cpos >= n) {
								break;
							}
							else if (epos > n) {
								epos = n;
							}
							
							qsort(&(data[cpos]), (size_t)epos - cpos,
								  sizeof(int), int_comparator);
							
							cpos += m_actual[omp_get_thread_num()] * P;
							epos = cpos + m_actual[omp_get_thread_num()];
						} while (1);
					}
				}
			}
			#pragma omp taskwait
			
			#pragma omp single
			{
				int current_P = P;
				int united_P = 1;
				
				int merge_starts[P];
				int merge_ends[P];
				
				while (current_P > 0) {
					if (current_P > 1) {
						// FIXME: Sometimes P reduction should not be performed, but chunks should be merged
						// Reduce the number of merge segments
						current_P /= 2;
						united_P *= 2;
						
						// Assing new merge starts and ends
						merge_starts[0] = 0;
						for (i = 0; i < current_P; i++) {
							merge_ends[i] = merge_starts[i];
							int j;
							for (j = i * united_P; j < i * united_P + united_P; j++) {
								if (j >= P) {
									break;
								}
 								merge_ends[i] += m_actual[j];
							}
							if (i + 1 < current_P) {
								merge_starts[i + 1] = merge_ends[i];
							}
						}
						merge_ends[current_P - 1] = n;
					}
					else {
						merge_starts[0] = 0;
						merge_ends[0] = n;
					}
					
					// Do merge sort
					for (i = 0; i < current_P; i++) {
						#pragma omp task shared(data, data_sorted, merge_starts, merge_ends)
						{
							int dest_start = merge_starts[omp_get_thread_num()];
							int dest_end = merge_ends[omp_get_thread_num()];
							int src_L = dest_start;
							int src_end_L = (dest_start + dest_end) / 2;
							int src_R = src_end_L;
							int src_end_R = dest_end;
							
							int dest_curr = dest_start;
							while (dest_curr < dest_end) {
								if (src_R == src_end_R) {
									data_sorted[dest_curr] = data[src_L];
									src_L += 1;
									dest_curr += 1;
								}
								else if (src_L == src_end_L) {
									data_sorted[dest_curr] = data[src_R];
									src_R += 1;
									dest_curr += 1;
								}
								else if (data[src_L] < data[src_R]) {
									data_sorted[dest_curr] = data[src_L];
									src_L += 1;
									dest_curr += 1;
								}
								else {
									data_sorted[dest_curr] = data[src_R];
									src_R += 1;
									dest_curr += 1;
								}
							}
							dest_curr = dest_start;
							while (dest_curr < dest_end) {
								data[dest_curr] = data_sorted[dest_curr];
								dest_curr += 1;
							}
						}
					}
					#pragma omp taskwait
					
					if (current_P == 1) {
						current_P = 0;
					}
				}
			}
		}
		
		gettimeofday(&t_end, NULL);
	}
	
	// Output //
	{
		double elapsed_time = ((double)t_end.tv_sec +
				((double)t_end.tv_usec / 1000000.0) -
				(double)t_start.tv_sec -
				((double)t_start.tv_usec / 1000000.0));

		FILE* o_stream;
		
		o_stream = fopen("data.txt", "a");
		if (o_stream == NULL) {
			fprintf(stderr, "Cannot open data file\n");
		}
		else {
			for (i = 0; i < n; i++) {
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
			fprintf(o_stream, "%.5fs %s %s %s\n",
					elapsed_time, argv[1], argv[2], argv[3]);
			fclose(o_stream);
		}
		
		free(data);
		free(data_sorted);
	}
	
	return 0;
}



