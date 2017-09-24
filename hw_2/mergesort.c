#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
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
	int n;  // Array length
	
	int m;  // Max chunk size
	
	int P;  // Number of OpenMP threads
	
	int *data;
	int *data_sorted;
	
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
		for (int i = 0; i < n; i++) {
			if (fscanf(i_stream, "%d", &(data[i])) == EOF) {
				free(data);
				free(data_sorted);
				fprintf(stderr, "Incorrect file content\n");
				return -1;
			}
		}
		fclose(i_stream);
		
		// Prepare OpenMP
		omp_set_num_threads(P);
	}
	
	
	// Sort //
	
	{
		int m_actual = ((n / P >= m) ? (m) : (n / P));
		
		int divisors_n;
		// Use 1 extra chunk if there are some elements left after division
		if (n % m_actual == 0) {
			divisors_n = n / m_actual;
		}
		else {
			divisors_n = n / m_actual + 1;
		}
		divisors_n += 1;  // Add an extra divisor at the beginning
		
		int *divisors = malloc(sizeof(int) * divisors_n);
		if (divisors == NULL) {
			free(data);
			free(data_sorted);
			fprintf(stderr, "Malloc error\n");
			return -1;
		}
		
		divisors[0] = 0;
		for (int i = 1; i < divisors_n - 1; i++) {
			divisors[i] = divisors[i-1] + m_actual;
		}
		divisors[divisors_n - 1] = n;
		
		int chunk_step = 1;
		
		
		gettimeofday(&t_start, NULL);

		#pragma omp parallel default(none) shared(divisors_n, P, data, divisors)
		{
			#pragma omp single nowait
			{
				#pragma omp taskgroup
				{
					for (int i = 0; i < P; i++) {
						#pragma omp task
						{
							int thread_id = omp_get_thread_num();
							int chunk_id;
							
							int j = 0;
							while (1) {
								chunk_id = thread_id + j * P;
								j += 1;
								if (chunk_id + 1 >= divisors_n) {
									break;
								}
								qsort(&(data[divisors[chunk_id]]),
									  (size_t) divisors[chunk_id + 1] -
									  divisors[chunk_id],
									  sizeof(int), int_comparator);
							}
						}
					}
				}
			}
		}
		
		
		memcpy(data_sorted, data, sizeof(int) * n);
		
		while (chunk_step < divisors_n - 1) {
			#pragma omp parallel default(none) shared(divisors_n, chunk_step, P, data, data_sorted, divisors)
			{
				#pragma omp single nowait
				{
					#pragma omp taskgroup
					{
						for (int i = 0; i < P; i++) {
							#pragma omp task
							{
								int thread_id = omp_get_thread_num();
								
								int chunk_L_id;
								int chunk_R_end_id;
								
								int dest_start;
								int dest_end;
								int dest_curr;
								int src_L;
								int src_end_L;
								int src_R;
								int src_end_R;
								
								int req_processing_last_chunk = 0;
								
								int j = 0;
								while (1) {
									if (req_processing_last_chunk == 0) {
										chunk_L_id =
												chunk_step * 2 * thread_id +
												chunk_step * 2 * P * j;
									}
									else {
										// Use from previous iteration
										chunk_L_id =
												chunk_step * 2 * thread_id +
												chunk_step * 2 * P * (j - 1);
									}
									j += 1;
									
									// If this chunk is the last one or invalid
									if ((chunk_L_id + 2 >= divisors_n) ||
										(chunk_L_id < 0)) {
										break;
									}
									
									dest_start = divisors[chunk_L_id];
									
									if (req_processing_last_chunk == 0) {
										chunk_R_end_id =
												chunk_L_id + chunk_step * 2;
										if (chunk_R_end_id >= divisors_n) {
											chunk_R_end_id = divisors_n - 1;
										}
									}
									else {
										chunk_R_end_id = divisors_n - 1;
									}
									dest_end = divisors[chunk_R_end_id];
									
									dest_curr = dest_start;
									
									src_L = dest_start;
									
									if (req_processing_last_chunk == 0) {
										src_end_L = divisors[chunk_L_id +
															 chunk_step];
									}
									else {
										// Use chunk_R_end_id for previous iteration
										src_end_L = divisors[chunk_L_id +
															 chunk_step * 2];
									}
									src_R = src_end_L;
									src_end_R = dest_end;
									
									while (dest_curr < dest_end) {
										if (src_R == src_end_R) {
											data_sorted[dest_curr] = data[src_L];
											src_L += 1;
										}
										else if (src_L == src_end_L) {
											data_sorted[dest_curr] = data[src_R];
											src_R += 1;
										}
										else if (data[src_L] < data[src_R]) {
											data_sorted[dest_curr] = data[src_L];
											src_L += 1;
										}
										else {
											data_sorted[dest_curr] = data[src_R];
											src_R += 1;
										}
										dest_curr += 1;
									}
									
									if (req_processing_last_chunk == 1) {
										break;
									}
									if ((chunk_R_end_id + chunk_step + 1 >=
										 divisors_n) &&
										(chunk_R_end_id + 1 < divisors_n)) {
										memcpy(&data[dest_start],
											   &data_sorted[dest_start],
											   sizeof(int) *
											   (dest_end - dest_start));
										req_processing_last_chunk = 1;
									}
								}
							}
						}
					}
				}
			}
			
			chunk_step *= 2;
			memcpy(data, data_sorted, sizeof(int) * n);
		}
		
		gettimeofday(&t_end, NULL);
		
		free(divisors);
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
			fprintf(o_stream, "%.5fs %s %s %s\n",
					elapsed_time, argv[1], argv[2], argv[3]);
			fclose(o_stream);
		}
		
		free(data);
		free(data_sorted);
	}
	
	return 0;
}



