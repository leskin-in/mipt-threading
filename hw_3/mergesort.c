#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <errno.h>




int n;  // Array length
int m;  // Max chunk size
int P;  // Number of threads

int *data;

int *divisors;  // Chunk borders
int divisors_n;

int chunk_step;  // Number of continuous sorted chunks

int *slaves_parameters;  // pthread ID shared storage

pthread_cond_t merge_wait_point;
pthread_cond_t merge_master_wait_point;
pthread_mutex_t merge_mutex;
int merge_finished;  // Number of threads which have already finished merge iteration




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




void* chunk_qsort_thread(void* params) {
	int thread_id = *((int*)params);
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
			  sizeof(int), qsort_compare);
	}
	
	return NULL;
}




void* chunk_merge_thread(void* params) {
	int thread_id = *((int*)params);
	
	int *data_sorted = malloc(sizeof(int) * n);
	if (data_sorted == NULL) {
		fprintf(stderr, "In-thread malloc error\n");
		exit(-1);
	}
	
	int chunk_step_l;
	
	int chunk_L_id;
	int chunk_R_end_id;
	
	int dest_start;
	int dest_end;
	int dest_curr;
	int src_L;
	int src_end_L;
	int src_R;
	int src_end_R;
	
	int req_processing_last_chunk;
	
	while (1) {
		chunk_step_l = chunk_step;
		if (chunk_step_l >= divisors_n - 1) {
			break;
		}
		
		req_processing_last_chunk = 0;
		
		int j = 0;
		while (1) {
			if (req_processing_last_chunk == 0) {
				chunk_L_id =
						chunk_step_l * 2 * thread_id +
						chunk_step_l * 2 * P * j;
			}
			else {
				// Use from previous iteration
				chunk_L_id =
						chunk_step_l * 2 * thread_id +
						chunk_step_l * 2 * P * (j - 1);
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
						chunk_L_id + chunk_step_l * 2;
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
									 chunk_step_l];
			}
			else {
				// Use chunk_R_end_id for previous iteration
				src_end_L = divisors[chunk_L_id +
									 chunk_step_l * 2];
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
			
			memcpy(&data[dest_start], &data_sorted[dest_start],
				   sizeof(int) * (dest_end - dest_start));
			
			if (req_processing_last_chunk == 1) {
				break;
			}
			if ((chunk_R_end_id + chunk_step_l + 1 >=
				 divisors_n) &&
				(chunk_R_end_id + 1 < divisors_n)) {
				req_processing_last_chunk = 1;
			}
		}
		
		pthread_mutex_lock(&merge_mutex);
		merge_finished += 1;
		if (merge_finished == P) {
			pthread_cond_broadcast(&merge_master_wait_point);
			while (chunk_step_l == chunk_step) {
				pthread_cond_wait(&merge_wait_point, &merge_mutex);
			}
			pthread_mutex_unlock(&merge_mutex);
		}
		else {
			while (chunk_step_l == chunk_step) {
				pthread_cond_wait(&merge_wait_point, &merge_mutex);
			}
			pthread_mutex_unlock(&merge_mutex);
		}
	}
	
	return NULL;
}




int main(int argc, char** argv) {
	pthread_t *slaves;
	
	// Launch program and read input //
	
	{
		// Check argc
		if (argc != 4) {
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
		m = (int) strtol(argv[2], NULL, 10);
		if (m < 2) {
			fprintf(stderr, "Incorrect parameters.\n");
			print_help();
			return -1;
		}
		P = (int) strtol(argv[3], NULL, 10);
		if (P < 1) {
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
			free(data);
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
	
	// Prepare data //
	
	{
		int m_actual = ((n / P >= m) ? (m) : (n / P));
		
		// Use one extra chunk if there are some elements left after division
		if (n % m_actual == 0) {
			divisors_n = n / m_actual;
		}
		else {
			divisors_n = n / m_actual + 1;
		}
		divisors_n += 1;  // Add an extra divisor at the beginning
		
		divisors = malloc(sizeof(int) * divisors_n);
		if (divisors == NULL) {
			free(data);
			fprintf(stderr, "Malloc error\n");
			return -1;
		}
		
		divisors[0] = 0;
		for (int i = 1; i < divisors_n - 1; i++) {
			divisors[i] = divisors[i - 1] + m_actual;
		}
		divisors[divisors_n - 1] = n;
		
		// Prepare slave threads structures
		slaves = malloc(sizeof(pthread_t) * P);
		if (slaves == NULL) {
			free(data);
			free(divisors);
			fprintf(stderr, "Malloc error\n");
			return -1;
		}
		slaves_parameters = malloc(sizeof(int) * P);
		if (slaves_parameters == NULL) {
			free(data);
			free(divisors);
			free(slaves);
			fprintf(stderr, "Malloc error\n");
			return -1;
		}
		for (int i = 0; i < P; i++) {
			slaves_parameters[i] = i;
		}
	}
	
	// Quick sort of every chunk //
	
	{
		for (int i = 0; i < P; i++) {
			if (pthread_create(&slaves[i], NULL,
							   chunk_qsort_thread, &slaves_parameters[i]) !=
				0) {
				fprintf(stderr, "Thread launch error\n");
				return -1;
			}
		}
		for (int i = 0; i < P; i++) {
			if (pthread_join(slaves[i], NULL) != 0) {
				fprintf(stderr, "Thread join error\n");
				return -1;
			}
		}
	}
	
	// Chunk merge //
	
	{
		chunk_step = 1;
		merge_finished = 0;
		pthread_mutex_init(&merge_mutex, NULL);
		pthread_cond_init(&merge_wait_point, NULL);
		pthread_cond_init(&merge_master_wait_point, NULL);
		
		for (int i = 0; i < P; i++) {
			if (pthread_create(&slaves[i], NULL,
							   chunk_merge_thread, &slaves_parameters[i]) !=
				0) {
				fprintf(stderr, "Thread launch error\n");
				return -1;
			}
		}
		
		while (1) {
			pthread_mutex_lock(&merge_mutex);
			if (chunk_step >= divisors_n - 1) {
				pthread_mutex_unlock(&merge_mutex);
				break;
			}
			
			while (merge_finished < P) {
				pthread_cond_wait(&merge_master_wait_point, &merge_mutex);
			}
			
			chunk_step *= 2;
			merge_finished = 0;
			pthread_cond_broadcast(&merge_wait_point);
			pthread_mutex_unlock(&merge_mutex);
		}
		
		for (int i = 0; i < P; i++) {
			if (pthread_join(slaves[i], NULL) != 0) {
				fprintf(stderr, "Thread join error\n");
				return -1;
			}
		}
	}
	
	free(divisors);
	free(slaves);
	free(slaves_parameters);
	
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
			fprintf(o_stream, "%.5fs %s %s %s\n",
					elapsed_time, argv[1], argv[2], argv[3]);
			fclose(o_stream);
		}
	}
	
	free(data);
	
	return 0;
}



