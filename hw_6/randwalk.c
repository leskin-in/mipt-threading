#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <mpi.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>

#include "randgen.h"




// Custom types //

/// Particle type
typedef struct particle_t {
	long x;
	long y;
	long origin_sector;
	long steps_to_live;
} PARTICLE;

/// Message buffer (for internal communications) for particles
struct msgbuf_particle {
	long type;
	PARTICLE particle;
};




// Global service variables //

const size_t SIZEOF_MSGBUF_PARTICLE = 5 * sizeof(long);

const double DOUBLE_ERROR = 1e-7;

/// Interval (in ticks) at which the receive operation is performed
const int RECEIVE_TICK_INTERVAL = 64;

const int QUEUE_RAND_NORMAL_LENGTH = 512;


int mpi_size;

int mpi_rank;


/// Global MPI_Datatype for the PARTICLE struct
MPI_Datatype PARTICLE_mpi_t;


/// System V queue ID for the sender thread
int qid_send;

pthread_t thread_sender_tid;


/// System V queue ID for the receiver thread
int qid_receive;

pthread_t thread_receiver_tid;


/// System V message type and IPC tag: normal message
#define CLASSIFIER_NORM 1

/// System V message type and IPC tag: cease operations
#define CLASSIFIER_STOP 2

/// IPC tag: particle death report
#define CLASSIFIER_REPORT 3

/// IPC tag: message containing elapsed time of a node
#define CLASSIFIER_TIME 4


/// MPI ID of master unit (it performs a few extra operations)
#define MPI_MASTER_RANK 0




// Global variables for particle processing //

/// Sector side
long l;
/// Sector number at X axis
long a;
/// Sector number at Y axis
long b;
/// Total number of steps a particle should do
long n;
/// Total number of particles
long N;
/// Total number of dead particles
long N_dead;

double p_left;
double p_right;
double p_up;
double p_down;

/// Sector offset length
long offset;

/// Maximum value of 'offset'
const long OFFSET_MAX = 64;

/// Divisor of 'offset' calculation. 'offset' = 'l' / 'OFFSET_DIVISOR'
const long OFFSET_DIVISOR = 5;

long field_width;
long field_height;

/// Depot for the active particles (currently processed)
PARTICLE* depot_a;
/// Total number of the particles currently being processed
int depot_a_length;

/// Depot for the dead particles (which died in current thread's action zone)
PARTICLE* depot_d;
/// Total number of the dead particles in current thread's action zone
int depot_d_length;




// Service & processing functions //

void print_help(void) {
	printf("\x1B[4mUsage\x1B[0m: randwalk l a b n N p_left p_right p_up p_down\n");
	printf("Monte-Carlo random walk modelling MPI-based program\n");
	printf("IMPORTANT: Total number of squares must be equal to a number of MPI processes\n");
	
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



/**
 * @return Sector the 'unit' actually belongs to
 */
long get_sector(const PARTICLE* unit) {
	return (unit ->y / l) * a + (unit ->x / l);
}



/**
 * @param current_sector ID of the sector which holds 'unit' at the moment
 * @return Sector that should process 'unit'
 */
long get_sector_with_offset(const PARTICLE* unit, long current_sector) {
	long ruled_sector = get_sector(unit);
	if (ruled_sector == current_sector) {
		return current_sector;
	}
	
	long x_min = l * (current_sector % a);
	long y_min = l * (current_sector / a);
	long x_max = x_min + l - 1;
	long y_max = y_min + l - 1;
	x_min -= offset;
	y_min -= offset;
	x_max += offset;
	y_max += offset;
	
	if (x_min < 0) {
		x_min += field_width;
	}
	else if (x_max > field_width) {
		x_max -= field_width;
	}
	if ((unit ->x < x_min) && (unit ->x > x_max)) {
		return ruled_sector;
	}
	
	if (y_min < 0) {
		y_min += field_height;
	}
	else if (y_max > field_height) {
		y_max -= field_height;
	}
	if ((unit ->y < y_min) && (unit ->y > y_max)) {
		return ruled_sector;
	}
	
	return current_sector;
}




/**
 * @brief Add a patricle to 'depot_a'
 * @return 0 if successful
 * @return -1 if there is not enough memory available for the 'realloc()'
 */
int depot_add(const PARTICLE* unit) {
	depot_a_length += 1;
	depot_a = realloc(depot_a, depot_a_length * sizeof(PARTICLE));
	if (errno != 0) {
		return -1;
	}
	depot_a[depot_a_length - 1] = *unit;
	return 0;
}



/**
 * @brief Drop a particle from 'depot_a'
 * @noexcept
 */
void depot_remove(int depot_a_id) {
	if ((depot_a_id < 0) || (depot_a_id >= depot_a_length)) {
		return;
	}
	depot_a[depot_a_id] = depot_a[depot_a_length - 1];
	depot_a_length -= 1;
	depot_a = realloc(depot_a, depot_a_length * sizeof(PARTICLE));
	errno = 0;
}



/**
 * @brief Kill a particle from 'depot_a' and remove it from there
 * @param depot_a_id Sequence number of a particle in 'depot_a'
 * @return 0 if successful
 * @return -1 if the 'realloc()' failed or a given ID is incorrect
 */
int depot_die(int depot_a_id) {
	if ((depot_a_id < 0) || (depot_a_id >= depot_a_length)) {
		return -1;
	}
	depot_d_length += 1;
	depot_d = realloc(depot_d, depot_d_length * sizeof(PARTICLE));
	if (errno != 0) {
		return -1;
	}
	depot_d[depot_d_length - 1] = depot_a[depot_a_id];
	depot_remove(depot_a_id);
	return 0;
}



/**
 * @brief Thread function for the sender thread
 */
void* thread_sender(void* context) {
	// Check launch
	{
		int flag;
		MPI_Initialized(&flag);
		if ((context != NULL) || (flag == 0)) {
			fprintf(stderr, "Sender thread is not launched properly\n");
			return NULL;
		}
	}
	
	/// Message buffer for the System V message queue
	struct msgbuf_particle buffer;
	MPI_Request dummy;
	
	while (1) {
		if (msgrcv(qid_send,
				   &buffer, SIZEOF_MSGBUF_PARTICLE, 0, MSG_NOERROR) < 0) {
			// Something bad happened, but do nothing
			errno = 0;
			continue;
		}
		
		if (buffer.type == CLASSIFIER_STOP) {
			// Cease operations
			break;
		}
		
		if (buffer.type == CLASSIFIER_REPORT) {
			// Send report to everyone (including self)
			for (int i = 0; i < mpi_size; i++) {
				MPI_Isend(&buffer.particle,
						  1, PARTICLE_mpi_t, i, CLASSIFIER_REPORT,
						  MPI_COMM_WORLD, &dummy);
				MPI_Request_free(&dummy);
			}
			continue;
		}
		
		MPI_Isend(&buffer.particle,
				  1, PARTICLE_mpi_t,
				  (int)get_sector(&buffer.particle), CLASSIFIER_NORM,
				  MPI_COMM_WORLD, &dummy);
		MPI_Request_free(&dummy);
	}
	
	return NULL;
}



/**
 * @brief Thread function for the receiver thread
 */
void* thread_receiver(void* context) {
	// Check launch
	{
		int flag;
		MPI_Initialized(&flag);
		if ((context != NULL) || (flag == 0)) {
			fprintf(stderr, "Receiver thread is not launched properly\n");
			return NULL;
		}
	}
	
	/// Message buffer for the System V message queue
	struct msgbuf_particle buffer;
	MPI_Status status;
	
	while (1) {
		MPI_Recv(&buffer.particle,
				 1, PARTICLE_mpi_t, MPI_ANY_SOURCE, MPI_ANY_TAG,
				 MPI_COMM_WORLD, &status);
		
		if (status.MPI_TAG == CLASSIFIER_REPORT) {
			// Some particle has died
			N_dead += 1;
			if (N_dead == N) {
				// Cease operations: Send message both to main() and sender
				buffer.type = CLASSIFIER_STOP;
				msgsnd(qid_send, &buffer, SIZEOF_MSGBUF_PARTICLE, 0);
				msgsnd(qid_receive, &buffer, SIZEOF_MSGBUF_PARTICLE, 0);
				break;
			}
		}
		
		// Send the particle received for processing
		buffer.type = CLASSIFIER_NORM;
		msgsnd(qid_receive, &buffer, SIZEOF_MSGBUF_PARTICLE, 0);
	}
	
	return NULL;
}




int main(int argc, char** argv) {
	// Initialize MPI //
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	
	
	// Read arguments and initialize key variables //
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
		N_dead = 0;
		p_left = strtod(argv[6], NULL);
		p_right = strtod(argv[7], NULL);
		p_up = strtod(argv[8], NULL);
		p_down = strtod(argv[9], NULL);
		
		// Check conversions
		if (errno != 0) {
			print_help();
			return -1;
		}
		
		if (a * b != mpi_size) {
			fprintf(stderr, "A number of MPI units is incorrect\n");
			print_help();
			return -1;
		}
		
		offset = ((l / OFFSET_DIVISOR) > OFFSET_MAX) ?
				 (OFFSET_MAX) :
				 (l / OFFSET_DIVISOR);
		
		field_width = a * l;
		field_height = b * l;
	}
	
	struct timeval t_start;
	gettimeofday(&t_start, NULL);
	
	
	// Initialize PARTICE_mpi_t //
	{
		PARTICLE example;
		
		MPI_Datatype sequence[3] = {MPI_LONG, MPI_LONG, MPI_LONG};
		
		int sequence_amounts[3] = {2, 1, 1};
		
		MPI_Aint sequence_displacements[3];
		MPI_Get_address(&example, &sequence_displacements[0]);
		MPI_Get_address(&example.origin_sector, &sequence_displacements[1]);
		MPI_Get_address(&example.steps_to_live, &sequence_displacements[2]);
		for (int i = 1; i < 3; i++) {
			sequence_displacements[i] -= sequence_displacements[0];
		}
		sequence_displacements[0] = 0;
		
		MPI_Type_create_struct(1,
							   sequence_amounts, sequence_displacements, sequence, &PARTICLE_mpi_t);
		MPI_Type_commit(&PARTICLE_mpi_t);
	}
	
	
	// Prepare and launch sender and receiver threads //
	{
		if ((qid_receive = msgget(IPC_PRIVATE, 0666)) < 0) {
			fprintf(stderr, "Receive message queue cannot be created\n");
			return -1;
		}
		if ((qid_send = msgget(IPC_PRIVATE, 0666)) < 0) {
			fprintf(stderr, "Send message queue cannot be created\n");
			return -1;
		}
		
		if (pthread_create(&thread_receiver_tid, NULL,
						   thread_receiver, NULL) < 0) {
			fprintf(stderr, "Receiver cannot be launched\n");
			return -1;
		}
		if (pthread_create(&thread_sender_tid, NULL,
						   thread_sender, NULL) < 0) {
			fprintf(stderr, "Sender cannot be launched\n");
			return -1;
		}
		
		if (randgen_init(QUEUE_RAND_NORMAL_LENGTH) < 0) {
			fprintf(stderr, "Random number generator cannot be launched\n");
			return -1;
		}
	}
	
	
	// Prepare particles of this sector //
	{
		srand48((long)time(NULL));
		
		depot_a_length = (int)N;
		if ((depot_a = malloc(depot_a_length * sizeof(PARTICLE))) == NULL) {
			fprintf(stderr, "malloc() / realloc() error\n");
			return -1;
		}
		for (int i = 0; i < N; i++) {
			depot_a[i].origin_sector = mpi_rank;
			depot_a[i].steps_to_live = n;
			depot_a[i].x = lrand48() % l;
			depot_a[i].x += (mpi_rank % a) * l;
			depot_a[i].y = lrand48() % l;
			depot_a[i].y += (mpi_rank / a) * l;
		}
		
		depot_d_length = 0;
		depot_d = NULL;
	}
	
	
	// Run main cycle of the program //
	{
		struct msgbuf_particle particle_buffer;
		
		/// If set to nonzero, then main() should cease operations
		int stop_flag = 0;
		
		int tick = 0;
		while (stop_flag == 0) {
			tick = (tick + 1) % RECEIVE_TICK_INTERVAL;
			
			if (tick == 0) {
				// Receive all messages in the queue
				
				while (msgrcv(qid_receive, &particle_buffer, SIZEOF_MSGBUF_PARTICLE, 0, IPC_NOWAIT) >= 0) {
					if (particle_buffer.type == CLASSIFIER_STOP) {
						stop_flag = 1;
						break;
						// Note that the outer 'while' cycle is not broken.
						// This is correct as there may be some particles
						// send to this sector from another one
						// after their actual death.
						// They should be put into 'depot_d' of this sector.
					}
					
					depot_add(&particle_buffer.particle);
				}
				
				// Remove ENOMSG error
				errno = 0;
			}
			
			if (depot_a_length == 0) {
				// Require an immediate particle receive
				tick = -1;
				continue;
			}
			
			for (int i = 0; i < depot_a_length; i++) {
				if (depot_a[i].steps_to_live == 0) {
					if (get_sector(&depot_a[i]) == mpi_rank) {
						depot_die(i);
					}
					else {
						particle_buffer.particle = depot_a[i];
						particle_buffer.type = CLASSIFIER_NORM;
						msgsnd(qid_send,
							   &particle_buffer, SIZEOF_MSGBUF_PARTICLE, 0);
						depot_remove(i);
					}
				}
				else {
					/// Random number; see also 'randgen_get'
					double d_rand;
					
					/// Cumulative sum of p_*. Used to determine the way a particle should go
					double p_cursor = 0.0;
					
					while ((d_rand = randgen_get()) < -DOUBLE_ERROR) {}
					
					p_cursor += p_left;
					if (d_rand - p_cursor < DOUBLE_ERROR) {
						depot_a[i].x = (depot_a[i].x == 0) ?
									   (field_width - 1) :
									   (depot_a[i].x - 1);
						d_rand = 10.0;
					}
					
					p_cursor += p_up;
					if (d_rand - p_cursor < DOUBLE_ERROR) {
						depot_a[i].y = (depot_a[i].y == 0) ?
									   (field_height - 1) :
									   (depot_a[i].y - 1);
						d_rand = 10.0;
					}
					
					p_cursor += p_right;
					if (d_rand - p_cursor < DOUBLE_ERROR) {
						depot_a[i].x = (depot_a[i].x + 1) % field_width;
						d_rand = 10.0;
					}
					
					p_cursor += p_down;
					if (d_rand - p_cursor < DOUBLE_ERROR) {
						depot_a[i].y = (depot_a[i].y + 1) % field_height;
						// d_rand = 10.0;
					}
					
					depot_a[i].steps_to_live -= 1;
					
					if (get_sector_with_offset(&depot_a[i], mpi_rank) != mpi_rank) {
						particle_buffer.particle = depot_a[i];
						particle_buffer.type = CLASSIFIER_NORM;
						msgsnd(qid_send,
							   &particle_buffer, SIZEOF_MSGBUF_PARTICLE, 0);
						depot_remove(i);
					}
				}
			}
		}
		
		// Check 'depot_a'
		if (depot_a_length != 0) {
			fprintf(stderr, "A processing error occurred: Some particles may have stayed alive\n");
		}
	}
	
	
	// Completion actions //
	{
		pthread_join(thread_receiver_tid, NULL);
		pthread_join(thread_sender_tid, NULL);
		randgen_stop();
		
		struct timeval t_end;
		gettimeofday(&t_end, NULL);
		
		double elapsed_time = ((double) t_end.tv_sec +
							   ((double) t_end.tv_usec / 1000000.0) -
							   (double) t_start.tv_sec -
							   ((double) t_start.tv_usec / 1000000.0));
		
		if (mpi_rank == MPI_MASTER_RANK) {
			// Write statistics into file
			double elapsed_time_total = elapsed_time;
			for (int i = 1; i < mpi_size; i++) {
				MPI_Recv(&elapsed_time,
						 1, MPI_DOUBLE,
						 MPI_ANY_SOURCE, CLASSIFIER_TIME, MPI_COMM_WORLD, NULL);
				elapsed_time_total += elapsed_time;
			}
			
			MPI_File stats;
			MPI_File_open(MPI_COMM_SELF,
						  "stats.txt", MPI_MODE_APPEND | MPI_MODE_WRONLY,
						  MPI_INFO_NULL, &stats);
			
			for (int i = 0; i < argc; i++) {
				char buffer_to_print = ' ';
				MPI_File_write(stats,
							   argv[i], (int)strlen(argv[i]),
							   MPI_CHAR, MPI_STATUS_IGNORE);
				MPI_File_write(stats,
							   &buffer_to_print, 1,
							   MPI_CHAR, MPI_STATUS_IGNORE);
			}
			
			char buffer_to_print[3] = "s\n";
			MPI_File_write(stats,
						   &elapsed_time_total, 1,
						   MPI_DOUBLE, MPI_STATUS_IGNORE);
			MPI_File_write(stats,
						   &buffer_to_print, 2,
						   MPI_CHAR, MPI_STATUS_IGNORE);
			
			MPI_File_close(&stats);
		}
		else {
			MPI_Send(&elapsed_time, 1, MPI_DOUBLE, 0, CLASSIFIER_TIME, MPI_COMM_WORLD);
		}
	}
	
	
	// Output //
	{
		long particle_absolute_distribution[mpi_size];
		for (int i = 0; i < mpi_size; i++) {
			particle_absolute_distribution[i] = 0;
		}
		
		for (int i = 0; i < depot_d_length; i++) {
			particle_absolute_distribution[depot_d[i].origin_sector] += 1;
		}
		
		MPI_File datafile;
		MPI_File_open(MPI_COMM_WORLD,
					  "data.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY,
					  MPI_INFO_NULL, &datafile);
		
		MPI_Aint offset;
		MPI_Type_extent(MPI_LONG, &offset);
		
		MPI_File_write_at_all(datafile,
							  offset * mpi_rank * mpi_size,
							  particle_absolute_distribution,
							  mpi_size, MPI_LONG,
							  MPI_STATUS_IGNORE);
	}
	
	free(depot_a);
	free(depot_d);
	
	MPI_Finalize();
	return 0;
}



