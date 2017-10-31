#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <mpi.h>
#include <pthread.h>




// Custom types //

/// Particle type
typedef struct particle_t {
	long x;
	long y;
	long origin_sector;
	long steps_to_live;
} PARTICLE;

/// Message buffer (for internal communications) for particle
struct msgbuf_particle {
	long type;
	PARTICLE particle;
};




// Global service variables //

const size_t SIZEOF_MSGBUF_PARTICLE = 5 * sizeof(long);

const double DOUBLE_ERROR = 1e-7;

/// Interval (in ticks) at which the receive operation is performed
const int RECEIVE_TICK_INTERVAL = 64;

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
#define CLASSIFIER_NORM 0

/// System V message type and IPC tag: cease operations
#define CLASSIFIER_STOP 1

/// IPC tag: particle death report
#define CLASSIFIER_REPORT 3




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
inline long get_sector(const PARTICLE* unit) {
	return (unit ->y / l) * a + (unit ->x / l);
}



/**
 * @param current_sector ID of the sector which holds 'unit' at the moment
 * @return Sector that should process 'unit'
 */
inline long get_sector_with_offset(const PARTICLE* unit, long current_sector) {
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
 * @return Position of 'unit' in 1-dimension space
 */
inline long get_through_position(const PARTICLE* unit) {
	return unit ->y * field_width + unit ->x;
}



/**
 * @brief Add a patricle to 'depot_a'
 * @return 0 if successful
 * @return -1 if there is not enough memory available for the 'realloc()'
 */
inline int depot_add(const PARTICLE* unit) {
	depot_a_length += 1;
	depot_a = realloc(depot_a, depot_a_length * sizeof(PARTICLE));
	if (errno != 0) {
		return -1;
	}
	depot_a[depot_a_length - 1] = *unit;
	return 0;
}



/**
 * @brief Remove a particle from 'depot_a'
 * @noerror
 */
inline void depot_remove(int depot_a_id) {
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
inline int depot_die(int depot_a_id) {
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



void* thread_rand_supplier(void* context) {
	// Check launch
	{
		if (context != NULL) {
			fprintf(stderr, "Random supplier thread is not launched properly\n");
			return NULL;
		}
	}
	
	
}








int main(int argc, char** argv) {
	// Initialize MPI //
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	
	
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
	}
	
	
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
	int tick = 0;
	while (1) {
		tick = (tick + 1) % RECEIVE_TICK_INTERVAL;
		
		if (tick == 0) {
			// Receive particles
		}
		
		
		break;
	}
	
	MPI_Finalize();
	return 0;
}




