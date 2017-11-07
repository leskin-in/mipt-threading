#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include <unistd.h>

#include "lock_queue.h"
#include "lock_queue_2.h"




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

const size_t SIZEOF_MSGBUF_PARTICLE = 4 * sizeof(long);

const double DOUBLE_ERROR = 1e-7;

/// Interval (in ticks) at which the receive operation is performed
const int RECEIVE_TICK_INTERVAL = 64;


int Mpi_size;

int Mpi_rank;


/// Global MPI_Datatype for the PARTICLE struct
MPI_Datatype PARTICLE_mpi_t;


pthread_t Thread_processor_tid;

// When this is locked, thread_processor cannot act
pthread_mutex_t Thread_processor_lock;


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
long Base_l;
/// Sector number at X axis
long Base_a;
/// Sector number at Y axis
long Base_b;
/// Total number of steps a particle should do
long Base_n;
/// Total number of particles
long Base_N;
/// Total number of dead particles
long N_dead;

double P_left;
double P_right;
double P_up;
double P_down;

/// Sector offset length
long Offset;

/// Maximum value of 'offset'
const long OFFSET_MAX = 64;

/// Divisor of 'offset' calculation. 'offset' = 'l' / 'OFFSET_DIVISOR'
const long OFFSET_DIVISOR = 5;

long Field_width;
long Field_height;

/// Depot for the active particles (currently processed)
PARTICLE* Depot_a;
/// Total number of the particles currently being processed
int Depot_a_length;
/// Total size of allocated memory for 'Depot_a'
int Depot_a_memory_length;

/// Depot for the dead particles (which died in current thread's action zone)
PARTICLE* Depot_d;
/// Total number of the dead particles in current thread's action zone
int Depot_d_length;
/// Total size of allocated memory for 'Depot_d'
int Depot_d_memory_length;




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
    return (unit ->y / Base_l) * Base_a + (unit ->x / Base_l);
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
    
    int field_borders_reached;
    
    long x_min = Base_l * (current_sector % Base_a);
    long y_min = Base_l * (current_sector / Base_a);
    long x_max = x_min + Base_l - 1;
    long y_max = y_min + Base_l - 1;
    x_min -= Offset;
    y_min -= Offset;
    x_max += Offset;
    y_max += Offset;
    
    field_borders_reached = 0;
    if (x_min < 0) {
        x_min += Field_width;
        field_borders_reached = 1;
    }
    else if (x_max >= Field_width) {
        x_max -= Field_width;
        field_borders_reached = 1;
    }
    if (field_borders_reached == 1) {
        if ((unit ->x < x_min) && (unit ->x > x_max)) {
            return ruled_sector;
        }
    }
    else {
        if ((unit ->x < x_min) || (unit ->x > x_max)) {
            return ruled_sector;
        }
    }
    
    field_borders_reached = 0;
    if (y_min < 0) {
        y_min += Field_height;
        field_borders_reached = 1;
    }
    else if (y_max >= Field_height) {
        y_max -= Field_height;
        field_borders_reached = 1;
    }
    if (field_borders_reached == 1) {
        if ((unit ->y < y_min) && (unit ->y > y_max)) {
            return ruled_sector;
        }
    }
    else {
        if ((unit ->y < y_min) || (unit ->y > y_max)) {
            return ruled_sector;
        }
    }
    
    return current_sector;
}



/**
 * @brief Add a patricle to 'depot_a'
 * @return 0 if successful
 */
int depot_add(const PARTICLE* unit) {
    Depot_a_length += 1;
    if (Depot_a_length > Depot_a_memory_length) {
        Depot_a_memory_length *= 2;
        Depot_a = realloc(Depot_a, Depot_a_memory_length * sizeof(PARTICLE));
    }
    Depot_a[Depot_a_length - 1] = *unit;
    return 0;
}



/**
 * @brief Drop a particle from 'depot_a'
 * @noexcept
 */
void depot_remove(int depot_a_id) {
    if ((depot_a_id < 0) || (depot_a_id >= Depot_a_length)) {
        return;
    }
    Depot_a[depot_a_id] = Depot_a[Depot_a_length - 1];
    Depot_a_length -= 1;
}



/**
 * @brief Kill a particle from 'depot_a' and place it to 'depot_d'
 * @param depot_a_id Sequence number of a particle in 'depot_a'
 * @return 0 if successful
 * @return -1 if the 'realloc()' failed or a given ID is incorrect
 */
int depot_die(int depot_a_id) {
    if ((depot_a_id < 0) || (depot_a_id >= Depot_a_length)) {
        return -1;
    }
    Depot_d_length += 1;
    if (Depot_d_length > Depot_d_memory_length) {
        Depot_d_memory_length *= 2;
        Depot_d = realloc(Depot_d, Depot_d_memory_length * sizeof(PARTICLE));
    }
    Depot_d[Depot_d_length - 1] = Depot_a[depot_a_id];
    depot_remove(depot_a_id);
    return 0;
}



/**
 * @brief Thread function for the particle processing thread
 */
void* thread_processor(void* context) {
    // Check launch
    {
        if (context != NULL) {
            fprintf(stderr, "Receiver thread is not launched properly\n");
            return NULL;
        }
        
        pthread_mutex_lock(&Thread_processor_lock);
    }
    
    struct msgbuf_particle particle_buffer;
    
    int tick = 0;
    while (1) {
        tick = (tick + 1) % RECEIVE_TICK_INTERVAL;
        
        if (tick == 0) {
            int stop_flag = 0;
            
            // Receive all messages in the queue
            while (queue2_pop(&particle_buffer) == 0) {
                if (particle_buffer.type == CLASSIFIER_STOP) {
                    stop_flag = 1;
                    break;
                }
                
                depot_add(&(particle_buffer.particle));
            }
            // Remove ENOMSG error
            errno = 0;
            
            // Break the outer 'while' if the inner one was broken
            if (stop_flag == 1) {
                break;
            }
        }
        
        if (Depot_a_length == 0) {
            // Perform an immediate particle receive
            tick = -1;
            continue;
        }
        
        for (int i = 0; i < Depot_a_length; i++) {
            if (Depot_a[i].steps_to_live == 0) {
                if (get_sector(&Depot_a[i]) == Mpi_rank) {
                    particle_buffer.particle = Depot_a[i];
                    particle_buffer.type = CLASSIFIER_REPORT;
                    queue_push(&particle_buffer);
                    depot_die(i);
                    i -= 1;
                }
                else {
                    particle_buffer.particle = Depot_a[i];
                    particle_buffer.type = CLASSIFIER_NORM;
                    queue_push(&particle_buffer);
                    depot_remove(i);
                    i -= 1;
                }
            }
            else {
                /// Random number; see also 'randgen_get'
                double d_rand;
                
                /// Cumulative sum of p_*. Used to determine the way a particle should go
                double p_cursor = 0.0;
                
                d_rand = drand48();
                
                p_cursor += P_left;
                if (d_rand - p_cursor < DOUBLE_ERROR) {
                    Depot_a[i].x = (Depot_a[i].x == 0) ?
                                   (Field_width - 1) :
                                   (Depot_a[i].x - 1);
                    d_rand = 10.0;
                }
                
                p_cursor += P_up;
                if (d_rand - p_cursor < DOUBLE_ERROR) {
                    Depot_a[i].y = (Depot_a[i].y == 0) ?
                                   (Field_height - 1) :
                                   (Depot_a[i].y - 1);
                    d_rand = 10.0;
                }
                
                p_cursor += P_right;
                if (d_rand - p_cursor < DOUBLE_ERROR) {
                    Depot_a[i].x = (Depot_a[i].x + 1) % Field_width;
                    d_rand = 10.0;
                }
                
                p_cursor += P_down;
                if (d_rand - p_cursor < DOUBLE_ERROR) {
                    Depot_a[i].y = (Depot_a[i].y + 1) % Field_height;
                    // d_rand = 10.0;
                }
                
                Depot_a[i].steps_to_live -= 1;
                
                if (get_sector_with_offset(&Depot_a[i], Mpi_rank) != Mpi_rank) {
                    particle_buffer.particle = Depot_a[i];
                    particle_buffer.type = CLASSIFIER_NORM;
                    queue_push(&particle_buffer);
                    depot_remove(i);
                    i -= 1;
                }
            }
        }
    }
    
    // Check 'Depot_a'
    if (Depot_a_length != 0) {
        fprintf(stderr, "A processing error occurred\n");
    }
    
    pthread_mutex_unlock(&Thread_processor_lock);
    
    sleep(5);
    
    return NULL;
}




int main(int argc, char** argv) {
    // Initialize MPI //
    {
        int mpi_thread_env_provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_env_provided);
        if (mpi_thread_env_provided == MPI_THREAD_SINGLE) {
            fprintf(stderr, "This program requires pthread support\n");
            return -1;
        }
        
        MPI_Comm_size(MPI_COMM_WORLD, &Mpi_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &Mpi_rank);
    }
    
    
    // Read arguments and initialize key variables //
    {
        // Check argc
        if (argc != 10) {
            print_help();
            return -1;
        }
        
        // Convert arguments
        Base_l = strtol(argv[1], NULL, 10);
        Base_a = strtol(argv[2], NULL, 10);
        Base_b = strtol(argv[3], NULL, 10);
        Base_n = strtol(argv[4], NULL, 10);
        Base_N = strtol(argv[5], NULL, 10);
        
        N_dead = 0;
        P_left = strtod(argv[6], NULL);
        P_right = strtod(argv[7], NULL);
        P_up = strtod(argv[8], NULL);
        P_down = strtod(argv[9], NULL);
        
        if (Base_a * Base_b != Mpi_size) {
            fprintf(stderr, "A number of MPI units is incorrect\n");
            return -1;
        }
        
        Offset = ((Base_l / OFFSET_DIVISOR) > OFFSET_MAX) ?
                 (OFFSET_MAX) :
                 (Base_l / OFFSET_DIVISOR);
        
        Field_width = Base_a * Base_l;
        Field_height = Base_b * Base_l;
    }
    
    struct timeval t_start;
    gettimeofday(&t_start, NULL);
    
    
    // Initialize PARTICE_mpi_t //
    {
        MPI_Datatype sequence[1] = {MPI_LONG};
        int sequence_amounts[1] = {4};
        MPI_Aint sequence_displacements[1];
        sequence_displacements[0] = 0;
        
        MPI_Type_create_struct(1,
                               sequence_amounts, sequence_displacements,
                               sequence,
                               &PARTICLE_mpi_t);
        
        MPI_Type_commit(&PARTICLE_mpi_t);
    }
    
    
    // Prepare and launch sender and receiver threads //
    {
        srand48(time(NULL) + Mpi_rank * 10);
        
        // Prevent 'thread_processor' from doing anything
        pthread_mutex_init(&Thread_processor_lock, NULL);
        pthread_mutex_lock(&Thread_processor_lock);
        
        if (pthread_create(&Thread_processor_tid, NULL,
                           thread_processor, NULL) < 0) {
            fprintf(stderr, "Processor thread cannot be launched\n");
            return -1;
        }
    }
    
    
    // Prepare particles of this sector //
    {
        srand48((long)time(NULL));
        
        Depot_a_length = (int)Base_N;
        Depot_a_memory_length = 1;
        while (Depot_a_memory_length < Depot_a_length) {
            Depot_a_memory_length *= 2;
        }
        if ((Depot_a = malloc(Depot_a_memory_length * sizeof(PARTICLE))) == NULL) {
            fprintf(stderr, "malloc() / realloc() error\n");
            return -1;
        }
        
        for (int i = 0; i < Base_N; i++) {
            Depot_a[i].origin_sector = Mpi_rank;
            Depot_a[i].steps_to_live = Base_n;
            Depot_a[i].x = lrand48() % Base_l;
            Depot_a[i].x += (Mpi_rank % Base_a) * Base_l;
            Depot_a[i].y = lrand48() % Base_l;
            Depot_a[i].y += (Mpi_rank / Base_a) * Base_l;
        }
        
        Depot_d_length = 0;
        Depot_d_memory_length = Depot_a_memory_length;
        if ((Depot_d = malloc(Depot_d_memory_length * sizeof(PARTICLE))) == NULL) {
            fprintf(stderr, "malloc() / realloc() error\n");
            free(Depot_a);
            return -1;
        }
    }
    
    
    // Run main cycle of the program. 'main()' works with MPI //
    {
        pthread_mutex_unlock(&Thread_processor_lock);
        
        PARTICLE mpi_recv_holder;
        struct msgbuf_particle buffer_recv;
        struct msgbuf_particle buffer_send;
        MPI_Request request_recv;
        
        // Pre-cycle MPI_Irecv for the 'request_recv' initialization
        {
            MPI_Irecv(&mpi_recv_holder,
                      1, PARTICLE_mpi_t, MPI_ANY_SOURCE, MPI_ANY_TAG,
                      MPI_COMM_WORLD, &request_recv);
        }
        
        while (1) {
            // Receive
            {
                MPI_Status status;
                int flag = 0;
                
                MPI_Test(&request_recv, &flag, &status);
                
                if (flag) {
                    buffer_recv.particle = mpi_recv_holder;
                    if (status.MPI_TAG == CLASSIFIER_REPORT) {
                        // Some particle has died
                        N_dead += 1;
                        if (N_dead == Base_N * Mpi_size) {
                            // Cease operations: Send message to 'thread_processor'
                            buffer_recv.type = CLASSIFIER_STOP;
                            queue2_push(&buffer_recv);
                            break;
                        }
                    }
                    else {
                        buffer_recv.type = CLASSIFIER_NORM;
                        queue2_push(&buffer_recv);
                    }
                    MPI_Irecv(&mpi_recv_holder,
                              1, PARTICLE_mpi_t, MPI_ANY_SOURCE, MPI_ANY_TAG,
                              MPI_COMM_WORLD, &request_recv);
                }
            }
            
            
            // Send
            {
                MPI_Request dummy;
                
                if (queue_pop(&buffer_send) < 0) {
                    continue;
                }
                
                if (buffer_send.type == CLASSIFIER_REPORT) {
                    // Send report to everyone (including self)
                    for (int i = 0; i < Mpi_size; i++) {
                        MPI_Isend(&(buffer_send.particle), 1, PARTICLE_mpi_t,
                                  i,
                                  CLASSIFIER_REPORT, MPI_COMM_WORLD, &dummy);
                        MPI_Request_free(&dummy);
                    }
                }
                else {
                    MPI_Isend(&(buffer_send.particle), 1, PARTICLE_mpi_t,
                              (int)get_sector(&(buffer_send.particle)),
                              CLASSIFIER_NORM, MPI_COMM_WORLD, &dummy);
                    MPI_Request_free(&dummy);
                }
            }
        }
    }
    struct timeval t_end;
    gettimeofday(&t_end, NULL);
    double elapsed_time = ((double) t_end.tv_sec +
                           ((double) t_end.tv_usec / 1000000.0) -
                           (double) t_start.tv_sec -
                           ((double) t_start.tv_usec / 1000000.0));
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    // Completion actions //
    {
        MPI_Type_free(&PARTICLE_mpi_t);
        pthread_join(Thread_processor_tid, NULL);
        queue_clean();
        queue2_clean();
        errno = 0;
        
        if (Mpi_rank == MPI_MASTER_RANK) {
            // Write statistics into file
            double elapsed_time_total = elapsed_time;
            for (int i = 1; i < Mpi_size; i++) {
                MPI_Recv(&elapsed_time,
                         1, MPI_DOUBLE,
                         MPI_ANY_SOURCE, CLASSIFIER_TIME, MPI_COMM_WORLD, NULL);
                elapsed_time_total += elapsed_time;
            }
            
            MPI_File stats;
            MPI_File_open(MPI_COMM_SELF,
                          "./stats.txt",
                          MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
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
            
            char double_to_print[64];
            sprintf(double_to_print, "%.3lfs\n", elapsed_time_total);
            MPI_File_write(stats,
                           double_to_print, (int)strlen(double_to_print),
                           MPI_CHAR, MPI_STATUS_IGNORE);
            
            MPI_File_close(&stats);
        }
        else {
            MPI_Send(&elapsed_time, 1, MPI_DOUBLE, 0, CLASSIFIER_TIME, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    // Output //
    {
        unsigned int particle_absolute_distribution[Mpi_size];
        for (int i = 0; i < Mpi_size; i++) {
            particle_absolute_distribution[i] = 0;
        }
        
        for (int i = 0; i < Depot_d_length; i++) {
            particle_absolute_distribution[Depot_d[i].origin_sector] += 1;
        }
        
        MPI_File datafile;
        MPI_File_open(MPI_COMM_WORLD,
                      "./data.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY,
                      MPI_INFO_NULL, &datafile);
        
        MPI_Aint dummy;
        MPI_Aint offset;
        MPI_Type_get_extent(MPI_UNSIGNED, &dummy, &offset);
        
        MPI_File_write_ordered(datafile,
                               particle_absolute_distribution,
                               Mpi_size, MPI_UNSIGNED,
                               MPI_STATUS_IGNORE);
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File_close(&datafile);
    }
    
    
    free(Depot_a);
    free(Depot_d);
    
    MPI_Finalize();
    return 0;
}



