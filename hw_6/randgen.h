#ifndef RANDGEN_H
#define RANDGEN_H

#include <errno.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/ipc.h>
#include <sys/msg.h>




/**
 * @brief Start randgen thread
 * @param buffer_size Total number of pre-generated numbers stored in memory
 * @param eseed Special value which will make random generator more random
 *
 * @return 0 if successful
 * @return -1 for errors
 */
int randgen_init(unsigned int buffer_size, unsigned int eseed);


/**
 * @return Pre-generated uniformly distributed t [0.0; 1.0] double random number
 * @return -1.0 if an error occured
 */
double randgen_get(void);


/**
 * @brief Stop randgen thread
 */
void randgen_stop(void);




#endif //RANDGEN_H
