#ifndef Q_H
#define Q_H

#include <stdlib.h>
#include <pthread.h>




// Insert your element type here
struct queue_preelement {
	long x;
	long y;
	long origin_sector;
	long steps_to_live;
};
typedef struct {
	long type;
	struct queue_preelement preelement;
} Q_ELEMENT;




/**
 * @brief Insert an element into a queue
 * @param elem
 * @return 0 if successful
 * @return -1 for errors
 */
int queue_push(const Q_ELEMENT* elem);




/**
 * @brief Get an element from a queue
 * @param result
 * @return 0 if an element was retrieved
 * @return -1 if the queue is empty
 */
int queue_pop(Q_ELEMENT* result);




/**
 * @brief Remove all the elements from the queue
 */
void queue_clean(void);




#endif // RANDGEN_H
