#ifndef Q2_H
#define Q2_H

#include <stdlib.h>
#include <pthread.h>




// Insert your element type here
struct queue_preelement_2 {
	long x;
	long y;
	long origin_sector;
	long steps_to_live;
};
typedef struct {
	long type;
	struct queue_preelement_2 preelement;
} Q2_ELEMENT;




/**
 * @brief Insert an element into a queue
 * @param elem
 * @return 0 if successful
 * @return -1 for errors
 */
int queue2_push(const Q2_ELEMENT *elem);




/**
 * @brief Get an element from a queue
 * @param result
 * @return 0 if an element was retrieved
 * @return -1 if the queue is empty
 */
int queue2_pop(Q2_ELEMENT *result);




/**
 * @brief Remove all the elements from the queue
 */
void queue2_clean(void);




#endif // RANDGEN_H
