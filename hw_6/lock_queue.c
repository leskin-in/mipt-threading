#include "lock_queue.h"




/// Queue node element
typedef struct Q_NODE_t {
	struct Q_NODE_t* next;
	Q_ELEMENT elem;
} Q_NODE;




// Queue structures //

pthread_mutex_t Q_mutex = PTHREAD_MUTEX_INITIALIZER;

Q_NODE* Q_head = NULL;

Q_NODE* Q_tail = NULL;




int queue_push(const Q_ELEMENT* elem) {
	pthread_mutex_lock(&Q_mutex);
	
	Q_NODE* new = malloc(sizeof(Q_NODE));
	if (new == NULL) {
		pthread_mutex_unlock(&Q_mutex);
		return -1;
	}
	new ->elem = *elem;
	new ->next = NULL;
	
	if (Q_tail == NULL) {
		Q_tail = new;
		Q_head = new;
	}
	else {
		Q_tail ->next = new;
		Q_tail = new;
	}
	
	pthread_mutex_unlock(&Q_mutex);
	return 0;
}




int queue_pop(Q_ELEMENT* result) {
	pthread_mutex_lock(&Q_mutex);
	
	if (Q_head == NULL) {
		pthread_mutex_unlock(&Q_mutex);
		return -1;
	}
	
	Q_NODE* old_ptr = Q_head;
	*result = Q_head ->elem;
	Q_head = Q_head ->next;
	
	if (Q_head == NULL) {
		Q_tail = NULL;
	}
	
	free(old_ptr);
	pthread_mutex_unlock(&Q_mutex);
	return 0;
}




void queue_clean(void) {
	Q_ELEMENT dummy;
	while (queue_pop(&dummy) != -1) {}
}
