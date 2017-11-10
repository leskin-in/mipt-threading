#include "lock_queue_2.h"




/// Queue node element
typedef struct Q_NODE_t_2 {
	struct Q_NODE_t_2* next;
	Q2_ELEMENT elem;
} Q2_NODE;




// Queue structures //

pthread_mutex_t Q2_mutex = PTHREAD_MUTEX_INITIALIZER;

Q2_NODE* Q2_head = NULL;

Q2_NODE* Q2_tail = NULL;




int queue2_push(const Q2_ELEMENT *elem) {
	pthread_mutex_lock(&Q2_mutex);
	
	Q2_NODE* new = malloc(sizeof(Q2_NODE));
	if (new == NULL) {
		pthread_mutex_unlock(&Q2_mutex);
		return -1;
	}
	new ->elem = *elem;
	new ->next = NULL;
	
	if (Q2_tail == NULL) {
		Q2_tail = new;
		Q2_head = new;
	}
	else {
		Q2_tail ->next = new;
		Q2_tail = new;
	}
	
	pthread_mutex_unlock(&Q2_mutex);
	return 0;
}




int queue2_pop(Q2_ELEMENT *result) {
	pthread_mutex_lock(&Q2_mutex);
	
	if (Q2_head == NULL) {
		pthread_mutex_unlock(&Q2_mutex);
		return -1;
	}
	
	Q2_NODE* old_ptr = Q2_head;
	*result = Q2_head ->elem;
	Q2_head = Q2_head ->next;
	
	if (Q2_head == NULL) {
		Q2_tail = NULL;
	}
	
	free(old_ptr);
	pthread_mutex_unlock(&Q2_mutex);
	return 0;
}




void queue2_clean(void) {
	Q2_ELEMENT dummy;
	while (queue2_pop(&dummy) != -1) {}
}
