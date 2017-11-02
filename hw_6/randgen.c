#include "randgen.h"




// Custom structures //

/**
 * @brief Message buffer for random numbers
 */
struct randgen_msgbuf {
	long type;  /// Message type (always set to 0)
	double number;  /// Uniformly distributed at [0.0; 1.0) random number
};




// Global variables //

int randgen_buffer_size = -1;

/// System V queue identifier
int randgen_qid;

pthread_t randgen_supplier_tid;

pthread_mutex_t randgen_supplier_lock;

int randgen_supplier_status;




void* randgen_supplier(void* context) {
	struct msqid_ds qid_rand_supply_info;
	struct randgen_msgbuf buffer;
	buffer.type = 1;
	
	srand48((long)time(NULL));
	
	while (1) {
		pthread_mutex_lock(&randgen_supplier_lock);
		if (randgen_supplier_status == -1) {
			pthread_mutex_unlock(&randgen_supplier_lock);
			break;
		}
		pthread_mutex_unlock(&randgen_supplier_lock);
		
		msgctl(randgen_qid, IPC_STAT, &qid_rand_supply_info);
		for (int i = 0;
			 i < randgen_buffer_size - qid_rand_supply_info.msg_qnum;
			 i++) {
			buffer.number = drand48();
		}
	}
	
	return NULL;
}




int randgen_init(unsigned int buffer_size) {
	if ((randgen_buffer_size != -1) || (buffer_size == 0)) {
		return -1;
	}
	
	randgen_buffer_size = (int)buffer_size;
	
	if ((randgen_qid = msgget(IPC_PRIVATE, 0666)) < 0) {
		randgen_buffer_size = -1;
		return -1;
	}
	
	randgen_supplier_status = 0;
	
	pthread_mutex_init(&randgen_supplier_lock, NULL);
	
	if ((pthread_create(&randgen_supplier_tid, NULL, randgen_supplier, NULL) != 0)) {
		msgctl(randgen_qid, IPC_RMID, NULL);
		randgen_buffer_size = -1;
		return -1;
	}
	
	return 0;
}




double randgen_get(void) {
	if (randgen_buffer_size == -1) {
		return -1.0;
	}
	struct randgen_msgbuf buffer;
	if (msgrcv(randgen_qid, &buffer, sizeof(struct randgen_msgbuf), 0, 0) < 0) {
		return -1.0;
	}
	return buffer.number;
}




void randgen_stop(void) {
	if (randgen_buffer_size == -1) {
		return;
	}
	
	pthread_mutex_lock(&randgen_supplier_lock);
	randgen_supplier_status = -1;
	pthread_mutex_unlock(&randgen_supplier_lock);
	
	pthread_join(randgen_supplier_tid, NULL);
	msgctl(randgen_qid, IPC_RMID, NULL);
	
	randgen_buffer_size = -1;
}



