#include "vector_list.h"

void mcmclib_vector_list_add(gsl_vector* v, vector_list* last) {
	vector_list* item = (vector_list*) malloc(sizeof(vector_list));
	item->v = v;
	item->next = NULL;
	last->next = item;
}

int mcmclib_vector_list_length(vector_list* first) {
	int n = 0;
	while(first) {
		first = first->next;
		n++;
	}
	return n;
}

void mcmclib_vector_list_free(vector_list* first) {
	vector_list* save;
	while(first) {
		save = first;
		first = first->next;
		gsl_vector_free(save->v);
		free(save);
	}
}
