#include "vector_list.h"

vector_list* mcmclib_vector_list_alloc() {
	vector_list* ans = (vector_list*) malloc(sizeof(vector_list));
	ans->next = NULL;
	return ans;
}

vector_list* mcmclib_vector_list_last(vector_list* i) {
	while(i->next) i = i->next;
	return i;
}

vector_list* mcmclib_vector_list_append(gsl_vector* v, vector_list* element) {
	vector_list* item = (vector_list*) malloc(sizeof(vector_list));
	item->v = v;
	item->next = NULL;
	mcmclib_vector_list_last(element)->next = item;
	return item;
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
