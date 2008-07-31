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

vector_list* mcmclib_vector_list_transpose(vector_list* first) {
	int d = first->v->size;
	int n = mcmclib_vector_list_length(first);
	vector_list* head = mcmclib_vector_list_alloc();
	vector_list *last, *current;
	head->v = gsl_vector_alloc(n);
	last = head;
	for(int i=0; i<d; i++)
		last = mcmclib_vector_list_append(gsl_vector_alloc(n), last);

	int i=0;
	while(first->next) {
		current = head;
		int j=0;
		while(current->next) {
			gsl_vector_set(current->v, i, gsl_vector_get(first->v, j));
			current = current->next;
			j++;
		}
		first = first->next;
		i++;
	}

	return head;
}
