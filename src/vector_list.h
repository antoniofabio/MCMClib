#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

struct vector_list_str {
	gsl_vector* v;
	struct vector_list_str* next;
};

/** List of gsl_vector objects
*/
typedef struct vector_list_str vector_list;

/** alloc a new vector_list object
*/
vector_list* mcmclib_vector_list_alloc();

/** get last list element
*/
vector_list* mcmclib_vector_list_last(vector_list* i);

/** append a new vector to the list
@param v pointer to the vector to be added
@param last pointer to last list element
@return the new appended list element
*/
vector_list* mcmclib_vector_list_append(gsl_vector* v, vector_list* last);

/** compute list length
@param first	pointer to list head
*/
int mcmclib_vector_list_length(vector_list* first);

/** de-allocates the full list from the heap
@param first pointer to list head
*/
void mcmclib_vector_list_free(vector_list* first);

/** 'transpose' a list of vectors
*/
vector_list* mcmclib_vector_list_transpose(vector_list* first);

/** convert a list of vectors into a matrix
each list element is interpreted as a row
*/
gsl_matrix* mcmclib_vector_list_asmatrix(vector_list* first);
