typedef struct {
  gsl_vector* num; /*numerators*/
  gsl_vector* den; /*denumerators*/
  int B0; /*block length*/
  int t; /*current time step*/
  mcmclib_rapt* s;
} block_info;

block_info* block_info_alloc(mcmclib_rapt* s, int B0) {
  int K = s->K;
  block_info* ans = (block_info*) malloc(sizeof(block_info));
  ans->num = gsl_vector_alloc(K);
  ans->den = gsl_vector_alloc(K);
  ans->s = s;
  ans->B0 = B0;
  ans->t = 0;
  return ans;
}

void block_info_free(block_info* p) {
  gsl_vector_free(p->num);
  gsl_vector_free(p->den);
  free(p);
}

void block_info_update(block_info* p) {
  mcmclib_rapt* s = p->s;
  p->t = s->t % p->B0;
  if(p->t == 1) {
    gsl_vector_set_all(p->num, 0.0);
    gsl_vector_set_all(p->den, 0.0);
  }
  int r = s->which_region_x;
  gsl_vector_set(p->den, r, gsl_vector_get(p->den, r) + 1.0);
  if(! s->accepted)
    return;
  gsl_vector_set(p->num, r, gsl_vector_get(p->num, r) + 1);
}

double block_info_score(block_info* p) {
  double ans = 0.0;
  gsl_vector* num = p->num;
  gsl_vector* den = p->den;
  for(int k=0; k< p->s->K; k++)
    if(gsl_vector_get(den, k) > 0)
      ans += pow(gsl_vector_get(num, k) / gsl_vector_get(den, k) - 0.234, 2);
    else
      printf("still no visits to region %d!", k);
  return ans;
}
