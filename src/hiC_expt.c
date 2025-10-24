#include <Rinternals.h>
#include <R.h>
#include <pthread.h>

#define MAX_THREADS 96

// some experimental functions for looking at stats of hiC data


// To multithread the (otherwise) horrendously CPU intensive assessment
// we need a struct that contains suitable parameters. The simplest way
// to multithread this is for each thread to handle one column. 
// However, since the time taken to complete different columns can vary
// it makes sense to have some form of a cueing system. This can be implemented
// within the thread, where the thread locks a mutex, finds a job, and
// allocates that, and then unlocks the mutex.
typedef struct {
  size_t pos_n;
  int *left;
  int *right;
  // pos_n: the number of positions; i.e. the length of left and right
  // The thread will;
  // 1. lock the job_mutex
  // 2. find an available job (jobs[i] == 0)
  // 3. reserve the job (jobs[i] = 1)
  // 4. unlock the mutex
  pthread_mutex_t *job_mutex;
  size_t jobs_n;
  // note; there is one more break than jobs;
  // each job should run from breaks[i] -> breaks[i+1]
  int *jobs;
  double *breaks;
  // positions should be in the range:
  // min_pos <= pos < max_pos
  int max_left;
  double max_right_log;
  // values to be incremented; there should be nrow of them
  // these will be assumed to have been set to 0.
  int thinner;
  size_t nrow;
  double *obs_table;
  double *exp_table;
} thread_args;

thread_args init_thread_args(size_t pos_n, int *left, int *right,
			     pthread_mutex_t *job_mutex, size_t jobs_n,
			     int *jobs, double *breaks,
			     int max_left, double max_right_log,
			     int thinner, size_t nrow,
			     double *obs_table, double *exp_table){
  thread_args args;
  args.pos_n = pos_n;
  args.left = left;
  args.right = right;
  args.job_mutex = job_mutex;
  args.jobs_n = jobs_n;
  args.jobs = jobs;
  args.breaks = breaks;
  args.max_left = max_left;
  args.max_right_log = max_right_log;
  args.thinner = thinner;
  args.nrow = nrow;
  args.obs_table = obs_table;
  args.exp_table = exp_table;
  return(args);
}

void *assess_column(void *v_args){
  thread_args *args = (thread_args*)v_args;
  size_t job_i = 0;
  size_t i = 0; // should only be incremented from here on
  size_t pos_n = args->pos_n;
  int *left = args->left;
  int *right = args->right;
  int thinner = args->thinner;
  double max_right_log = args->max_right_log;
  // first select a job to do;
  while(job_i < args->jobs_n){
    pthread_mutex_lock(args->job_mutex);
    while( job_i < args->jobs_n && args->jobs[job_i] )
      ++job_i;
    if(job_i < args->jobs_n)
      args->jobs[job_i] = 1;
    pthread_mutex_unlock(args->job_mutex);
    if(job_i >= args->jobs_n)
      break;
    double min_pos = args->breaks[job_i];
    double max_pos = args->breaks[job_i+1];
    size_t nrow = args->nrow;
    double *obs = args->obs_table + job_i * nrow;
    double *exp = args->exp_table + job_i * nrow;
    while(i < pos_n && (double)left[i] < min_pos)
      ++i;
    while(i < pos_n && (double)left[i] < max_pos){
      size_t row = (size_t)(nrow * log((double)(right[i] - left[i])) / max_right_log);
      obs[row]++;
      for(size_t j=i+thinner; j < pos_n; j += thinner){
	double d_l = (double)( left[j] - left[i] );
	double d_r = (double)( right[j] - left[i] );
	exp[(size_t)(nrow * log(d_l) / max_right_log)] += 1/d_l;
	exp[(size_t)(nrow * log(d_r) / max_right_log)] += 1/d_r;
      }
      ++i;
    }
  }
  pthread_exit((void*)args);
}

// Given a matrix with columns:
// 1. left position
// 2. right position
// And 
// 3. h_breaks_n: no of linear horizontal breaks (left position)
// 4. v_breaks_n : vertical breaks, exponential (right - left)
//
// the function will count the number of links (left_position -> right_position)
// found; and estimate the expected number given the the number of positions
// lying to the right of each left position.
//
// The expected proportion of links in each cell will be calculated as:
// sum[ (1 / (right[(i+1):n]-left[i]) ) ]
//
// \sum_{j=i+1..n}\over{1}{ r_j - l_i }
//
// where r is the right position and l the left position.
// both arguments should be integers,
// pos_r should be a matrix of length n
// breaks_n_r should be a vector of length 3
//       giving horizontal breaks, vertical breaks and a thinning
//       parameter to save time..
//
SEXP assess_region_links(SEXP pos_r, SEXP breaks_n_r){
  if(TYPEOF(pos_r) != INTSXP || TYPEOF(breaks_n_r) != INTSXP)
    error("Both arguments should be integer vectors");
  SEXP pos_dr = PROTECT( getAttrib( pos_r, R_DimSymbol ));
  if(length(pos_dr) != 2){
    UNPROTECT(1);
    error("pos_r should be a matrix");
  }
  int nrow = INTEGER(pos_dr)[0];
  int ncol = INTEGER(pos_dr)[1];
  UNPROTECT(1);
  if(length(breaks_n_r) != 3)
    error("breaks_n_r should have three elements");
  int h_breaks_n = INTEGER(breaks_n_r)[0];
  int v_breaks_n = INTEGER(breaks_n_r)[1];
  int thinner = INTEGER(breaks_n_r)[2];
  // Sanity check the sizes of elements; should specify the thresholds
  // somewhere else.
  if(nrow < 2 * h_breaks_n || nrow < 2 * v_breaks_n)
    error("not enough rows");
  if(ncol != 2 || h_breaks_n < 3 || v_breaks_n < 3 || thinner <= 0 || nrow < thinner)
    error("unreasonable parameters");

  int *left = INTEGER(pos_r);
  int *right = left + nrow;

  // the integer is for bounds checking without having to cast
  int max_right = right[0];
  int max_left = left[0];
  for(int i=1; i < nrow; ++i){
    if(left[i] < left[i-1] || left[i] >= right[i])
      error("left positions must be sorted and right positions must be > left");
    max_right = max_right < right[i] ? right[i] : max_right;
    max_left = max_left < left[i] ? left[i] : max_left;
  }
    
  double max_right_log = log((double)max_right);
  
  // The return data should be two tables;
  // 1. one of counts,
  // 2. expected frequencies
  // In addition to the breaks of those tables.
  size_t nr = (size_t)v_breaks_n - 1;
  size_t nc = (size_t)h_breaks_n - 1;
  SEXP ret_data = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ret_data, 0, allocMatrix(REALSXP, nr, nc));
  SET_VECTOR_ELT(ret_data, 1, allocMatrix(REALSXP, nr, nc));
  SET_VECTOR_ELT(ret_data, 2, allocVector(REALSXP, h_breaks_n));
  SET_VECTOR_ELT(ret_data, 3, allocVector(REALSXP, v_breaks_n));
  double *counts = REAL(VECTOR_ELT(ret_data, 0));
  double *expected = REAL(VECTOR_ELT(ret_data, 1));
  double *h_breaks = REAL(VECTOR_ELT(ret_data, 2));
  double *v_breaks = REAL(VECTOR_ELT(ret_data, 3));
  memset( counts, 0, sizeof(double) * nr * nc );
  memset( expected, 0, sizeof(double) * nr * nc);

  for(int i=0; i < h_breaks_n; ++i)
    h_breaks[i] = ((double)i * max_left) / nc;
  for(int i=0; i < v_breaks_n; ++i)
    v_breaks[i] = (double)i * max_right_log / nr;
  
  Rprintf("max pos: %d  max pos_log %lf  nr: %ld  nc: %ld thinner: %d\n", max_right, max_right_log, nr, nc, thinner);
  
  for(int i=0; i < nrow; ++i){
    if(left[i] >= max_right)
      break;
    size_t k = (nc * (size_t)left[i] / max_left);
    size_t row = (size_t)(nr * log((double)right[i] - left[i]) / max_right_log );
    double *counts_c = counts + k * nr;
    double *expected_c = expected + k * nr;
    counts_c[row]++;
    /* Rprintf("left: %d  right: %d  delta: %d  max_right: %d  max_right_log %lf  k: %ld row: %ld\n", */
    /* 	    left[i], right[i], right[i] - left[i], max_right, max_right_log, k, row); */
    for(int j=i + thinner; j < nrow; j += thinner){
      // I am not sure as to whether we should considerin left[j] - left[i]
      // or right[j] - left[i]
      // We are trying to estimate the density of mapped reads at different distances from
      // left[i]; it should not really matter whether this are left or right, so the best
      // thing may be to consider both.
      double d_l = (double)( left[j] - left[i] );
      double d_r = (double)( right[j] - left[i] );
      expected_c[(size_t)(nr * log(d_l) / max_right_log)] += 1/d_l;
      expected_c[(size_t)(nr * log(d_r) / max_right_log)] += 1/d_r;
    }
  }
  UNPROTECT(1);
  return(ret_data);
}

/// temporarily copied and pasted; eventually keep only the multithreaded
/// version
/// arguments as above; except breaks_n_r should have one additional
/// argument giving the number of threads;
SEXP assess_region_links_mt(SEXP pos_r, SEXP breaks_n_r){
  if(TYPEOF(pos_r) != INTSXP || TYPEOF(breaks_n_r) != INTSXP)
    error("Both arguments should be integer vectors");
  SEXP pos_dr = PROTECT( getAttrib( pos_r, R_DimSymbol ));
  if(length(pos_dr) != 2){
    UNPROTECT(1);
    error("pos_r should be a matrix");
  }
  int nrow = INTEGER(pos_dr)[0];
  int ncol = INTEGER(pos_dr)[1];
  UNPROTECT(1);
  if(length(breaks_n_r) != 4)
    error("breaks_n_r should have four elements");
  int h_breaks_n = INTEGER(breaks_n_r)[0];
  int v_breaks_n = INTEGER(breaks_n_r)[1];
  int thinner = INTEGER(breaks_n_r)[2];
  int nthreads = INTEGER(breaks_n_r)[3];
  // Sanity check the sizes of elements; should specify the thresholds
  // somewhere else.
  if(nrow < 2 * h_breaks_n || nrow < 2 * v_breaks_n)
    error("not enough rows");
  if(ncol != 2 || h_breaks_n < 3 || v_breaks_n < 3 || thinner <= 0 || nrow < thinner)
    error("unreasonable parameters");
  if(nthreads <= 0)
    nthreads = 1;
  if(nthreads > MAX_THREADS)
    nthreads = MAX_THREADS;
  // but given the design, we cannot have fewer threads than columns of the results
  // table;
  if(nthreads > (h_breaks_n-1)){
    warning("more threads than columns requested; nthreads reduced to column number");
    nthreads = h_breaks_n - 1;
  }
  
  int *left = INTEGER(pos_r);
  int *right = left + nrow;

  // the integer is for bounds checking without having to cast
  int max_right = right[0];
  int max_left = left[0];
  for(int i=1; i < nrow; ++i){
    if(left[i] < left[i-1] || left[i] >= right[i])
      error("left positions must be sorted and right positions must be > left");
    max_right = max_right < right[i] ? right[i] : max_right;
    max_left = max_left < left[i] ? left[i] : max_left;
  }
    
  double max_right_log = log((double)max_right);
  
  // The return data should be two tables;
  // 1. one of counts,
  // 2. expected frequencies
  // In addition to the breaks of those tables.
  int nr = v_breaks_n - 1;
  int nc = h_breaks_n - 1;
  SEXP ret_data = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ret_data, 0, allocMatrix(REALSXP, nr, nc));
  SET_VECTOR_ELT(ret_data, 1, allocMatrix(REALSXP, nr, nc));
  SET_VECTOR_ELT(ret_data, 2, allocVector(REALSXP, h_breaks_n));
  SET_VECTOR_ELT(ret_data, 3, allocVector(REALSXP, v_breaks_n));
  double *counts = REAL(VECTOR_ELT(ret_data, 0));
  double *expected = REAL(VECTOR_ELT(ret_data, 1));
  double *h_breaks = REAL(VECTOR_ELT(ret_data, 2));
  double *v_breaks = REAL(VECTOR_ELT(ret_data, 3));
  memset( counts, 0, sizeof(double) * nr * nc );
  memset( expected, 0, sizeof(double) * nr * nc);

  for(int i=0; i < h_breaks_n; ++i)
    h_breaks[i] = ((double)i * max_left) / nc;
  for(int i=0; i < v_breaks_n; ++i)
    v_breaks[i] = (double)i * max_right_log / nr;
  
  pthread_mutex_t job_mutex = PTHREAD_MUTEX_INITIALIZER;
  size_t jobs_n = (size_t)nc;  // one job per column
  int *thread_jobs = calloc(jobs_n, sizeof(int));
  
  thread_args t_args = init_thread_args(nrow, left, right, &job_mutex, jobs_n, thread_jobs,
					h_breaks, max_left, max_right_log, thinner,
					nr, counts, expected);

  pthread_t *threads = malloc( sizeof(pthread_t) * nthreads );
  for(int i=0; i < nthreads; ++i){
    pthread_create( &threads[i], NULL, &assess_column, (void*)&t_args );
  }
  void *status;
  for(int i=0; i < nthreads; ++i)
    pthread_join(threads[i], &status);
  free(thread_jobs);
  free(threads);
  
  UNPROTECT(1);
  return(ret_data);
}


static const R_CallMethodDef callMethods[] =
  {
   {"assess_region_links", (DL_FUNC)&assess_region_links, 2},
   {"assess_region_links_mt", (DL_FUNC)&assess_region_links_mt, 2},
   {NULL, NULL, 0}
  };

void R_init_hiC_expt(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

