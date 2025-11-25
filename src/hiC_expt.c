#include <Rinternals.h>
#include <R.h>
#include <pthread.h>
#include "kvec.h"
#include "ksort.h"

#define MAX_THREADS 96

// some experimental functions for looking at stats of hiC data

// minimal information about an alignment pair;
typedef struct {
  int i;
  int tg;
  int pos;
  int m_tg;
  int m_pos;
} align_pair;

#define align_pair_lt(a, b) ( (a).pos < (b).pos )
#define align_pair_gt(a, b) ( (a).pos > (b).pos )
KSORT_INIT(pair_sort, align_pair, align_pair_lt);
KSORT_INIT(pair_sort_rev, align_pair, align_pair_gt);

align_pair init_al_pair(int i, int tg, int pos, int m_tg, int m_pos){
  align_pair ap;
  ap.i = i;
  ap.tg = tg;
  ap.pos = pos;
  ap.m_tg = m_tg;
  ap.m_pos = m_pos;
  return(ap);
}

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

#define IN_REGION(pos_id, pos, reg_id, reg_beg, reg_end) ( ((pos_id) == (reg_id)) && ((pos) >= (reg_beg)) && ((pos) < (reg_end)) )
#define IN_RANGE(pos, reg_beg, reg_end) ( ((pos) >= (reg_beg)) && ((pos) < (reg_end)) )
// IS_LEFT returns true if reg_id and m_reg (mate region) are different OR if pos < m_pos
// 
#define IS_LEFT(reg_id, m_reg, pos, m_pos) ( ((reg_id) != (m_reg)) || ((pos) < (m_pos)) )
// Like IS_LEFT, but pos > m_pos
#define IS_RIGHT(reg_id, m_reg, pos, m_pos) ( ((reg_id) != (m_reg)) || ((pos) > (m_pos)) )



// assess_separation of two regions
// Given read pair information for two regions (A,B)  provided as two matrices with
// pairs in rows and columns:
// index, pos, m.tg (mate target), m.pos (mate.pos)
//
// The regions to be assessed: This should be an integer vector giving
// 1. The index of region A (as used in the table; may be 0 based)
// 2. A min, A max, rotate A (0/1)
// 3. The index of region B
// 4. B min, B max, rotate B (0/1)
// 5. The minimum separation between adjacent reads
// 6. The strategy for selecting alignments:
//    0: Only alignments where the mate pairs are internal to either
//       or span the region and where the read IS_LEFT or IS_RIGHT
//       for region A and B respectively.
//    1: Any alignments that are present in either region.
// seps_r: a series of region separations to model
//
// The proportion of A internal to A-B spanning links and their positions
// will be determined as the sum( 1/d ) for all potential mate pairs
// for the separations indicated by seps_r
SEXP assess_separation(SEXP A_r, SEXP B_r, SEXP reg_limits_r, SEXP seps_r){
  if(TYPEOF(A_r) != INTSXP || TYPEOF(B_r) != INTSXP || TYPEOF(reg_limits_r) != INTSXP
     || TYPEOF(seps_r) != INTSXP)
    error("All arguments should be integer structures");
  const int regl_n = 10;
  if(length(reg_limits_r) != regl_n)
    error("The region information (reg_limits_r) must have %d elements", regl_n);
  int *reg_limits = INTEGER(reg_limits_r);
  // A_tg, A_min, A_max, A_rot, B_tg, B_min, B_rot;
  int A_tg = reg_limits[0]; int B_tg = reg_limits[4];
  int A_min = reg_limits[1]; int A_max = reg_limits[2]; int A_rot = reg_limits[3];
  int B_min = reg_limits[5]; int B_max = reg_limits[6]; int B_rot = reg_limits[7];
  int min_sep = reg_limits[8];
  int strategy = reg_limits[9];
  // sanity check;
  if(A_tg < 0 || B_tg < 0 || A_min < 0 || A_max <= A_min || B_min < 0 || B_max <= B_min)
    error("Unreasonable limits information given. Try again");
  if(strategy != 0)
    error("Only strategy 0 allowed");
    
  int seps_n = length(seps_r);
  if(seps_n < 1)
    error("Specify at least one separation distance");
  int *seps = INTEGER(seps_r);
    
  SEXP A_dims_r = PROTECT( getAttrib(A_r, R_DimSymbol) );
  SEXP B_dims_r = PROTECT( getAttrib(B_r, R_DimSymbol) );
  if(length(A_dims_r) != 2 || length(B_dims_r) != 2){
    UNPROTECT(2);
    error("Both A_r and B_r must be matrices");
  }
  int *A_dim = INTEGER(A_dims_r);
  int *B_dim = INTEGER(B_dims_r);
  if(A_dim[1] != 4 || B_dim[1] != 4 || A_dim[0] < 1 || A_dim[0] < 1){
    UNPROTECT(2);
    error("Both A_r and B_r must have four columns and some rows");
  }
  // for convenience define pointers into the matrix:
  int *A = INTEGER(A_r);
  int *B = INTEGER(B_r);
  int *A_pos = A + A_dim[0]; int *A_m_tg = A + A_dim[0] * 2; int *A_m_pos = A + A_dim[0] * 3;
  int *B_pos = B + B_dim[0]; int *B_m_tg = B + B_dim[0] * 2; int *B_m_pos = B + B_dim[0] * 3;
  
  // We will consider B to be to the right of A;
  // To reduce the complexity of the procedure we will consider only:
  // left reads on A that have a mate in A or B
  // right reads on B that have a mate in B or A
  //
  // This should be reasonable as we are concerned with the fraction of
  // reads that are either internal or that bridge the two regions

  // First determine the set that we keep
  // integers so that we can pass back to R
  // keep the indices in kvectors
  /* kvec_t(int) A_i;  kv_init(A_i); */
  /* kvec_t(int) B_i;  kv_init(B_i); */
    // copy the alignment information into vectors of align_pair structs
  // change coordinate system such that we use 0 -> distance for each one
  // This is in order to allow changing of the order if rot is true;
  kvec_t(align_pair) A_pairs; kv_init(A_pairs);
  kvec_t(align_pair) B_pairs; kv_init(B_pairs);

  // if the region is rotated iterate backwards through the region
  // but 
  int a_i[3] = {0, A_dim[0]-1, 1};
  int b_i[3] = {0, B_dim[0]-1, 1};

  if(A_rot)
    a_i[2] = -1;
  if(B_rot)
    b_i[2] = -1;

  // The positions in B should be a continuum with those in A; that means
  // that they should start from A_max+1 and run upwards
  int B_off = 1 + A_max-A_min;
  
  // keep left reads within the regions
  // This should clearly be refactored to a function of some sort.
  for(int i=a_i[A_rot]; i >= a_i[0] && i <= a_i[1]; i += a_i[2]){
    int in_range_A = IN_RANGE(A_pos[i], A_min, A_max);
    int m_in_range_A = IN_REGION(A_m_tg[i], A_m_pos[i], A_tg, A_min, A_max);
    int m_in_range_B = IN_REGION(A_m_tg[i], A_m_pos[i], B_tg, B_min, B_max);
    int is_left = IS_LEFT(A_tg, A_m_tg[i], A_pos[i], A_m_pos[i]);
    if(!in_range_A || !(m_in_range_A || m_in_range_B))
      continue;

    switch(A_rot){
    case 0:
      if(is_left)
	kv_push(align_pair, A_pairs,
		m_in_range_A ? 
		init_al_pair(i, A_tg, A_pos[i]-A_min, A_m_tg[i], A_m_pos[i]-A_min) :
		B_rot ?
		init_al_pair(i, A_tg, A_pos[i]-A_min, A_m_tg[i], B_off + B_max-A_m_pos[i]) :
		init_al_pair(i, A_tg, A_pos[i]-A_min, A_m_tg[i], B_off + A_m_pos[i]-B_min));
      break;
    case 1:
      if(!is_left)
	break;
      // this is more difficult as the mate may be in the other region
      // we should swap the positions over as well.
      if(m_in_range_A){     // easy
	kv_push(align_pair, A_pairs, init_al_pair(i, A_m_tg[i], A_max - A_m_pos[i], A_tg, A_max - A_pos[i]));
	//	kv_push(align_pair, A_pairs, init_al_pair(i, A_tg, A_max - A_pos[i], A_m_tg[i], A_max - A_m_pos[i]));
      }else{
	// 
	kv_push(align_pair, A_pairs, B_rot ?
		init_al_pair(i, A_tg, A_max - A_pos[i], A_m_tg[i], B_off + B_max - A_m_pos[i]) :
		init_al_pair(i, A_tg, A_max - A_pos[i], A_m_tg[i], B_off + A_m_pos[i] - B_min));
      }
      break;
    default :
      warning("A_rot has a value of %d. Should be between 0 and 1", A_rot);
    }
  }
  // REFACTOR requirement! This is BAD
  for(int i=b_i[B_rot]; i >= b_i[0] && i <= b_i[1]; i += b_i[2]){
    int in_range_B = IN_RANGE(B_pos[i], B_min, B_max);
    int m_in_range_A = IN_REGION(B_m_tg[i], B_m_pos[i], A_tg, A_min, A_max);
    int m_in_range_B = IN_REGION(B_m_tg[i], B_m_pos[i], B_tg, B_min, B_max);
    int is_right = IS_RIGHT(B_tg, B_m_tg[i], B_pos[i], B_m_pos[i]);
    if(!in_range_B || !(m_in_range_A || m_in_range_B))
      continue;

    switch(B_rot){
    case 0:
      if(is_right)
	kv_push(align_pair, B_pairs,
		m_in_range_B ?
		init_al_pair(i, B_tg, B_off + B_pos[i]-B_min, B_m_tg[i], B_off + B_m_pos[i]-B_min) :
		A_rot ?
		init_al_pair(i, B_tg, B_off + B_pos[i]-B_min, B_m_tg[i], A_max-B_m_pos[i]) :
		init_al_pair(i, B_tg, B_off + B_pos[i]-B_min, B_m_tg[i], B_m_pos[i]-A_min));
      break;
    case 1:
      if(!is_right)
	break;
      // this is more difficult as the mate may be in the other region
      if(m_in_range_B){     // easy
	kv_push(align_pair, B_pairs, init_al_pair(i, B_m_tg[i], B_off + B_max - B_m_pos[i], B_tg, B_off + B_max - B_pos[i]));
      }else{
	// The mate is in A, and may also need to be rotated
	kv_push(align_pair, B_pairs, A_rot ?
		init_al_pair(i, B_tg, B_off + B_max - B_pos[i], B_m_tg[i], A_max - B_m_pos[i]) :
		init_al_pair(i, B_tg, B_off + B_max - B_pos[i], B_m_tg[i], B_m_pos[i] - A_min));
      }
      break;
    default :
      warning("B_rot has a value of %d. Should be between 0 and 1", B_rot);
    }
  }
  // if we have rotated either A or B, then we need to sort again as we
  // do not have any guarantees of order. I'm not completely sure right now
  // if this is needed but try doing it first.
  if(A_rot || B_rot){
    ks_introsort( pair_sort, A_pairs.n, A_pairs.a );
    ks_introsort( pair_sort, B_pairs.n, B_pairs.a );
  }
  
  // For every left read (i) in A:
  // 1. Increment the expected internal link counter for all reads with an
  //    index > i.
  // 2. Increment expected counters for each specified separation for all
  //    right reads in B
  //
  // Then do the same for B -> A, except swap left and right definitions

  // The expected counts should be returned to R, so we will use R memory
  // allocation here:
  // 0: A internal expectation
  // 1: A -> B expect
  // 2: B internal expect
  // 3: B -> A expect
  // 4: A_i
  // 5: B_i
  // 6: A_obs information about linked reads
  // 7: B_obs
  // A_obs and B_obs are integer matrices with 3 columns:
  //   1. Flag: bit 1: read and mate within region
  //            bit 2: mate outside region
  //            bit 3: read is left or right (0,1 -> left, right)
  //   2. Pos (this may be modified to be region specific or reordered)
  //   3. Mate pos. If outside region, not included.
  SEXP ret_data = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(ret_data, 0, allocVector(REALSXP, A_pairs.n));
  SET_VECTOR_ELT(ret_data, 1, allocMatrix(REALSXP, seps_n, B_pairs.n));
  SET_VECTOR_ELT(ret_data, 2, allocVector(REALSXP, B_pairs.n));
  SET_VECTOR_ELT(ret_data, 3, allocMatrix(REALSXP, seps_n, A_pairs.n));
  SET_VECTOR_ELT(ret_data, 4, allocMatrix(INTSXP, sizeof(align_pair) / sizeof(int),  A_pairs.n));
  SET_VECTOR_ELT(ret_data, 5, allocMatrix(INTSXP, sizeof(align_pair) / sizeof(int), B_pairs.n));
  SET_VECTOR_ELT(ret_data, 6, allocMatrix(INTSXP, A_pairs.n, 3));
  SET_VECTOR_ELT(ret_data, 7, allocMatrix(INTSXP, B_pairs.n, 3));

  memcpy(INTEGER(VECTOR_ELT(ret_data, 4)), A_pairs.a, sizeof(align_pair) * A_pairs.n);
  memcpy(INTEGER(VECTOR_ELT(ret_data, 5)), B_pairs.a, sizeof(align_pair) * B_pairs.n);
  
  double *A_int = REAL(VECTOR_ELT(ret_data, 0));
  double *A_to_B = REAL(VECTOR_ELT(ret_data, 1));
  double *B_int = REAL(VECTOR_ELT(ret_data, 2));
  double *B_to_A = REAL(VECTOR_ELT(ret_data, 3));

  memset(A_int, 0, sizeof(double) * A_pairs.n);
  memset(A_to_B, 0, sizeof(double) * B_pairs.n * seps_n);
  memset(B_int, 0, sizeof(double) * B_pairs.n);
  memset(B_to_A, 0, sizeof(double) * A_pairs.n * seps_n);

  // I would be better to have a struct of some sort here. 
  int *A_obs_flag = INTEGER(VECTOR_ELT(ret_data, 6));
  int *B_obs_flag = INTEGER(VECTOR_ELT(ret_data, 7));
  int *A_obs_pos = A_obs_flag + A_pairs.n;
  int *A_obs_mpos = A_obs_pos + A_pairs.n;
  int *B_obs_pos = B_obs_flag + B_pairs.n;
  int *B_obs_mpos = B_obs_pos + B_pairs.n;

  memset(A_obs_flag, 0, sizeof(int) * A_pairs.n * 3);
  memset(B_obs_flag, 0, sizeof(int) * B_pairs.n * 3);

  Rprintf("A_pairs.n: %ld  B_pairs.n: %ld\n", A_pairs.n, B_pairs.n);
  //  size_t ii=0; size_t jj=0;
  for(size_t i=0; i < A_pairs.n; ++i){
    A_obs_pos[i] = A_pairs.a[i].pos; 
    A_obs_mpos[i] = A_pairs.a[i].m_pos;
    A_obs_flag[i] = A_pairs.a[i].m_pos < B_off | ((A_pairs.a[i].m_pos >= B_off) << 1) |
      ((A_pairs.a[i].pos > A_pairs.a[i].m_pos) << 2);
    // internal links first:
    // For internal links, we need to consider only ones that have a distance
    // of more than min_sep; this is because we only count such links for
    // the observation; the probability of such links are thus 0.
    for(size_t j=i+1; j < A_pairs.n; ++j){
      if( A_pairs.a[j].pos - A_pairs.a[i].pos >= min_sep)
	A_int[i] += 1 / (double)(A_pairs.a[j].pos - A_pairs.a[i].pos);
    }
    for(size_t j=0; j < B_pairs.n; ++j){
      for(size_t k=0; k < seps_n; ++k){ // we should possibly check for consecutive reads that are too close to each other
	if(seps[k] + B_pairs.a[j].pos - A_pairs.a[i].pos)
	  A_to_B[ j * seps_n + k ] += 1 / (double)( (B_pairs.a[j].pos - A_pairs.a[i].pos) + seps[k] );
      }
    }
  }
  // Then do the same for B,
  // Not sure that reverse order is needed given the refactoring to using A_pairs
  // and B_pairs. 
  // ssize_t necessary; or we could check with i < B_i.n-1
  for(ssize_t i=(B_pairs.n-1); i >=0; --i){
    B_obs_pos[i] = B_pairs.a[i].pos;
    B_obs_mpos[i] = B_pairs.a[i].m_pos;
    B_obs_flag[i] = B_pairs.a[i].m_pos >= B_off | ((B_pairs.a[i].m_pos < B_off) << 1) |
      ((B_pairs.a[i].pos > B_pairs.a[i].m_pos) << 2);
    for(ssize_t j=i-1; j >= 0; --j){
      if( B_pairs.a[i].pos - B_pairs.a[j].pos >= min_sep )
	B_int[i] += 1 / (double)(B_pairs.a[i].pos - B_pairs.a[j].pos);
    }
    for(size_t j=0; j < A_pairs.n; ++j){
      for(size_t k=0; k < seps_n; ++k){
	if(seps[k] + B_pairs.a[i].pos - A_pairs.a[j].pos > min_sep)
	  B_to_A[j * seps_n + k] += 1 / (double)( seps[k] + B_pairs.a[i].pos - A_pairs.a[j].pos );
      }
    }
  }
  kv_destroy(A_pairs);
  kv_destroy(B_pairs);
  UNPROTECT(3);
  return(ret_data);
}

static const R_CallMethodDef callMethods[] =
  {
   {"assess_region_links", (DL_FUNC)&assess_region_links, 2},
   {"assess_region_links_mt", (DL_FUNC)&assess_region_links_mt, 2},
   {"assess_separation", (DL_FUNC)&assess_separation, 4},
   {NULL, NULL, 0}
  };

void R_init_hiC_expt(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

