#include <iostream>
#include <cstring>
using namespace std;
#include <vector>
#include <valarray>
#include <math.h>

/* set the structure set - just use the non-OO versions - maybe they'll change
*/
/* define characteristics for each dimension */
struct dim_info_struc
  {
  int32_t nvals;  /* # values in this grid dim */
  double min;  /* min and max of this dim */
  double max;
  double delta;  // if evenly distributed, type 0, point separation
  int32_t type; /* type of vals along dim: 0 - evenly distrib, 1 - custom */
  double *dim_coords;  /* if type = 1, an array of the grid point values */
  };

/* data info structure - one for each data group of any # dimensions */
struct data_info_struc
  {
  int32_t n_intervals;  /* # of bottom intervals pointing to this data */
  int32_t dat_status;  /* status of this data point: -1 - not set yet,
                           0 - no data available, 1 - data is there */
  int32_t *ix_arr;  /* size n_dim array of indicies for this point */
  void *dat;   /* points to the data storage area */
  };

struct interval_struc
  {
  int32_t access_id;  /* A 'last time of reference' index for this interval
                         if it is a bottom interval - identifies old
                         intervals that haven't been used before a certain
                         value */
  int32_t dat_filled;  /*  status of the data pointed to: 0 - not filled,
                                                          1 - filled */
  data_info_struc  **dat_info;
  };

// the hash_entry_struc is the hash table item for a managed intervals 
struct hash_entry_struc
  {
  valarray<int> ix_arr;
  interval_struc *int_str;
  };

//  Similarly, the grid point hash table item is a convenient way to get 
//  to the individual grid points
struct gpt_hash_struc
  {
  valarray<int> ix_arr;
  data_info_struc *dat_info;
  };

  // The pt_info_struc is a convenient collector for the info about a 
  //   found point, and it keeps all that together.  We'll have the 
  //   caller set up the storage for it
struct pt_info_struc
  {
  int32_t interval_needs_data;  // flag saying the interval needs data
  int32_t *pt_status;   //  point status for the 2^ndim corner points: 
                        //    -1 - not set yet, 0 - no data 
                        //    available, 1 - data is there
  int32_t *pt_base_loc; // size ndim base location of the found point in 
                        //   the grid
  double *wt;           // weight for each dimension
  double *wt_pt;        // size 2^ndim weight for that point
  void **dat_ptrs;    // size 2^ndim pointers to the data blobs for that 
                        //    corner
  };

/* Looks like declaring the classes here is best, with defining in the 
   .cpp of the same name - I don't like separating the implementation from 
   the definition, but it seems the only way if you have multi-file 
   programs
*/
class dim_mgr {
  public:
  // 2 different constructors
  dim_mgr(int32_t);
  dim_mgr(int32_t, int32_t);
  // just 1 destructor
  ~dim_mgr();
  // for init of dims
  int32_t init_dim( int32_t, int32_t, double * );
  int32_t init_dim( int32_t, int32_t, double, double );
  // for setting the data blob address
  void set_store( void *, int32_t );
  void update_new_data( );
  // for the getting of the point and recording it in the table
  pt_info_struc *  mng_pt( double *, int32_t, int32_t * );
  // to clear all point entries or selected ones
  int32_t purge();
  int32_t prune(int32_t);
  
  // info getting
  void dump_last_interval();
  void dump_mgr( double * );

  // manage total # data blobs
  int32_t get_ndim();
  int32_t get_hdim();
  void add_pts( int32_t );
  int32_t get_tot_blobs();

  private:
  // helper for the constructors
  void const_hlp(int32_t);
  // other routines
  // help for mng_pt
  int32_t sparse_get_loc( dim_info_struc **, int32_t, double *,
    int32_t *, double *, double * );
  interval_struc *access( valarray<int>&, int32_t );
  int32_t gpt_add( valarray<int>, data_info_struc * );
  int hash_func(valarray<int>, int32_t );
  int32_t share_gpt( interval_struc *, int32_t * );
  void dump_interval( interval_struc * );

  // remember prior state of point sent in to look for and search status
  double old_pt[10] = { NAN };
  int32_t old_status = 2;
  //  for hash table
  int32_t hash_tbl_siz;
  int32_t hash_h_mult, hash_ifirst = 0;

  int32_t ndim, n_tot_dat_blobs;
  pt_info_struc pt_info;
  interval_struc *prev_int;  // remember the previous interval here
  dim_info_struc **dim_info;
  vector <hash_entry_struc> *hash_tbl;  // array of vectors for each hash
                                      // table entry
  vector <gpt_hash_struc> *gpt_hash_tbl;  // Similar to above but for the 
                                          // grid points

  };
//  The non-member function declarations
  int32_t linear_to_offset( int32_t, int32_t, int32_t * );
