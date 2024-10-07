/********************************************************************
   dim_mgr.cpp - logic for the dimension manager
     for a n dimensional array where you need only m (< n) dimensions
     manage the storege for just the m-n dimensions for only the grid boxes
     needed - this saves reading a big block in when you only need some
     kind of slice.

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       8 Jan 2020      Original development
      W. Robinson       Nov 2020        make as a class in c++

********************************************************************/
#include "dim_mgr.hpp"

// for class dim_mgr
// Separate the constructor out

  // ndim is # of dimensions to be managed
/********************************************************************
   dim_mgr - 1st constructor form for dim_mgr class

   ndim   I  # of dimensions under management
********************************************************************/
  dim_mgr::dim_mgr(int32_t ndim) {
    hash_tbl_siz = 7;
    this->hash_tbl_siz = hash_tbl_siz;
    this->const_hlp(ndim);
    //cout << __FILE__ << ": 1-arg dim_mgr\n";
  }

/********************************************************************
   dim_mgr - 2nd constructor form for dim_mgr class

   ndim   I  # of dimensions under management
   hash_tbl_siz  I  # elements in the hash tablew
********************************************************************/
  dim_mgr::dim_mgr(int32_t ndim, int32_t hash_tbl_sz) {
    hash_tbl_siz = hash_tbl_sz;
    this->const_hlp(ndim);
    this->hash_tbl_siz = hash_tbl_siz;
    //cout << __FILE__ << ": 2-arg dim_mgr\n";
  }

/********************************************************************
   const_hlp  - helper for the 2 constructors
     This does the real set-up of the dimension manager

   ndim  I  # of dimensions under management
********************************************************************/
  void dim_mgr::const_hlp(int32_t ndim) {
    this->ndim = ndim;
    prev_int = NULL;
    n_tot_dat_blobs = 0;
    dim_info = ( dim_info_struc ** )malloc
      ( ndim * sizeof( dim_info_struc * ) );
    //cout << __FILE__ << ": Setting the # dims\n";
    // allocate all dim_info structs
    for( int i = 0; i < ndim; i++ )
      {
      dim_info[i] = (dim_info_struc *) malloc( sizeof(dim_info_struc) );
      dim_info[i]->dim_coords = NULL;  // so we know if not filled
      dim_info[i]->nvals = 0;  //  ditto
      dim_info[i]->type = 0;  //  ditto
      }
    // to set the hash table
    hash_tbl = new vector <hash_entry_struc>[hash_tbl_siz];
    gpt_hash_tbl = new vector <gpt_hash_struc>[hash_tbl_siz];

    // to set up the point information structure
    int ncorner = pow( 2, ndim );
    pt_info.pt_status = (int32_t *)malloc( ncorner * sizeof(int32_t ) );
    pt_info.pt_base_loc = (int32_t *)malloc( ndim * sizeof(int32_t ) );
    pt_info.wt = (double *)malloc( ndim * sizeof(double) );
    pt_info.wt_pt = (double *)malloc( ncorner * sizeof(double) );
    pt_info.dat_ptrs = (void **)malloc( ncorner * sizeof(void *) );
  }

/********************************************************************
   ~dim_mgr - destructor - needed for clean up of the class
********************************************************************/
  dim_mgr::~dim_mgr() {
    delete[] hash_tbl;
   delete[] gpt_hash_tbl;
    for( int i = 0; i < ndim; i++ )
    {
      if( dim_info[i]->type == 1 )
        free( dim_info[i]->dim_coords );
      free(dim_info[i]);
    }
    //  dispense with point info
    free( pt_info.pt_status );
    free( pt_info.pt_base_loc );
    free( pt_info.wt );
    free( pt_info.wt_pt );
    free( pt_info.dat_ptrs );
    // 
    free( dim_info );
  }

  // record an addition of points
/********************************************************************
   add_pts - add to the total count of data blobs under management

    Returns:
     void - none
   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t           nptadd           I      # grid points added
********************************************************************/
  void dim_mgr::add_pts( int32_t nptadd ) {
    n_tot_dat_blobs += nptadd;
  }

/********************************************************************
   get_tot_blobs - get total # data blobs being managed

   Returns:
     int32_t # data blobs being managed
********************************************************************/
   int32_t dim_mgr::get_tot_blobs() { return n_tot_dat_blobs; }

/********************************************************************
   init_dim - initialize information for a dimension.  First form is for 
     the dimension coordinates specified in an array

    Returns:
     int32_t - status 0 = good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t           dim_num          I      Dimension to set
      int32_t           nvals            I      # grid pts in this dim
      double            min              I      min of dimension range
      double            max              I      max of dimension range
********************************************************************/
  int32_t dim_mgr::init_dim( int32_t dim_num, int32_t nvals, 
    double min, double max ) {
    // set the dim_info_struc if not set yet
    if( dim_info == NULL ) {
      cout << __FILE__ << __LINE__ << " dim_info array is null\n";
      return 0;
    }
    if( dim_info[dim_num]->nvals != 0 ) {
      cout << __FILE__ << __LINE__ << " dim_info array entry " << dim_num <<
       " is filled already\n";
      return 0; 
    }
    dim_info[dim_num]->nvals = nvals;
    dim_info[dim_num]->min = min;
    dim_info[dim_num]->max = max;
    dim_info[dim_num]->type = 0;
    dim_info[dim_num]->delta = ( max - min ) / ( nvals - 1 );
    return 0;
  }

/********************************************************************
   init_dim - initialize information for a dimension.  Second form is for 
      the dimension coordinates specified in an array

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t           dim_num          I      Dimension to set
      int32_t           nvals            I      # grid pts in this dim
      double *          dim_coords       I      set of coordinates for this 
                                                  dimension
********************************************************************/
  int32_t dim_mgr::init_dim( int32_t dim_num, int32_t nvals, 
    double *dim_coords ) {
  
    if( dim_info == NULL ) {
      cout << __FILE__ << __LINE__ << " dim_info array is null\n";
      return 0;
    }
    if( dim_info[dim_num]->dim_coords != NULL ) {
      cout << __FILE__ << __LINE__ << " dim_info array entry " << dim_num <<
      " is filled already\n";
      return 0;
    }
    dim_info[dim_num]->nvals = nvals;
    dim_info[dim_num]->min = dim_coords[0];
    dim_info[dim_num]->max = dim_coords[ nvals - 1 ];
    dim_info[dim_num]->type = 1;
    dim_info[dim_num]->dim_coords = 
      (double *)malloc( nvals * sizeof( double ) );
    for( int i = 0; i < nvals; i++ )
      dim_info[dim_num]->dim_coords[i] = dim_coords[i];
    return 0;
  }

/********************************************************************
   update_new_data - Update info about new data under management

    This will refresh the interval's data status and data pointer using
    the pt_info structure the user updated.

********************************************************************/
  void dim_mgr::update_new_data( ) {
    int32_t npt = pow( 2, ndim );

    //  just may as well update all points
    //  WDR 27 Aug 2024 do not update if prev_int is NULL
    // to test use of prev_int
    if( prev_int != (interval_struc *)NULL )
      {
      for( int i = 0; i < npt; i++ )
        {
        prev_int->dat_info[i]->dat_status = pt_info.pt_status[i];
        prev_int->dat_info[i]->dat = pt_info.dat_ptrs[i];
        }
      }
  }

/********************************************************************
   mng_pt - for a point in a grid, find or set up the location for the 
     associated data and return that and the weights to apply

   Returns:
     pt_info_struc *   the structure with all point information
                       needed to interpolagte the data

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      double *          pt               I      coordinates of the point to
                                                locate
      int32_t           access_id        I      A last time of access to
                                                assign to the bottom interval
                                                the lower the value, the
                                                more time since last access
      int32_t *         status           O      return status 0 = good, 
                                                  1 = code failure, you 
                                                  shoiuld fail at this point
                                                  2 = not valid point (not in
                                                  the grid bounds) treat as 
                                                  unable to do the point
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       1 Sep 2020      based loosly on sparse_mng_pt()
********************************************************************/
  pt_info_struc * dim_mgr::mng_pt( double *pt, int32_t access_id, 
    int32_t *status ) {

    int32_t idim, idat, ndat, same_interval;
    int32_t *ix, *ix_real, *offset_lcl, acc_mode, ret;
    interval_struc *found_int;

    *status = 0;  // start with good status
   /*
    *  If the new point is the same as the old, just return the same point info
    */
    //  not in std c++ so... if( isnormal( old_pt[0] )
    if( old_pt[0] == old_pt[0] )
      {
      int32_t match = 1;
      for( idim = 0; idim < ndim; idim++ )
        {
        if( *( pt + idim ) != *( old_pt + idim ) )
          {
          match = 0;
          break;
          }
        }
      if( match == 1 )
        {
        *status = old_status;
        return &pt_info;
        }
      }
    for( idim = 0; idim < ndim; idim++ )
      *( old_pt + idim ) = *( pt + idim );
   /*
    *  set up some local storage
    */
    if( ( ( ix = (int32_t *) malloc( ndim * sizeof( int32_t ) ) ) == NULL ) ||
        ( ( ix_real = (int32_t *)
            malloc( ndim * sizeof( int32_t ) ) ) == NULL ) ||
        ( ( offset_lcl = (int32_t *)
            malloc( ndim * sizeof( int32_t ) ) ) == NULL ) )
      {
      printf( "%s, %d: Unable to allocate the dim index array\n", __FILE__,
        __LINE__ );
      *status = 1;
      return &pt_info;
      }
   /*
    * also set the weight and point weight if needed
    */
    ndat = pow( 2., ndim );
   /*
    *  get the index of the point in all dims and the weight to it
    */
    if( sparse_get_loc( dim_info, ndim, pt, ix, pt_info.wt, pt_info.wt_pt )
      != 0 )
      {
      *status = 2;
      old_status = 2;
      return &pt_info;
      }
    /*
    printf( "\n\n%s - index: %d,%d, weight: %f,%f\n", __FILE__,
        ix[0], ix[1], pt_info.wt[0], pt_info.wt[1] );
    printf( "data point: (%f, %f)\n", pt[0], pt[1] );
    */
   /*
    *  Its likely that the previous point is in the same grid box as the current
    *  point.  In this case, just return the previously found interval
    */
    same_interval = 0;
    pt_info.interval_needs_data = 0;
    if( prev_int != NULL )
      {
     /* see if we are in the same interval */
      same_interval = 1;
      for( idim = 0; idim < ndim; idim++ )
        {
        if( ix[idim] != prev_int->dat_info[0]->ix_arr[idim] )
          {
          same_interval = 0;
          break;
          }
        }
      }
  
    if( same_interval == 1 )
      {
      found_int = prev_int;
      found_int->access_id = access_id;
      }
    else
      {
     /*
      *  get/set the interval for this grid point
      */
      acc_mode = 0;
      valarray<int> ix_val(ndim);
      for( int i = 0; i < ndim; i++ )
        ix_val[i] = ix[i];
      found_int = access( ix_val, acc_mode );
  
     /*
      *  Set the found interval's access id to the one passed in
      */
      found_int->access_id = access_id;
     /*
      *  if interval is new, we need to donate any data areas shared
      *  to the new interval
      */
      if( found_int->dat_filled == 0 )
        {
        idim = -1;  // set this way to start at the top
  
        if( share_gpt( found_int, ix ) != 0 ) {
          *status = 1;
          return &pt_info;
        }
        //  this should be obsolete
       /*  set the interval filled value to yes */
        found_int->dat_filled = 1;
       /*  loop through all data pointers in the new interval */
       /*  and set up the data structures, less the pointer to the data */
        for( idat = 0; idat < ndat; idat++ )
          {
         /*  if the data for this corner is not there, create and fill it */
          if( found_int->dat_info[idat] == NULL )
            {
            valarray<int> ixv_off(ndim);
            pt_info.interval_needs_data = 1;  // flag that data needed 
                                               // for this interval / grid box
           /* set up each data info struct */
            found_int->dat_info[idat] =
              (data_info_struc *) malloc( sizeof( data_info_struc ) );
           /*  get index for this point */
            linear_to_offset( ndim, idat, offset_lcl );
            found_int->dat_info[idat]->ix_arr = (int32_t *) malloc( ndim *
              sizeof( int32_t ) );
            for( idim = 0; idim < ndim; idim++ )
              {
              ix_real[idim] = ix[idim] + offset_lcl[idim];
              ixv_off[idim] = ix_real[idim];
              found_int->dat_info[idat]->ix_arr[idim] = ix_real[idim];
              }
           /*
            *  set data status to not filled and # intervals pointing to it to 1
            */
            found_int->dat_info[idat]->dat_status = -1; /* not filled here */
            found_int->dat_info[idat]->n_intervals = 1;
           /*
            *  add new data_info to the grid point hash table
            */
            ret = gpt_add( ixv_off, found_int->dat_info[idat] );
            if( ret != 0 ) 
              {
              *status = 1;
              return &pt_info;
              }
            }
          else /* for data already there, incriment access count */
            {
            // WDR keep in case some functions needed here
            }
          }
        }
     /* remember the found interval and its dim location */
      prev_int = found_int;
      }
  
  /* end set up of (possible) new grid interval */
  
    free( ix ); free( ix_real ); free( offset_lcl );
  
   //  put the needed info into the pt_info struct
    for( int i = 0; i < ndim; i++ )
      pt_info.pt_base_loc[i] = found_int->dat_info[0]->ix_arr[i];
    for( int i = 0; i < ndat; i++ ) {
      pt_info.pt_status[i] = found_int->dat_info[i]->dat_status;
      // pt_info.wt_pt, and wt are already set
      pt_info.dat_ptrs[i] = found_int->dat_info[i]->dat;
    }
   /*
    *  return success
    */
    *status = 0;
    old_status = 0;
    return &pt_info;
  }
  // all the rest of routines that go with mng_pt

  int32_t dim_mgr::sparse_get_loc( dim_info_struc **dim_info, int32_t ndim, 
     double *pt, int32_t *ix, double *wt, double *wt_pt )
/********************************************************************
   sparse_get_loc - get the location of a point in all grid dimensions
     and the weights

   Returns:
     int32_t - status 0 = good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      dim_info_struc ** dim_info         I      array of descriptions of each
                                                dimension
                                                that is handled
      int32_t           ndim             I      # dimensions handled
      double *          pt               I      coordinates of the point to
                                                locate
      int32_t *         ix               O      index of the grid interval
                                                for all dimensions
      double *          wt               O      for each dimension, the weight
                                                to apply to the index point
                                                (1 - wt for the next index
                                                point)
      double *          wt_pt            O      set of consolidated weights to
                                                apply to the data blobs pointed
                                                to by the dat_info_struc
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       26 Nov 2019     Original development

********************************************************************/
  {
  int32_t idim, iint, full_len, msk;
  /*
   *  we start with all composite weights at 1
   */
  full_len = pow( 2, ndim );
  for( idim = 0; idim < full_len; idim++ )
    *( wt_pt + idim ) = 1.;
 /*
  loop through all dimensions and get the point and weight
  */
  for( idim = 0; idim < ndim; idim++ )
    {
   /* for now, assume the grid is uneven and we look through each interval */
    if( ( pt[idim] < dim_info[idim]->min ) ||
        ( pt[idim] > dim_info[idim]->max ) )
      {
      //printf( "%s, %d - point is outside grid dimension range\n",
      //  __FILE__, __LINE__ ); 
      return 1;
      }
    for( iint = 1; iint < dim_info[idim]->nvals; iint++ )
      {
      if( pt[idim] <= dim_info[idim]->dim_coords[iint] )
        { 
        ix[idim] = iint - 1;

        wt[idim] = ( pt[idim] - dim_info[idim]->dim_coords[iint-1] ) /
          ( dim_info[idim]->dim_coords[iint] -
            dim_info[idim]->dim_coords[iint-1] );

        break;
        }
      }
   /*
    *  This will modify the composite weight so it applies to each
    *  corner's data - apply the weight to the end point
    *  of that grid else 1 - weight
    */ 
    msk = pow( 2, idim );
    for( iint = 0; iint < full_len; iint++ )
      {
      if( ( iint & msk ) == 0 )  /* apply 1 - wt */
        *( wt_pt + iint ) *= 1. - wt[idim]; 
      else 
        *( wt_pt + iint ) *= wt[idim];
      }
    }
  return 0;
  }

interval_struc * dim_mgr::access( valarray<int>& s, int32_t acc_mode )
/********************************************************************
   access  search / insert interval defining the data at the grid point

   Returns:
     interval_struc *  pointer to the interval information

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      valarray<int>     s                I      set of coordinates of the point
      int32_t           acc_mode         I      action: 0 - find the interval,
                                                if not there, make it.  return
                                                the interval pointer, 1 -
                                                find interval and return NULL
                                                if not found, 2 - remove an
                                                entry based on coordinates s
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       31 Aug 2020     Original development
********************************************************************/
  {
  int32_t ndat = pow( 2, ndim );

  //Compute the index by using the hash function
  int index = hash_func( s, hash_tbl_siz );
  //Search the linked list at that specific index for the index array
  for( long unsigned int i = 0; i <  hash_tbl[index].size(); i++)
    {
    //  the == for the 2 valarrays make a valarray of booleans, 
    //  basically all dims match (all 1), so min is not = 0
    //  if all dimensions match...
    if( (  hash_tbl[index][i].ix_arr == s ).min() != 0 )
      {
      // interval is found, we can either return interval info or erase it
      if( acc_mode == 2 )  // remove NOTE that underlying structures
                           // not handled here!!! so fix before use
        {
        hash_tbl[index].erase( hash_tbl[index].begin() + i );
        cout << "Set is removed" << endl;
        return (interval_struc *) NULL;
        }
      else      // return interval struc for the interval
        {
        //cout << "Set is found!" << endl;
        //cout << "position: " << i << endl;
        return hash_tbl[index][i].int_str;
        }
      }
    }
  //  if not found,
  //cout << "Set is not found!" << endl;
  if( acc_mode == 0 )  // add a new entry
    {
    //  create the new entry
    hash_entry_struc hash_entry;
    interval_struc *insert;

    insert = (interval_struc *)malloc( sizeof(interval_struc) );
    insert->dat_filled = 0;
    insert->access_id = 0; //  tenp setting, done later
    insert->dat_info = (data_info_struc **)
      calloc( ndat, sizeof(data_info_struc *) );
    hash_entry.ix_arr = s;
    hash_entry.int_str = insert;
    // Insert the element in the linked list at the particular hash index

    hash_tbl[index].push_back(hash_entry);
    //cout << "Inserting new interval, hash item\n";
    return insert;
    }
  else
    {
    return (interval_struc *) NULL;
    }
  cout << "Code should not get here!/n";
  return (interval_struc *)NULL;
  }

//  New - add a data info structure to the grid point hash table

int32_t dim_mgr::gpt_add( valarray<int> s, data_info_struc *dat_info )
/*-------------------------------------------------------------------
   gpt_add - add a data info structure to the grid point hash table

   Returns:  int 0 all OK, 1 trying to insert to a existing coordinate

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      valarray<int>     s                I      set of coordinates of the point
      data_info_struc * dat_info         I      descriptor for the data at the 
                                                grid point

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       30 Nov 2021     replace interval searching with this

-------------------------------------------------------------------*/
  {
  // as with access, Compute the index by using the hash function
  // and check for this coordinate already (there should not be one)
  int index = hash_func( s, hash_tbl_siz );

  for( long unsigned int i = 0; i < gpt_hash_tbl[index].size(); i++)
    {
    if( ( gpt_hash_tbl[index][i].ix_arr == s ).min() != 0 )
      {
      // we have a match, which should not happen as this is insert
      cout << __FILE__ << __LINE__ << 
        "Error: grid point is already in hash table, Exiting" << endl;
      return 1;
      }
    }
  //  Normal branch, insert the new hash entry
  gpt_hash_struc hash_entry;

  hash_entry.ix_arr = s;
  hash_entry.dat_info = dat_info;

  gpt_hash_tbl[index].push_back( hash_entry );
  /*
   cout << __FILE__ << ", " << __LINE__ << ", " << 
   "Just added a grid point to hash" << endl;
  */
  return 0;
  }

int dim_mgr::hash_func(valarray<int> s, int32_t hash_tbl_siz )
/********************************************************************
   hash_func - from the grid coordinates, make a hash entry number
   This is one of the more odd parts of using a hashtable - get a good
   distribution of table entries or the hash is not as useful

   Returns:
     interval_struc *  pointer to the interval information

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      valarray<int>     s                I      set of coordinates of the point
      int32_t           hash_tbl_siz     I      size of hash table

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       31 Aug 2020     Original development
********************************************************************/
  {
  int retval;
  long sum = 0;

  if( hash_ifirst == 0 ) {
    int32_t ndim = s.size();
    hash_ifirst = 1;
    hash_h_mult = hash_tbl_siz / ( 2 * ndim );
    hash_h_mult = hash_h_mult - ndim / 2;
    }

  for( long unsigned int i = 0; i < s.size(); i++ )
    sum += ( hash_h_mult + i ) * s[i];
  sum = sum  % hash_tbl_siz;

  retval = sum;
  return retval;
  }

int32_t dim_mgr::share_gpt( interval_struc *intvl, int32_t *ix )
/********************************************************************
   share_gpt  will find grid points that are already set up and share 
     them with the current interval

   Returns:
     int32_t - 0 if all good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      interval_struc *  intvl           I/O     a new interval structure to
                                                find data for
      int32_t *         ix               I      start grid location of the 
                                                interval
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       1 Dec 2021      to hopefully simplify the grid point 
                                        data sharing process, replacing the 
                                        iterative int_share_data
********************************************************************/
  {
  int32_t *ix_off, i;
  valarray<int> s(ndim);
  int32_t npt = pow( 2, ndim );
  ix_off = ( int32_t * ) malloc( ndim * sizeof(int32_t) );
  //
  // look at all the corners of this interval
  for( int ipt = 0; ipt < npt; ipt++ )
    {
    // get the grid point location from base and offset to this point
    linear_to_offset( ndim, ipt, ix_off );
    for( i = 0; i < ndim; i++ )
      {
      ix_off[i] += ix[i];
      s[i] = ix_off[i];
      }

    // see if the point is in the grid point hash already
    int index = hash_func( s, hash_tbl_siz );
    for( i = 0; i < (int32_t)gpt_hash_tbl[index].size(); i++ )
      {
      if( ( gpt_hash_tbl[index][i].ix_arr == s ).min() != 0 )
        {
        // found matching pt it in hash table, point to it from the interval
        intvl->dat_filled = 1;
        intvl->dat_info[ipt] = gpt_hash_tbl[index][i].dat_info;
        intvl->dat_info[ipt]->n_intervals++;  // one more ref to the point
        break;
        }
      }
    //  new points are handled in mng_pt
    }
  free( ix_off );
  return 0;
  } 

/********************************************************************
   purge  will completely remove all entries from the dimension manager 

   Returns:
     int32_t - 0 all good

********************************************************************/
int32_t dim_mgr::purge()
  {
  int32_t nhash = hash_tbl_siz;
  int32_t ndat = pow( 2, ndim );
  //  go to hash entry
  //  loop over each hash bin
  for( int ihash = 0; ihash < nhash; ihash++ )
    {
    int32_t nentry = hash_tbl[ihash].size();
    for( int iint = 0; iint < nentry; iint++ )
      {
      // go to interval struc
      hash_entry_struc hash_ent = hash_tbl[ihash][0];

      interval_struc *int_str = hash_ent.int_str;

      // go to data structures
      for( int idat = 0; idat < ndat; idat++ )
        {
        // rid data if # intervals accessing is 1
        if( int_str->dat_info[idat]->n_intervals > 1 )
          {
          int_str->dat_info[idat]->n_intervals--;
          }
        else
          {
          // no more intervals point to this data struct,
          // so free the data info structure
          free( int_str->dat_info[idat]->dat );
          free( int_str->dat_info[idat]->ix_arr );
          free( int_str->dat_info[idat] );
          }
        }

        // free the array with data info structure addresses
        free( int_str->dat_info );

        // free the interval struc
        free( int_str );
        // rid the top hash entry for bin i
        hash_tbl[ihash].erase(hash_tbl[ihash].begin() );
      }
    }

  //  Clean up all the grid point hash table entries
  //  the dat_info is already gone from the interval purge above, so
  //  the only thing left is all the entries in each gpt_hash
  for( int ihash = 0; ihash < nhash; ihash++ )
    {
    int32_t nentry = gpt_hash_tbl[ihash].size();
    for( int iint = 0; iint < nentry; iint++ )
      {
      gpt_hash_tbl[ihash].erase( gpt_hash_tbl[ihash].begin() );
      }
    }

  //  the previous interval pointer is invalid, say so
  prev_int = NULL;
  n_tot_dat_blobs = 0;
  //  The rest can stay - the destructor handles that (I think)
  return 0;
  }

/********************************************************************
  prune - Remove selected intervals from management

  Returns:  int32_t - 0 all good
  access_id  I   erase intervals with access_id lower than access_id
********************************************************************/
int32_t dim_mgr::prune(int32_t access_id ) {

  // start at dim_mgr
  int32_t nhash = hash_tbl_siz;
  int32_t ndat = pow( 2, ndim );

  //  go to hash entry,  loop over each hash bin
  for( int ihash = 0; ihash < nhash; ihash++ )
    {
    // prune from list end to start so next entry location is unchanged
    // even if an entry is deleted
    int32_t nentry = hash_tbl[ihash].size();
    for( int iint = ( nentry - 1 ); iint >= 0; iint-- )
      {
      // go to interval struc
      hash_entry_struc hash_ent = hash_tbl[ihash][iint];

      interval_struc *int_str = hash_ent.int_str;
      if( int_str->access_id < access_id )
        {
        // go to data structures
        for( int idat = 0; idat < ndat; idat++ )
          {
          // rid data if # intervals accessed is 1 else lower interval count
          if( int_str->dat_info[idat]->n_intervals > 0 )
            {
            int_str->dat_info[idat]->n_intervals--;
            }
          else
            {
            // no more intervals point to this data struct,
            // leave the data info struc around to be handled when
            // cleaning up the grid point hash
            }
          }
        // WDR *** MAY need to free( int_str->dat_info );
        free( int_str->dat_info );
        // free the interval struc
        free( int_str );
        // rid the hash entry for hash ihash, list entry iint
        hash_tbl[ihash].erase(hash_tbl[ihash].begin() + iint);
        }
      }
    }
  //  go through the grid point hash table and remove entries with 
  //  no accesses
    for( int ihash = 0; ihash < nhash; ihash++ )
    {
    // prune from list end to start as before
    int32_t nentry = gpt_hash_tbl[ihash].size();
    for( int iint = ( nentry - 1 ); iint >= 0; iint-- )
      {
      // go to data info struc
      gpt_hash_struc hash_ent = gpt_hash_tbl[ihash][iint];
      data_info_struc *lcl_info = hash_ent.dat_info;
      if( lcl_info->n_intervals == 0 )
        {
        //  as no intervals point to this data, remove it
        free( lcl_info->dat );
        free( lcl_info->ix_arr );
        free( lcl_info );
        n_tot_dat_blobs--;  // one less managed data element
        gpt_hash_tbl[ihash].erase( gpt_hash_tbl[ihash].begin() + iint );
        }
      }
    }
  //  This is much like purge, but should work a bit slower
  prev_int = NULL;  //  resetprevious interval in case we removed it
  return 0;
  }

/********************************************************************
   dump_mgr - report all info on all managed intervals

   double  pt  I  depth of reporting: 0  do all, 1 - leave out interval info
********************************************************************/
void dim_mgr::dump_mgr( double *pt ) {
  // print # bins in hash table
  int32_t idepth = (int32_t)pt[0];
  printf( "\nCurrent manager state:\n" );
  printf( "  # of hash tbl bins: %d,  # data blobs managed: %d\n", 
    hash_tbl_siz, n_tot_dat_blobs );

  // loop all hash bins
  for( int ibin = 0; ibin < hash_tbl_siz; ibin++ ) {
    // print bin # an entries inbin
    int32_t n_entry = hash_tbl[ibin].size();
    printf( "  bin # %d, # entries: %d\n", ibin, n_entry );
    // loop bin entries
    for( int ient = 0; ient < n_entry; ient++ ) {
      // print entry # and interval index associated with this entry
      printf( "    Entry %d, interval index for this entry:\n    ",
        ient );
      for( int ix = 0; ix < ndim; ix++ )
        printf( " %d,", hash_tbl[ibin][ient].ix_arr[ix] );
      printf( "\n" );
      if( idepth < 1 ) {
        printf( "-------------------------\n" );
        printf( "Interval Information:\n" );
        //  call the interval dumper
        dump_interval( hash_tbl[ibin][ient].int_str );
        printf( "-------------------------\n" );
      }
    }
  }
  // Add a report for the grid point hash bins
  printf( "\n\nGrid point hash table summary:\n" );
  for( int ibin = 0; ibin < hash_tbl_siz; ibin++ )
    {
    // print bin # an entries inbin
    int32_t n_entry = gpt_hash_tbl[ibin].size();
    printf( "-------------------------\n" );
    printf( "  bin # %d, # entries: %d\n", ibin, n_entry );
    // loop bin entries
    for( int ient = 0; ient < n_entry; ient++ ) {
      // print entry #
      printf( "    Entry #: %d,  Point coords follow\n", ient );
      for( int ix = 0; ix < ndim; ix++ )
        printf( " %d,", gpt_hash_tbl[ibin][ient].ix_arr[ix] );
      printf( "\n" );
      printf( "data info Point coords\n" );
      for( int ix = 0; ix < ndim; ix++ )
        printf( " %d,", gpt_hash_tbl[ibin][ient].dat_info->ix_arr[ix] );
      printf( "\n" );
      printf( "data info # intervals pointed to: %d\n", 
        gpt_hash_tbl[ibin][ient].dat_info->n_intervals );
      printf( "data info status: %d\n",
        gpt_hash_tbl[ibin][ient].dat_info->dat_status );
      printf( "data info address: %ld\n",
       (long int)gpt_hash_tbl[ibin][ient].dat_info->dat );
      }
    }
  }

/********************************************************************
   dump_last_interval - report information about last interval 
     = a grid box's info

  interval_struc *  interval         I      Interval structure for the

********************************************************************/
void dim_mgr::dump_last_interval() { dump_interval( prev_int ); }

/********************************************************************
   dump_interval - report information about an interval = a grid box's info

   interval_struc *  interval         I      Interval structure for the

********************************************************************/
void dim_mgr::dump_interval( interval_struc *interval ) {
    if( interval != NULL ) {
    int32_t npt = pow( 2, ndim );
    printf( "Interval address: %ld  access ID: %d\n",(long int)interval,
      interval->access_id );
    printf( "    Grid loc (%d dims): ", ndim );
    for( int i = 0; i < ndim; i++ )
      printf( " %d", interval->dat_info[0]->ix_arr[i] );
    printf( "\n" );
    for( int i = 0; i < npt; i++ )
      {
      printf( 
        "# %d data blob, address: %ld,  # intervals pointing to it: %d\n", i,
        (long int)interval->dat_info[i]->dat, 
        interval->dat_info[i]->n_intervals );
      printf( "    status: %d, grid point indexes: \n    ",
        interval->dat_info[i]->dat_status );
      for( int j = 0; j < ndim; j++ )
        printf( "%d,  ", interval->dat_info[i]->ix_arr[j] );
      printf( "\n" );
      }
    } else {
    printf( "Sorry, last interval is NULL (no intervals added yet, " );
    printf( "or purge or prune done)\n" );
    }
  }

  // information reporting
  int dim_mgr::get_ndim() { return ndim; }
  int dim_mgr::get_hdim() { return hash_tbl_siz; }

// The linear_to_offset does not need to be in the class, in fact, it 
// himders use from outside.  So make it a non-member function here

/********************************************************************
   linear_to_offset will convert a linear location in the data pointer
      array to a offset array

   Returns: int32_t - 0 if all good

      int32_t           n_dim            I      # dimensions to use
      int32_t           dat_ix           I      linear offset value
      int32_t *         offset           O      offset array

********************************************************************/
int32_t linear_to_offset( int32_t n_dim, int32_t dat_ix, int32_t *offset )
  {
  int32_t idim;

  for( idim = 0; idim < n_dim; idim++ )
    {
   /*  the offset for that dim is the modulo 2 of the linear,
       then next time divide te linear offset by 2 to do the next dim offset */
    offset[idim] = dat_ix % 2;
    dat_ix /= 2;
    }
  return 0;
  }
