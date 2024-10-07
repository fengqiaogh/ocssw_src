#include "cgal_interp.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triangulation;
typedef CGAL::Interpolation_traits_2<K>    Traits;
typedef K::FT                              Coord_type;
typedef K::Point_2                         Point;
typedef std::map<Point, Coord_type, K::Less_xy_2>         Coord_map;
typedef CGAL::Data_access<Coord_map>                      Value_access;

typedef std::vector<std::pair<Point, Coord_type> > Point_coordinate_vector;

static const int pool_size = 10;
static bool cpool[pool_size]={false,false,false,false,false,false,false,false,false,false};
static std::vector< std::vector<std::pair<Point, Coord_type> >> cov[pool_size];
static Coord_type *nrms[pool_size];

int cgal_nnc( int nxy, float *x, float *y, int nxyq, float *xq, float *yq) 
{
  Delaunay_triangulation T;

  int nnc_tag = -1;
  for (int i=0 ; i<pool_size; i++) {
    if (cpool[i] == false) {
      cpool[i] = true;
      nnc_tag = i;
      //std::cout << "Assigning NNC tag: " << i << std::endl;
      break;
    }
  }
  if(nnc_tag == -1) {
      std::cout << "Error: no NNC tags left\n";
      return -1;
  }

  for (int i=0 ; i<nxy ; i++){
    K::Point_2 p(x[i],y[i]);
    T.insert(p);
  }
  
  nrms[nnc_tag] = new Coord_type[nxyq];
  
  // Generate interpolated values at xy gridpoints
  for (int i=0 ; i<nxyq ; i++) {
    //    if (i % 100000 == 0) std::cout << i << std::endl;

    K::Point_2 p(xq[i], yq[i]);
    std::vector<std::pair<Point, Coord_type> > coords;

    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>,
                 K::FT, bool> nn_result =
      CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords));

    bool status = nn_result.third;
    if (status) 
      nrms[nnc_tag][i] = nn_result.second; 
    else 
      nrms[nnc_tag][i] = -1;
    cov[nnc_tag].push_back(coords);
  }

  return nnc_tag;
}

int cgal_interp2( int nxy, float *x, float *y, float *v, int nxyq, float *vq,
                     int nnc_tag)
{
  Coord_map value_function;

  if(nnc_tag < 0 || nnc_tag >= pool_size) {
      std::cout << "Error: bad NNC tag: " << nnc_tag << std::endl;
      return -1;
  }

  for (int i=0 ; i<nxy ; i++){
    K::Point_2 p(x[i],y[i]);
    value_function.insert(std::make_pair(p, v[i]));
  }

  // Generate interpolated values at xy gridpoints
  for (int i=0 ; i<nxyq ; i++) {
    //    K::Point_2 p(xq[i], yq[i]);
    std::vector<std::pair<Point, Coord_type> > coords;

    Coord_type norm ;
    coords = cov[nnc_tag][i];
    norm = nrms[nnc_tag][i];

    if (coords.size() > 0) {
      Coord_type res =
        CGAL::linear_interpolation(coords.begin(), coords.end(), norm,
                                   Value_access(value_function));
      vq[i] = res;
    } else {
      vq[i] = -999;
    }
   }

  return EXIT_SUCCESS;
} 

int cgal_release_tag( int nnc_tag) {

  if(nnc_tag < 0 || nnc_tag >= pool_size) {
      std::cout << "Error: bad NNC tag: " << nnc_tag << std::endl;
      return -1;
  }

  if (cpool[nnc_tag] == true) {
    //std::cout << "Releasing NNC tag: " << nnc_tag << std::endl;
    cpool[nnc_tag] = false;
    delete[] nrms[nnc_tag];
    cov[nnc_tag].clear();
  } else {
    std::cout << "Improper nnc tag: " << nnc_tag << std::endl;
    return -1;
  }
  
  return EXIT_SUCCESS;
}
