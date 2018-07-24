/////////////////////////////////////////////////////////////////
// Copyright 2015 Gagarine Yaikhom (MIT License)
// https://github.com/gyaikhom/dbscan
// Modified by tjyang@fnal.gov
////////////////////////////////////////////////////////////////////
#ifndef DBSCAN3DALG_H
#define DBSCAN3DALG_H

#define UNCLASSIFIED -1
#define NOISE -2

#define CORE_POINT 1
#define NOT_CORE_POINT 0

#define SUCCESS 0
#define FAILURE -3

#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace recob { class SpacePoint; }

typedef struct point_s point_t;
struct point_s {
  art::Ptr<recob::SpacePoint> sp;
  int cluster_id;
};

typedef struct node_s node_t;
struct node_s {
  unsigned int index;
  node_t *next;
};

typedef struct epsilon_neighbours_s epsilon_neighbours_t;
struct epsilon_neighbours_s {
  unsigned int num_members;
  node_t *head, *tail;
};

namespace cluster{

  //--------------------------------------------------------------- 
  class DBScan3DAlg {
  public:
    
    
    DBScan3DAlg(fhicl::ParameterSet const& pset);
    virtual ~DBScan3DAlg();

    std::vector<point_t> points;

    void init(const std::vector<art::Ptr<recob::SpacePoint>>& sps);
    void dbscan();    
    
  private:

    double epsilon;
    unsigned int minpts;

    node_t *create_node(unsigned int index);
    int append_at_end(unsigned int index,
                      epsilon_neighbours_t *en);
    epsilon_neighbours_t *get_epsilon_neighbours(unsigned int index);
    void destroy_epsilon_neighbours(epsilon_neighbours_t *en);
    int expand(unsigned int index,
               unsigned int cluster_id);
    int spread(unsigned int index,
               epsilon_neighbours_t *seeds,
               unsigned int cluster_id);
    float dist(point_t *a, point_t *b);


  }; // class DBScan3DAlg
} // namespace

#endif // ifndef DBSCAN3DALG_H
