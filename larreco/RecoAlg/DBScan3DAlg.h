/////////////////////////////////////////////////////////////////
//
// The MIT License (MIT)
//
// Copyright (c) 2015 Gagarine Yaikhom
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// https://github.com/gyaikhom/dbscan
//
// Modified by tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////
#ifndef DBSCAN3DALG_H
#define DBSCAN3DALG_H

#define UNCLASSIFIED -1
#define NOISE -2

#define CORE_POINT 1
#define NOT_CORE_POINT 0

#define SUCCESS 0
#define FAILURE -3

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "fhiclcpp/fwd.h"

#include <map>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // for WireID

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

typedef struct point_s point_t;
struct point_s {
  art::Ptr<recob::SpacePoint> sp;
  unsigned int nbadchannels;
  int cluster_id;
};

typedef struct node_s node_t;
struct node_s {
  unsigned int index;
  node_t* next;
};

typedef struct epsilon_neighbours_s epsilon_neighbours_t;
struct epsilon_neighbours_s {
  unsigned int num_members;
  node_t *head, *tail;
};

namespace cluster {

  //---------------------------------------------------------------
  class DBScan3DAlg {
  public:
    DBScan3DAlg(fhicl::ParameterSet const& pset);

    std::vector<point_t> points;

    void init(const std::vector<art::Ptr<recob::SpacePoint>>& sps,
              art::FindManyP<recob::Hit>& hitFromSp,
              art::Timestamp ts);
    void dbscan();

  private:
    double epsilon;
    unsigned int minpts;
    double badchannelweight;
    unsigned int neighbors;
    std::map<geo::WireID, int> badchannelmap;

    node_t* create_node(unsigned int index);
    int append_at_end(unsigned int index, epsilon_neighbours_t* en);
    epsilon_neighbours_t* get_epsilon_neighbours(unsigned int index);
    void destroy_epsilon_neighbours(epsilon_neighbours_t* en);
    int expand(unsigned int index, unsigned int cluster_id);
    int spread(unsigned int index, epsilon_neighbours_t* seeds, unsigned int cluster_id);
    float dist(point_t* a, point_t* b) const;

  }; // class DBScan3DAlg
} // namespace

#endif // ifndef DBSCAN3DALG_H
