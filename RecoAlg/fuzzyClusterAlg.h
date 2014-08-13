/////////////////////////////////////////////////////////////////
/*===================================================================
  The standard implementation of FLAME data clustering algorithm.

  FLAME (Fuzzy clustering by Local Approximation of MEmberships)
  was first described in:
  "FLAME, a novel fuzzy clustering method for the analysis of DNA
  microarray data", BMC Bioinformatics, 2007, 8:3.
  Available from: http://www.biomedcentral.com/1471-2105/8/3
  
  Copyright(C) 2007, Fu Limin (phoolimin@gmail.com).
  All rights reserved.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
  2. The origin of this software must not be misrepresented; you must 
     not claim that you wrote the original software. If you use this 
     software in a product, an acknowledgment in the product 
     documentation would be appreciated but is not required.
  3. Altered source versions must be plainly marked as such, and must
     not be misrepresented as being the original software.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
===================================================================*/
////////////////////////////////////////////////////////////////////

#ifndef fuzzyClusterALG_H
#define fuzzyClusterALG_H

// Standard C/C++ libraries
#include <vector>
#include <cmath>
#include <iostream>
#include <stdint.h> // uint32_t

// ROOT libraries
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVectorF.h"
#include "TVector.h"
#include "TH1.h"

// ART and support libraries
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArSoft libraries
#include "Geometry/Geometry.h"
#include "RecoAlg/HoughBaseAlg.h"
#include "RecoAlg/DBScanAlg.h"


namespace recob { class Hit; }

namespace cluster{

  namespace fuzzy_cluster {
    
    // This stores information about a showerlike cluster
    class baseCluster
      {
        public:
          int clusterNumber=-999999;
          std::vector<protoTrack> clusterProtoTracks;
          baseCluster(const protoTrack& protoTrackTemp)
            {
              clusterNumber=protoTrackTemp.clusterNumber;
              clusterProtoTracks.push_back(protoTrackTemp);
            }

          void addProtoTracks(const std::vector<protoTrack>& tracksToAdd)
            {
              // prepare to expand the list of tracks
              clusterProtoTracks.reserve(clusterProtoTracks.size() + tracksToAdd.size());
              for(const protoTrack& track: tracksToAdd) {
                protoTrack new_track(track);
                new_track.clusterNumber = clusterNumber;
                clusterProtoTracks.push_back(std::move(new_track));
              } // for
            } // addProtoTracks()
          
          void clearProtoTracks()
            { clusterProtoTracks.clear(); }

      }; // class baseCluster
  } // namespace fuzzy_cluster
  
  
  /// This class stores information about a showerlike cluster
  class showerCluster: public fuzzy_cluster::baseCluster {
      public:
    showerCluster(const protoTrack& protoTrackTemp):
      fuzzy_cluster::baseCluster(protoTrackTemp) {}
  }; // class showerCluster
  

  /// This class stores information about a tracklike cluster
  class trackCluster: public fuzzy_cluster::baseCluster {
      public:
    trackCluster(const protoTrack& protoTrackTemp):
      fuzzy_cluster::baseCluster(protoTrackTemp) {}
  }; // class trackCluster
  

  //--------------------------------------------------------------- 
  class fuzzyClusterAlg {
  public:
    
    
    fuzzyClusterAlg(fhicl::ParameterSet const& pset);
    virtual ~fuzzyClusterAlg();
    
    void reconfigure(fhicl::ParameterSet const& p);
    void InitFuzzy(std::vector<art::Ptr<recob::Hit> >& allhits, std::set<uint32_t> badChannels);
    // Three differnt version of the clustering code
    void run_fuzzy_cluster(std::vector<art::Ptr<recob::Hit> >& allhits);     
   

    std::vector<std::vector<unsigned int> > fclusters;               ///< collection of something
    std::vector<std::vector<double> >       fps;                     ///< the collection of points we are working on     
    std::vector<unsigned int>               fpointId_to_clusterId;   ///< mapping point_id -> clusterId     
    std::vector<std::vector<double> >       fsim;                    ///<
    std::vector<std::vector<double> >       fsim2;            	     ///<
    std::vector<std::vector<double> >       fsim3;            	     ///<
    double fMaxWidth;

    //Needed for Ben's FLAME cluster
    int NNumOfRows;
    int MNumOfCols;

    // Get functions and structures from HoughBaseAlg
    //friend class HoughBaseAlg;

   



   //double **data = NULL;
   std::vector<std::vector<double>> data;




  private:
   

    // The distance metric chosen for clustering, you likely do not need to modify this
    int fDistanceMetric;
    // The number of hits to be considered per k-nearest neighbors cluster 
    unsigned int fKNN;
    // The maximum number of iterations to try, needed for FLAME clustering
    int fIterations;
    // The limit in the difference between memberships when FLAME clustering stops
    double fEpsilon;
    // Sets the threshold parameter in FLAME cluster, it effectively sets a lower limit on hit density for whether a hit is considered part of a cluster or an outlier
    double fThreshold;


    // Run the Hough line finder?
    bool    fRunHough;


    bool    fDoFuzzyRemnantMerge;           ///< Tell the algorithm to merge fuzzy cluster remnants into showers or tracks (0-off, 1-on)
    double  fFuzzyRemnantMergeCutoff;       ///< cut off on merging the fuzzy cluster remnants into the nearest shower or track 

    bool    fDoTrackClusterMerge;           ///< Turn on cut on product of charge asymmetry and sin of angle between slopes of lines
    double  fTrackClusterMergeCutoff;          ///< Max distance between Hough lines before two lines are merged (muon tracks), 
    double  fChargeAsymAngleCut;            ///< Cut on product of charge asymmetry and sin of angle between slopes of lines
    double  fSigmaChargeAsymAngleCut;       ///< Cut on product of charge asymmetry and sin of angle between slopes of lines
  
    bool    fDoShowerClusterMerge;          ///< Turns on shower Hough line merging (0-off, 1-on)
    double  fShowerClusterMergeAngle;       ///< Max angle between slopes before two lines are merged, for lines in shower line regions
    double  fShowerClusterMergeCutoff;    ///< Max distance between Hough lines before two lines are merged (electron showers),

    bool    fDoShowerTrackClusterMerge;     ///< Turn on cut on product of charge asymmetry and sin of angle between slopes of lines
    double  fShowerTrackClusterMergeCutoff;    ///< Max distance between Hough lines before two lines are merged (electron showers),
    double  fShowerTrackClusterMergeAngle;  ///< Max angle between slopes before two lines are merged, for lines in shower line regions
    
    double  fShowerLikenessCut;             ///< Cut on shower likeness (larger the more shower like, smaller the less shower like)

    int    fMaxVertexLines;                ///< Max number of line end points allowed in a Hough line merge region for a merge to happen
    double  fVertexLinesCutoff;             ///< Size of the vertex region to count up lines for fMaxVertexLines 







    void mergeHoughLinesBySegment(unsigned int k,
        std::vector<protoTrack> *protoTracks, 
        double xyScale,
        int mergeStyle,
        double wire_dist,
        double tickToDist);

    bool mergeTrackClusters(unsigned int k,
        std::vector<trackCluster> *trackClusters, 
        double xyScale,
        double wire_dist,
        double tickToDist);

    bool mergeShowerClusters(unsigned int k,
        std::vector<showerCluster> *showerClusters, 
        double xyScale,
        double wire_dist,
        double tickToDist);
    
    bool mergeShowerTrackClusters(showerCluster *showerClusterI, 
        trackCluster *trackClusterJ, 
        double xyScale,
        double wire_dist,
        double tickToDist);

    //std::vector<lineSlope> linesFound;
    double HoughLineDistance(double p0MinLine1, 
        double p1MinLine1, 
        double p0MaxLine1, 
        double p1MaxLine1, 
        double p0MinLine2, 
        double p1MinLine2, 
        double p0MaxLine2, 
        double p1MaxLine2);
    bool   HoughLineIntersect(double x11,
        double  y11,
        double  x12,
        double  y12,
        double  x21,
        double  y21,
        double  x22,
        double  y22);
    double PointSegmentDistance(double px,
        double  py,
        double  x1,
        double  y1,
        double  x2,
        double  y2);

    double DistanceBetweenHits(
        art::Ptr<recob::Hit> hit0,
        art::Ptr<recob::Hit> hit1,
        double wire_dist,
        double tickToDist);

    // noise vector
    std::vector<bool>      fnoise;	
    std::vector<bool>      fvisited;					     
    std::vector<double>    fWirePitch;     ///< the pitch of the wires in each plane
    std::set<uint32_t>     fBadChannels;   ///< set of bad channels in this detector
    std::vector<uint32_t>  fBadWireSum;    ///< running total of bad channels. Used for fast intervening 
                                           ///< dead wire counting ala fBadChannelSum[m]-fBadChannelSum[n]. 

    


    //void IntArray_Push( IntArray *self, int value );


    // Object used for Hough transforms
    HoughBaseAlg fHBAlg;        

    // Object used for DBScan
    DBScanAlg fDBScan;        

    art::ServiceHandle<geo::Geometry> fGeom; ///< handle to geometry service

  }; // class fuzzyClusterAlg
    




} // namespace

#endif // ifndef fuzzyClusterALG_H














