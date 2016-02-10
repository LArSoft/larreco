
#ifndef fuzzyClusterALG_H
#define fuzzyClusterALG_H

// Standard C/C++ libraries
#include <vector>
#include <set>
#include <stdint.h> // uint32_t

// ART and support libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Ptr.h"

// LArSoft libraries
#include "Geometry/Geometry.h"
#include "RecoAlg/HoughBaseAlg.h"
#include "RecoAlg/DBScanAlg.h"

namespace fhicl { class ParameterSet; }

namespace recob { class Hit; }

namespace cluster {
  
  //--------------------------------------------------------------- 
  class fuzzyClusterAlg {
  public:
    
    
    fuzzyClusterAlg(fhicl::ParameterSet const& pset);
    virtual ~fuzzyClusterAlg();
    
    void reconfigure(fhicl::ParameterSet const& p);
    void InitFuzzy(std::vector<art::Ptr<recob::Hit> >& allhits, std::set<uint32_t> badChannels);
    // Three differnt version of the clustering code
    void run_fuzzy_cluster(const std::vector<art::Ptr<recob::Hit> >& allhits);
   

    std::vector<std::vector<unsigned int> > fclusters;               ///< collection of something
    std::vector<std::vector<double> >       fps;                     ///< the collection of points we are working on     
    std::vector<unsigned int>               fpointId_to_clusterId;   ///< mapping point_id -> clusterId     
    std::vector<std::vector<double> >       fsim;                    ///<
    std::vector<std::vector<double> >       fsim2;                   ///<
    std::vector<std::vector<double> >       fsim3;                   ///<
    double fMaxWidth;

    // Get functions and structures from HoughBaseAlg
    //friend class HoughBaseAlg;

   

   //double **data = NULL;
   std::vector<std::vector<double>> data;




  private:
    
    // privately defined
    class baseCluster;
    class trackCluster;
    class showerCluster;


    double  fNumberTimeBoundaries;          ///< Number of boundaries in ticks for the drift window to be divided up to make the Hough line finder easier on memory
    double  fNumberWireBoundaries;          ///< Number of boundaries in wires for the drift window to be divided up to make the Hough line finder easier on memory

    // Run the Hough line finder?
    bool    fRunHough;
    bool    fGenerateHoughLinesOnly;            ///< Show only the results of the Hough line finder, hits not in a line will not be clustered

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














