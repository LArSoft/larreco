////////////////////////////////////////////////////////////////////////
//
// fuzzyClusterAlg.cxx
//
// Ben Carls, bcarls@fnal.gov
//
////////////////////////////////////////////////////////////////////////


// C/C++ standard library
#include <cmath> // std::sqrt(), std::atan(), std::pow()
#include <utility> // std::move()
#include <iterator> // std::back_inserter()
#include <algorithm> // std::remove_if()
#include <functional> // std::mem_fn()

// ART and support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/fuzzyClusterAlg.h"
#include "lardataobj/RecoBase/Hit.h"


/* Since data for clustering are usually noisy,
 * so it is not very necessary to have EPSILON extremely small.
 */
static constexpr double EPSILON = 1e-9;

namespace {
  template <typename T>
  inline constexpr T sqr(T v) { return v*v; }

  template <typename T>
  inline constexpr T sumsq(T a, T b) { return sqr(a) + sqr(b); }

  template <typename T>
  inline constexpr T norm(T a, T b) { return std::sqrt(sumsq(a, b)); }

} // local namespace



namespace cluster {
  
  // Define parameters that will tell us if we are doing a normal Hough line merge
  // or a shower Hough line merge
  enum class MergeMode: short int {
    Shower,
    Normal,
    ShowerIntercept,
    ChargeAsymAngle
  }; // MergeMode
  
  
  /// This stores information about a cluster
  class fuzzyClusterAlg::baseCluster {
    public:
    int clusterNumber=-999999;
    std::vector<protoTrack> clusterProtoTracks;
      
    baseCluster(protoTrack protoTrackTemp)
      {
        clusterNumber = protoTrackTemp.clusterNumber;
        clusterProtoTracks.emplace_back(std::move(protoTrackTemp));
      }
      
    void addProtoTracks(std::vector<protoTrack> tracksToAdd)
      {
        for(auto& trackToAdd: tracksToAdd)
          trackToAdd.clusterNumber = clusterNumber;
        clusterProtoTracks.reserve
          (clusterProtoTracks.size() + tracksToAdd.size());
        std::move(tracksToAdd.begin(),tracksToAdd.end(),
          std::back_inserter(clusterProtoTracks));
      } // addProtoTracks()

    void clearProtoTracks()
      {
        clusterProtoTracks.clear();
      } // clearProtoTracks()
    
  }; // class fuzzyCluster::baseCluster
  
  /// This stores information about a showerlike cluster
  class fuzzyClusterAlg::showerCluster: public fuzzyClusterAlg::baseCluster {
      public:
    showerCluster(protoTrack protoTrackTemp): baseCluster(protoTrackTemp) {}
  }; // class fuzzyClusterAlg::showerCluster
  
  /// This stores information about a tracklike cluster
  class fuzzyClusterAlg::trackCluster: public fuzzyClusterAlg::baseCluster {
      public:
    trackCluster(protoTrack protoTrackTemp): baseCluster(protoTrackTemp) {}
  }; // class fuzzyClusterAlg::trackCluster
  

} // namespace cluster




namespace cluster{
  const unsigned int kNO_CLUSTER    = UINT_MAX;
  const unsigned int kNOISE_CLUSTER = UINT_MAX-1;
}

//----------------------------------------------------------
// fuzzyClusterAlg stuff
//----------------------------------------------------------
cluster::fuzzyClusterAlg::fuzzyClusterAlg(fhicl::ParameterSet const& pset) 
   : fHBAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg")),
    fDBScan(pset.get< fhicl::ParameterSet >("DBScanAlg"))

{
  this->reconfigure(pset);
}

//----------------------------------------------------------
cluster::fuzzyClusterAlg::~fuzzyClusterAlg()
{
}

//----------------------------------------------------------
void cluster::fuzzyClusterAlg::reconfigure(fhicl::ParameterSet const& p)
{

  fNumberTimeBoundaries           = p.get< double >("NumberTimeBoundaries"           );
  fNumberWireBoundaries           = p.get< double >("NumberWireBoundaries"           );
  fDoFuzzyRemnantMerge            = p.get< bool   >("DoFuzzyRemnantMerge"            );
  fFuzzyRemnantMergeCutoff        = p.get< double >("FuzzyRemnantMergeCutoff"        );
  fDoTrackClusterMerge            = p.get< bool   >("DoTrackClusterMerge"            );
  fTrackClusterMergeCutoff        = p.get< double >("TrackClusterMergeCutoff"        );
  fChargeAsymAngleCut             = p.get< double >("ChargeAsymAngleCut"             );
  fSigmaChargeAsymAngleCut        = p.get< double >("SigmaChargeAsymAngleCut"        );
  fDoShowerClusterMerge           = p.get< bool   >("DoShowerClusterMerge"           );
  fDoShowerTrackClusterMerge      = p.get< bool   >("DoShowerTrackClusterMerge"      );
  fShowerClusterMergeCutoff       = p.get< double >("ShowerClusterMergeCutoff"       );
  fShowerClusterMergeAngle        = p.get< double >("ShowerClusterMergeAngle"        );
  fShowerTrackClusterMergeCutoff  = p.get< double >("ShowerTrackClusterMergeCutoff"  );
  fShowerTrackClusterMergeAngle   = p.get< double >("ShowerTrackClusterMergeAngle"   );
  fShowerLikenessCut              = p.get< double >("ShowerLikenessCut"              );
  fMaxVertexLines                 = p.get< int    >("MaxVertexLines"                 );
  fVertexLinesCutoff              = p.get< double >("VertexLinesCutoff"              );
  fRunHough                       = p.get< bool   >("RunHough"                       );
  fGenerateHoughLinesOnly         = p.get< bool   >("GenerateHoughLinesOnly"         );
  fHBAlg.reconfigure(p.get< fhicl::ParameterSet >("HoughBaseAlg"));
  fDBScan.reconfigure(p.get< fhicl::ParameterSet >("DBScanAlg"));
}


//----------------------------------------------------------
void cluster::fuzzyClusterAlg::InitFuzzy(std::vector<art::Ptr<recob::Hit> >& allhits, 
					 std::set<uint32_t>                  badChannels)
{
  // clear all the data member vectors for the new set of hits
  fps.clear();
  fpointId_to_clusterId.clear();
  fnoise.clear();
  fvisited.clear();
  fsim.clear();
  fsim2.clear();
  fsim3.clear();
  fclusters.clear();
  fWirePitch.clear();

  fBadChannels = badChannels;
  fBadWireSum.clear();

  const detinfo::DetectorProperties* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;


  unsigned int cs=allhits[0]->WireID().Cryostat;
  unsigned int t=allhits[0]->WireID().TPC;

  fWirePitch.resize(geom->Nplanes(t, cs), 0.);
  for(size_t p = 0; p < fWirePitch.size(); ++p) 
    fWirePitch[p] = geom->WirePitch(0,1,p);



  // Collect the hits in a useful form,
  // and take note of the maximum time width
  fMaxWidth=0.0;

  double tickToDist = detp->DriftVelocity(detp->Efield(),detp->Temperature());
  tickToDist *= 1.e-3 * detp->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
  int dims = 3;//our point is defined by 3 elements:wire#,center of the hit, and the hit width
  std::vector<double> p(dims);
  for (auto allhitsItr = allhits.begin(); allhitsItr < allhits.end(); allhitsItr++){
        
    p[0] = ((*allhitsItr)->WireID().Wire)*fWirePitch[(*allhitsItr)->WireID().Plane];
    p[1] = (*allhitsItr)->PeakTime()*tickToDist;
    p[2] = (*allhitsItr)->Integral();   //width of a hit in cm

    // check on the maximum width condition
    if ( p[2] > fMaxWidth ) fMaxWidth = p[2];
    
    fps.push_back(p);

  }

  mf::LogInfo("fuzzyCluster") << "Initfuzzy: hits vector size is " << fps.size();

  return;
}


//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
//
//  Ben Carls' implementation of fuzzyClusterAlg as much like examples as possible
void cluster::fuzzyClusterAlg::run_fuzzy_cluster(const std::vector<art::Ptr<recob::Hit> >& allhits) {

  // Don't attempt to run the algorithm if we have 1 or fewer hits
  if(allhits.size() <= 1)
    return;

  mf::LogInfo("fuzzyClusterAlg") << "Clustering " << allhits.size() << " hits";
  
  
  fpointId_to_clusterId.resize(fps.size(), kNO_CLUSTER); // Not zero as before!
  fnoise.resize(fps.size(), false);
  fvisited.resize(fps.size(), false);

  art::ServiceHandle<geo::Geometry> geom;
  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();


  //factor to make x and y scale the same units
//   uint32_t     channel = allhits[0]->Channel();
//   double wirePitch = geom->WirePitch(geom->View(channel));
    uint32_t     channel = allhits[0]->WireID().Wire;
//     double wirePitch[2];
//     wirePitch[0] = geom->WirePitch(0);
//     wirePitch[1] = geom->WirePitch(1);
//     wirePitch[2] = geom->WirePitch(2);

// saw this one

  double xyScale = .001*detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
  xyScale*=detprop->SamplingRate();
  //  double wire_dist = wirePitch;

  double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
  tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
  
  
  double indcolscaling = 0.;       //a parameter to account for the different 
  //characteristic hit width of induction and collection plane
  /// \todo: the collection plane's characteristic hit width's are, 
  /// \todo: on average, about 5 time samples wider than the induction plane's. 
  /// \todo: this is hard-coded for now.
  geo::SigType_t sigt = geom->SignalType(channel);
  if(sigt == geo::kInduction)
    indcolscaling = 0.;
  else
    indcolscaling = 1.;
  

  unsigned int nClusters = 1;
  unsigned int nClustersTemp = 1;

  // Divide up readout window into 9 sections that the Hough transform will be run over individually
  nClusters = fNumberTimeBoundaries*fNumberWireBoundaries;
  nClustersTemp = fNumberTimeBoundaries*fNumberWireBoundaries;
  unsigned int nWires = geom->Nwires(allhits[0]->WireID().Plane);
  unsigned int nTicks = detprop->ReadOutWindowSize();
  std::vector<unsigned int> wireBoundaries;
  std::vector<unsigned int> timeBoundaries;
  for(size_t i = 0; i < fNumberWireBoundaries+1; i++){
    wireBoundaries.push_back(i*nWires/fNumberWireBoundaries);
  }
  for(size_t i = 0; i < fNumberTimeBoundaries + 1; i++){
    timeBoundaries.push_back(i*nTicks/fNumberTimeBoundaries);
  }
  for(auto hitsItr = allhits.cbegin(); hitsItr != allhits.cend(); ++hitsItr){
    for(size_t i = 0; i < fNumberWireBoundaries; i++){
      if(( *hitsItr)->WireID().Wire >= wireBoundaries[i] && (*hitsItr)->WireID().Wire < wireBoundaries[i+1]){
        for(size_t j = 0; j < fNumberTimeBoundaries; j++){
          if((*hitsItr)->PeakTime() >= timeBoundaries[j] && (*hitsItr)->PeakTime() < timeBoundaries[j+1]){
            fpointId_to_clusterId[hitsItr-allhits.cbegin()] = j*fNumberTimeBoundaries+i;
          }
        }
      }
    }
  }






  // Loop over clusters with the Hough line finder to break the clusters up further
  // list of lines
  std::vector<protoTrack> protoTracksFound;
  if(nClustersTemp > 0 && fRunHough)
    for (unsigned int i = 0; i <= (unsigned int)nClustersTemp-1; ++i){
      LOG_DEBUG("fuzzyClusterAlg")
        << "Running Hough transform on protocluster " << i;
      fHBAlg.Transform(allhits, &fpointId_to_clusterId, i, &nClusters, &protoTracksFound);
    }

  // Determine the shower likeness of lines
  std::vector<showerCluster> showerClusters; 
  std::vector<trackCluster>  trackClusters; 
  double totalBkgDistCharge;
  double fMaxDistance;
  double distance;
  double peakTimePerpMin;
  double peakTimePerpMax;
  for(auto protoTracksFoundItr = protoTracksFound.begin(); protoTracksFoundItr < protoTracksFound.end(); ++protoTracksFoundItr){
    totalBkgDistCharge = 0;
    fMaxDistance = 0.1;
    for(auto hitsItr = allhits.cbegin(); hitsItr != allhits.cend(); ++hitsItr){
      /// Veto the hit if it already belongs to a line, proto tracks (Hough lines) are added after the fuzzy clusters
      //if(fpointId_to_clusterId.at(hitsItr-allhits.cbegin()) < nClustersTemp)
        //continue;
      distance = (std::abs((*hitsItr)->PeakTime()-protoTracksFoundItr->clusterSlope*(double)((*hitsItr)->WireID().Wire)-protoTracksFoundItr->clusterIntercept)/(std::sqrt(sqr(xyScale/fWirePitch[(*hitsItr)->WireID().Plane]*protoTracksFoundItr->clusterSlope)+1)));
      /// Sum up background hits, use smart distance
      peakTimePerpMin=-(1/protoTracksFoundItr->clusterSlope)*(double)((*hitsItr)->WireID().Wire)+allhits[protoTracksFoundItr->iMinWire]->PeakTime()+(1/protoTracksFoundItr->clusterSlope)*(allhits[protoTracksFoundItr->iMinWire]->WireID().Wire);
      peakTimePerpMax=-(1/protoTracksFoundItr->clusterSlope)*(double)((*hitsItr)->WireID().Wire)+allhits[protoTracksFoundItr->iMaxWire]->PeakTime()+(1/protoTracksFoundItr->clusterSlope)*(allhits[protoTracksFoundItr->iMaxWire]->WireID().Wire);
      if(distance > 1*(fMaxDistance+(*hitsItr)->RMS()+indcolscaling)
         && distance < 25*(fMaxDistance+(*hitsItr)->RMS()+indcolscaling)){
        if((protoTracksFoundItr->clusterSlope < 0 && (*hitsItr)->PeakTime() < peakTimePerpMin && (*hitsItr)->PeakTime() > peakTimePerpMax)
            || (protoTracksFoundItr->clusterSlope > 0 && (*hitsItr)->PeakTime() > peakTimePerpMin && (*hitsItr)->PeakTime() < peakTimePerpMax)){
          if((*hitsItr)->Integral() > 0.1){
            totalBkgDistCharge+=distance/(*hitsItr)->Integral();
          }
        }
      }
    }/// end loop over hits
    protoTracksFoundItr->showerLikeness = totalBkgDistCharge/(double)protoTracksFoundItr->hits.size();
    //std::cout << "showerLikeness: " << totalBkgDistCharge/(double)protoTracksFoundItr->hits.size() << std::endl;

    if(protoTracksFoundItr->showerLikeness > fShowerLikenessCut)
      showerClusters.push_back(showerCluster(*protoTracksFoundItr));
    else
      trackClusters.push_back(trackCluster(*protoTracksFoundItr));

  }/// end loop over lines found



  // Merge Hough lines
  bool trackMerged;
  bool showerMerged;
  bool showerTrackMerged;

  if(fDoTrackClusterMerge && trackClusters.size() > 1){
    unsigned int i = 0;
    while(i < trackClusters.size()-1){ 
      int ip = trackClusters[i].clusterProtoTracks[0].planeNumber;
      trackMerged = mergeTrackClusters(i,&trackClusters,xyScale/fWirePitch[ip],
				       fWirePitch[ip],tickToDist);
      if(trackMerged)
        continue;
      else
        i++;
    } 
  }

  if(fDoShowerClusterMerge && showerClusters.size() > 1){
    unsigned int i = 0;
    while(i < showerClusters.size()-1){ 
      int ip = showerClusters[i].clusterProtoTracks[0].planeNumber;
      showerMerged = mergeShowerClusters(i,&showerClusters,xyScale/fWirePitch[ip],
					 fWirePitch[ip],tickToDist);
      if(showerMerged)
        continue;
      else
        i++;
    } 
  }

  if(fDoShowerTrackClusterMerge && showerClusters.size() > 0 && trackClusters.size() >0){
    unsigned int i = 0;
    while(i < showerClusters.size()){ 
      int ip = showerClusters[i].clusterProtoTracks[0].planeNumber;
      unsigned int j = 0;
      while(j < trackClusters.size()){
	//      int ipt = trackClusters[j].clusterProtoTracks[0].planeNumber;
        showerTrackMerged = mergeShowerTrackClusters(&showerClusters[i],&trackClusters[j],
						     xyScale/fWirePitch[ip],
						     fWirePitch[ip],tickToDist);
        if(showerTrackMerged)
          continue;
        else
          j++;
      }
      i++;
    } 
  }

  if(fDoShowerClusterMerge && showerClusters.size() > 1){
    unsigned int i = 0;
    while(i < showerClusters.size()-1){ 
      int ip = showerClusters[i].clusterProtoTracks[0].planeNumber;
      showerMerged = mergeShowerClusters(i,&showerClusters,xyScale/fWirePitch[ip],
					 fWirePitch[ip],tickToDist);
      if(showerMerged)
        continue;
      else
        i++;
    } 
  }


  // Reassign the merged lines
  for(auto fpointId_to_clusterIdItr = fpointId_to_clusterId.begin(); fpointId_to_clusterIdItr != fpointId_to_clusterId.end(); ++fpointId_to_clusterIdItr){
    for(auto trackClustersItr = trackClusters.begin(); trackClustersItr != trackClusters.end(); ++trackClustersItr){
      for(auto protoTracksFoundItr = trackClustersItr->clusterProtoTracks.begin(); protoTracksFoundItr < trackClustersItr->clusterProtoTracks.end(); ++protoTracksFoundItr){
        if(*fpointId_to_clusterIdItr == (unsigned int)protoTracksFoundItr->oldClusterNumber)
          *fpointId_to_clusterIdItr = protoTracksFoundItr->clusterNumber;
      }
    }
    for(auto showerClustersItr = showerClusters.begin(); showerClustersItr != showerClusters.end(); ++showerClustersItr){
      for(auto protoTracksFoundItr = showerClustersItr->clusterProtoTracks.begin(); protoTracksFoundItr < showerClustersItr->clusterProtoTracks.end(); ++protoTracksFoundItr){
        if(*fpointId_to_clusterIdItr == (unsigned int)protoTracksFoundItr->oldClusterNumber){
          *fpointId_to_clusterIdItr = protoTracksFoundItr->clusterNumber;
        }
      }
    }
  }


  // Find sizes of all merged lines combined
  // For protoTracksFoundSizes, key is cluster number and size is the mapped value
  std::map<int,double> protoTracksFoundSizes;
  for(auto protoTracksFoundItr = protoTracksFound.begin(); protoTracksFoundItr < protoTracksFound.end(); ++protoTracksFoundItr){
    double d = ::norm(protoTracksFoundItr->pMin0-protoTracksFoundItr->pMax0, protoTracksFoundItr->pMin1-protoTracksFoundItr->pMax1);
    if(!protoTracksFoundSizes.count(protoTracksFoundItr->clusterNumber))
      protoTracksFoundSizes[protoTracksFoundItr->clusterNumber] = d;
    else 
      protoTracksFoundSizes[protoTracksFoundItr->clusterNumber]+= d;
  }



  // If only showing Hough lines, set the fuzzy cluster remnants to have no cluster
  if(fGenerateHoughLinesOnly){
    for(auto allhitsItr = allhits.cbegin(); allhitsItr != allhits.cend(); ++allhitsItr){
      if(fpointId_to_clusterId.at(allhitsItr-allhits.begin()) < (unsigned int) nClustersTemp )
        fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = kNO_CLUSTER;
      }
  }

  std::vector< art::Ptr<recob::Hit> > unclusteredhits;
  std::vector<unsigned int> unclusteredhitsToAllhits;
  int nDBClusters = 0;
  bool unclustered;
  double p0;
  double p1;
  double minDistanceShower;
  // double minDistanceTrack;
  if(fDoFuzzyRemnantMerge && !fGenerateHoughLinesOnly){
    for(auto allhitsItr = allhits.cbegin(); allhitsItr != allhits.cend(); ++allhitsItr){
      unclustered = true;
      // nClusters is the number of fuzzy clusters we found, we only assign hits to lines here
      // if they are not already part of hough lines
      if(fpointId_to_clusterId.at(allhitsItr-allhits.begin()) >= (unsigned int) nClustersTemp 
         && fpointId_to_clusterId.at(allhitsItr-allhits.begin()) < nClusters ){
        unclustered = false;
        continue;
      }
      int ip = (*allhitsItr)->WireID().Plane;
      p0 = ((*allhitsItr)->WireID().Wire)*fWirePitch[ip];
      p1 = (*allhitsItr)->PeakTime()*tickToDist;

      // First try to group it with a shower-like cluster
      minDistanceShower = 999999;
      for(auto showerClustersItr = showerClusters.begin(); showerClustersItr != showerClusters.end(); ++showerClustersItr){
        for(auto protoTracksItr = showerClustersItr->clusterProtoTracks.begin(); protoTracksItr < showerClustersItr->clusterProtoTracks.end(); ++protoTracksItr){
          
          distance = PointSegmentDistance( p0, p1, protoTracksItr->pMin0, protoTracksItr->pMin1, protoTracksItr->pMax0, protoTracksItr->pMax1);

          if(distance > fFuzzyRemnantMergeCutoff)
            continue;

        // This line used to have 1/4 instead of 0.25; wthat made it uneffective
        // (1/4 = 0); but replacing 0.25 sensibly changes performances;
        // commenting out the line for now [GP]
        //  distance /= std::pow(protoTracksFoundSizes[protoTracksItr->clusterNumber], 0.25);
          if(distance < minDistanceShower){
            fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = protoTracksItr->clusterNumber;
            minDistanceShower = distance;
            unclustered = false;
          }
        }
      }
      
      if(!unclustered)
        continue;
      
      // Failing to group it with a shower-like cluster, try with a track-like cluster
      // minDistanceTrack = 999999;
      // for(auto trackClustersItr = trackClusters.begin(); trackClustersItr != trackClusters.end(); ++trackClustersItr){
      //   for(auto protoTracksItr = trackClustersItr->clusterProtoTracks.begin(); protoTracksItr < trackClustersItr->clusterProtoTracks.end(); ++protoTracksItr){
          
      //     distance = PointSegmentDistance( p0, p1, protoTracksItr->pMin0, protoTracksItr->pMin1, protoTracksItr->pMax0, protoTracksItr->pMax1);

      //     if(distance > fFuzzyRemnantMergeCutoff)
      //       continue;

      //   // commenting out the line for now (see above) [GP]
      //   //  distance /= std::pow(protoTracksFoundSizes[protoTracksItr->clusterNumber], 0.25);
      //     if(distance < minDistanceTrack){
      //       fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = protoTracksItr->clusterNumber;
      //       minDistanceTrack = distance;
      //       unclustered = false;
      //     }
      //   }
      // }


      if(unclustered){
        unclusteredhitsToAllhits.push_back(allhitsItr-allhits.begin());
        unclusteredhits.push_back(*allhitsItr);
        fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = kNOISE_CLUSTER;
      }
      
    }

    // Setup DBSCAN for noise and extra hits
    // Start by getting the ChannelFilter
    // filter::ChannelFilter chanFilt;
    // fDBScan.InitScan(unclusteredhits, chanFilt.SetOfBadChannels());
    // fDBScan.run_cluster();
   
    // nDBClusters = fDBScan.fclusters.size();
    // for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){          
    //   if (fDBScan.fpointId_to_clusterId[j]== kNO_CLUSTER || fDBScan.fpointId_to_clusterId[j]==kNOISE_CLUSTER) {
    //     fpointId_to_clusterId.at(unclusteredhitsToAllhits[j]) = kNOISE_CLUSTER;
    //   } 
    //   else {
    //     fpointId_to_clusterId.at(unclusteredhitsToAllhits[j]) = fDBScan.fpointId_to_clusterId[j] + nClusters;
    //   }
    // }
  }

  size_t cid = nClusters + nDBClusters;
  
  //mf::LogInfo("fuzzyCluster") << "cid: " << cid ;

  //for(size_t j = 0; j < fpointId_to_clusterId.size(); ++j)
    //mf::LogInfo("fuzzyCluster") << "fpointId_to_clusterId[j]: " << fpointId_to_clusterId[j] << " j: " << j ;


  // Construct clusters, count noise, etc..
  int noise = 0;
  unsigned int c;
  fclusters.resize(cid);
  for(size_t y = 0; y < fpointId_to_clusterId.size(); ++y){
    if (fpointId_to_clusterId[y] == kNO_CLUSTER) {
      // This shouldn't happen...all points should be clasified by now!
      mf::LogWarning("fuzzyCluster") << "Unclassified point!";
    } 
    else if (fpointId_to_clusterId[y]==kNOISE_CLUSTER) {
      ++noise;
    } 
    else {
      c = fpointId_to_clusterId[y];
      if (c >= cid) {
	mf::LogWarning("fuzzyCluster") << "Point in cluster " << c 
			      << " when only " << cid 
			      << " clusters were found [0-" << cid-1
				 << "]";
      }
      fclusters[c].push_back(y);
    }
  }  
  mf::LogInfo("fuzzyCluster") << "DWM (R*-tree): Found " 
			   << cid << " clusters...";
  for (unsigned int c = 0; c < cid; ++c){
    mf::LogVerbatim("fuzzyCluster") << "\t" << "Cluster " << c << ":\t" 
			     << fclusters[c].size();
  }
  mf::LogVerbatim("fuzzyCluster") << "\t" << "...and " << noise << " noise points.";
}



// Merges based on the distance between line segments
bool cluster::fuzzyClusterAlg::mergeShowerTrackClusters(showerCluster *showerClusterI,
						        trackCluster *trackClusterJ,
						        double xyScale,
                                                        double wire_dist,
                                                        double tickToDist)
{

  // If we have zero or one Hough lines, move on 
  //if(trackCluster->size() == 0 || trackCluster->size() == 1)
    //return false;

  //// If we reach the last Hough line, move on 
  //if(trackCluster->size() == clusIndexStart+1)
    //return false;

  std::vector<unsigned int> toMerge; 
  std::vector<double> mergeSlope;
  std::vector<double> mergeTheta;

  // toMerge trackCluster index, toMerge trackCluster proto track index
  bool potentialBestMerge=false;
  bool performedBestMerge=false;
  unsigned int bestTrackClusterProtoTrack;
  unsigned int bestShowerClusterProtoTrack;
  // Did we merge left (0) or right (1)?
  int bestShowerRightLeft = -1;
  //int bestClusIndexStartRightLeft = -1;
  double bestToMergeTrackClusterProtoTrackDistance=999999;
  double x11; 
  double y11; 
  double x12; 
  double y12; 
  double x21; 
  double y21; 
  double x22; 
  double y22; 
  
  for(auto trackClusterProtoTrackItr = trackClusterJ->clusterProtoTracks.begin();
           trackClusterProtoTrackItr != trackClusterJ->clusterProtoTracks.end();
           trackClusterProtoTrackItr++){ 

    //for(auto trackClustersToMergeItr = trackClusters->begin()+clusIndexStart+1; trackClustersToMergeItr != trackClusters->end(); trackClustersToMergeItr++){
      //if(trackClusters->at(clusIndexStart).clusterNumber == trackClustersToMergeItr->clusterNumber)
        //continue;
      //std::cout << "Made it here" << std::endl;

      toMerge.clear();
      mergeSlope.clear();
      mergeTheta.clear();

      //Count up how many lines are in merging distance to clusIndexStartProtoTrackItr
      int nInDistanceTrackClusterLeft = 1;
      int nInDistanceTrackClusterRight = 1;



      for(auto showerClusterProtoTrackItr = showerClusterI->clusterProtoTracks.begin();
               showerClusterProtoTrackItr != showerClusterI->clusterProtoTracks.end();
               ++showerClusterProtoTrackItr){ 

        double segmentDistance = HoughLineDistance(trackClusterProtoTrackItr->pMin0,trackClusterProtoTrackItr->pMin1,
                                                   trackClusterProtoTrackItr->pMax0,trackClusterProtoTrackItr->pMax1, 
          					   showerClusterProtoTrackItr->pMin0,showerClusterProtoTrackItr->pMin1,
                                                   showerClusterProtoTrackItr->pMax0,showerClusterProtoTrackItr->pMax1);
        if(segmentDistance<fShowerTrackClusterMergeCutoff) 
        {
          toMerge.push_back(showerClusterProtoTrackItr-showerClusterI->clusterProtoTracks.begin());
          mergeSlope.push_back(trackClusterProtoTrackItr->clusterSlope*xyScale);
        
        
          // Sum up number of protoTracks at the vertex
          //distance between two segments in the plane:
          //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
          //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
          double x11 = showerClusterProtoTrackItr->pMin0; 
          double y11 = showerClusterProtoTrackItr->pMin1; 
          double x12 = showerClusterProtoTrackItr->pMax0; 
          double y12 = showerClusterProtoTrackItr->pMax1; 
          double x21 = trackClusterProtoTrackItr->pMin0; 
          double y21 = trackClusterProtoTrackItr->pMin1; 
          double x22 = trackClusterProtoTrackItr->pMax0; 
          double y22 = trackClusterProtoTrackItr->pMax1; 

          // Compare toMergerItr min with clusIndexStart max
          double mergeRightClusIndexStartDist = ::norm(x11-x22, y11-y22);
          // Compare toMergerItr max with clusIndexStart min
          double mergeLeftClusIndexStartDist = ::norm(x12-x21, y12-y21);
         
          // Are we inside the vertex distance? This is smaller than the merge cutoff
          if(segmentDistance < fVertexLinesCutoff){ 
            if( mergeRightClusIndexStartDist > mergeLeftClusIndexStartDist )
              ++nInDistanceTrackClusterLeft;
            else
              ++nInDistanceTrackClusterRight;
          }
        
        
        }

      }// End of loop over trackClustersToMergeItr->clusterProtoTracks.begin()


      mergeTheta.resize(toMerge.size());

      // Find the angle between the slopes
      for(auto mergeThetaItr = mergeTheta.begin(); mergeThetaItr != mergeTheta.end(); ++mergeThetaItr){
        double toMergeSlope = showerClusterI->clusterProtoTracks[toMerge[mergeThetaItr-mergeTheta.begin()]].clusterSlope*xyScale;
        mergeTheta[mergeThetaItr-mergeTheta.begin()] = std::atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
      }


      // Perform the merge
      for(auto toMergeItr = toMerge.begin(); toMergeItr != toMerge.end(); toMergeItr++){
        
        // Apply the angle cut
        if(mergeTheta[toMergeItr-toMerge.begin()] > fShowerTrackClusterMergeAngle)
          continue;

        // First check averages of charge and sigma charge for hits in lines closest to each other
        //int closestShower=-1;
        //int closestTrack=-1;
        double closestDistance=999999;
        for (auto showerClusterProtoTrackHitItr = showerClusterI->clusterProtoTracks[*toMergeItr].hits.begin(); showerClusterProtoTrackHitItr != showerClusterI->clusterProtoTracks[*toMergeItr].hits.end(); ++showerClusterProtoTrackHitItr) {
          for (auto trackClusterProtoTrackHitItr = trackClusterProtoTrackItr->hits.begin(); trackClusterProtoTrackHitItr != trackClusterProtoTrackItr->hits.end(); trackClusterProtoTrackHitItr++) {
            //double distance = ::norm(clusIndStHitItr->first-(*toMergeHitItr).first,
                      // clusIndStHitItr->second-toMergeHitItr->second);
            
            double distance = DistanceBetweenHits(*trackClusterProtoTrackHitItr,
                                                  *showerClusterProtoTrackHitItr,
                                                  wire_dist,
                                                  tickToDist);
            if(distance < closestDistance){
              closestDistance = distance;
              //closestShower=showerClusterProtoTrackHitItr-showerClusterI->clusterProtoTracks[*toMergeItr].hits.begin();
              //closestTrack=trackClusterProtoTrackHitItr-trackClusterProtoTrackItr->hits.begin();
            }
          }
        }


        // Veto the merge if the lines are not colinear 
      
        //distance between two segments in the plane:
        //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
        //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
        x11 = showerClusterI->clusterProtoTracks[*toMergeItr].pMin0; 
        y11 = showerClusterI->clusterProtoTracks[*toMergeItr].pMin1; 
        x12 = showerClusterI->clusterProtoTracks[*toMergeItr].pMax0; 
        y12 = showerClusterI->clusterProtoTracks[*toMergeItr].pMax1; 
        x21 = trackClusterProtoTrackItr->pMin0; 
        y21 = trackClusterProtoTrackItr->pMin1; 
        x22 = trackClusterProtoTrackItr->pMax0; 
        y22 = trackClusterProtoTrackItr->pMax1; 
        // BUG the double brace syntax is required to work around clang bug 21629
        // (https://bugs.llvm.org/show_bug.cgi?id=21629)
        std::array<double, 4> distances = {{
        // Compare toMergerItr min with clusIndexStart min, if this is the min distance, lines are not colinear, merge is vetoed
          ::norm(x11-x21, y11-y21),
          ::norm(x11-x22, y11-y22), // Compare toMergerItr min with clusIndexStart max
          ::norm(x12-x21, y12-y21), // Compare toMergerItr max with clusIndexStart min
        // Compare toMergerItr max with clusIndexStart max, if this is the min distance, lines are not colinear, merge is vetoed
          ::norm(x12-x22, y12-y22)
        }}; // distances

        double minDistance = 999999; 
        int minDistanceIndex = -1;
        for(unsigned int j = 0; j < distances.size(); j++){
          if (distances[j] < minDistance){
            minDistance = distances[j];
            minDistanceIndex = j;
          }
        }

        if(minDistanceIndex  == 0 || minDistanceIndex  == 3)
          continue;

        // How many lines do we have at the merging point? If too many, veto the merge
        //std::cout << nInDistanceTrackClusterLeft << " " << nInDistanceTrackClusterRight << std::endl;

        if(nInDistanceTrackClusterLeft > fMaxVertexLines) {
          trackClusterProtoTrackItr->mergedLeft = true;
          showerClusterI->clusterProtoTracks[*toMergeItr].mergedRight = true;
        }
        if(nInDistanceTrackClusterRight > fMaxVertexLines) {
          trackClusterProtoTrackItr->mergedRight = true;
          showerClusterI->clusterProtoTracks[*toMergeItr].mergedLeft = true;
        }
          

        // Check if we merged left or right already for clusIndexStart, we only do it once for each side
        //if(trackClusterProtoTrackItr->mergedLeft == true && minDistanceIndex == 2)
          //continue;
        //if(trackClusterProtoTrackItr->mergedRight == true && minDistanceIndex == 1)
          //continue;
        //if(showerClusterI->clusterProtoTracks[*toMergeItr].mergedLeft == true && minDistanceIndex == 1)
          //continue;
        //if(showerClusterI->clusterProtoTracks[*toMergeItr].mergedRight == true && minDistanceIndex == 2)
          //continue;

        //std::cout << "Potential merge" << std::endl;
        //std::cout << "main trackClustersItr slope: " << trackClusters->at(*toMergeItr).clusterSlope << " clusIndexStart slope: " << trackClusters->at(clusIndexStart).clusterSlope << std::endl;
        potentialBestMerge=true;     


        if(minDistance < bestToMergeTrackClusterProtoTrackDistance){
          bestShowerClusterProtoTrack=*toMergeItr;
          bestTrackClusterProtoTrack=trackClusterProtoTrackItr-trackClusterJ->clusterProtoTracks.begin();
          
          bestToMergeTrackClusterProtoTrackDistance=minDistance;
         
          // Did we merge left (0) or right (1)?
          if(minDistanceIndex == 1){
            bestShowerRightLeft = 0;
            //bestClusIndexStartRightLeft = 1;
          }
          if(minDistanceIndex == 2){
            bestShowerRightLeft = 1;
            //bestClusIndexStartRightLeft = 0;
          }

        }

      }// End of loop over toMerge
    //}// End of loop over trackClusters->begin()+clusIndexStart+1
  }//End of loop over trackClusters->at(clusIndexStart).clusterProtoTracks

  if(potentialBestMerge){
    showerClusterI->clusterProtoTracks[bestShowerClusterProtoTrack].merged=true;
    trackClusterJ->clusterProtoTracks[bestTrackClusterProtoTrack].merged=true;   
    if(bestShowerRightLeft == 0){
      showerClusterI->clusterProtoTracks[bestShowerClusterProtoTrack].mergedLeft = true;
      trackClusterJ->clusterProtoTracks[bestTrackClusterProtoTrack].mergedRight = true;   
      performedBestMerge=true;
    }
    if(bestShowerRightLeft == 1){
      showerClusterI->clusterProtoTracks[bestShowerClusterProtoTrack].mergedRight = true;
      trackClusterJ->clusterProtoTracks[bestTrackClusterProtoTrack].mergedLeft = true;   
      performedBestMerge=true;
    }
   
    if(performedBestMerge){ 
      showerClusterI->addProtoTracks(trackClusterJ->clusterProtoTracks);
      trackClusterJ->clearProtoTracks();
      //std::cout << "Merged shower-track" << std::endl;
    }

  }


  
  return performedBestMerge;

}



// Merges based on the distance between line segments
bool cluster::fuzzyClusterAlg::mergeTrackClusters(unsigned int clusIndexStart,
						     std::vector<trackCluster> *trackClusters,
						     double xyScale,
                                                     double wire_dist,
                                                     double tickToDist)
{



  // If we have zero or one Hough lines, move on 
  if(trackClusters->size() == 0 || trackClusters->size() == 1)
    return false;

  // If we reach the last Hough line, move on 
  if(trackClusters->size() == clusIndexStart+1)
    return false;

  std::vector<unsigned int> toMerge; 
  std::vector<double> mergeSlope;
  std::vector<double> mergeTheta;

  // toMerge trackCluster index, toMerge trackCluster proto track index
  bool potentialBestMerge=false;
  bool performedBestMerge=false;
  unsigned int bestToMergeTrackCluster;
  unsigned int bestTrackClustersClusIndexStartProtoTrack;
  unsigned int bestToMergeTrackClusterProtoTrack;
  // Did we merge left (0) or right (1)?
  int bestToMergeRightLeft = -1;
  //int bestClusIndexStartRightLeft = -1;
  double bestToMergeTrackClusterProtoTrackDistance=999999;

  for(auto trackClustersClusIndexStartProtoTrackItr = trackClusters->at(clusIndexStart).clusterProtoTracks.begin();
           trackClustersClusIndexStartProtoTrackItr != trackClusters->at(clusIndexStart).clusterProtoTracks.end();
           ++trackClustersClusIndexStartProtoTrackItr){ 

    //Count up how many lines are in merging distance to clusIndexStartProtoTrackItr
    int nInDistanceClusIndexStartLeft = 1;
    int nInDistanceClusIndexStartRight = 1;

    
    for(auto trackClustersToMergeItr = trackClusters->begin()+clusIndexStart+1; trackClustersToMergeItr != trackClusters->end(); trackClustersToMergeItr++){
      if(trackClusters->at(clusIndexStart).clusterNumber == trackClustersToMergeItr->clusterNumber)
        continue;
      //std::cout << "Made it here" << std::endl;

      toMerge.clear();
      mergeSlope.clear();
      mergeTheta.clear();

      for(auto trackClustersToMergeProtoTrackItr = trackClustersToMergeItr->clusterProtoTracks.begin();
               trackClustersToMergeProtoTrackItr != trackClustersToMergeItr->clusterProtoTracks.end();
               ++trackClustersToMergeProtoTrackItr){ 

        double segmentDistance = HoughLineDistance(trackClustersClusIndexStartProtoTrackItr->pMin0,trackClustersClusIndexStartProtoTrackItr->pMin1,
                                                   trackClustersClusIndexStartProtoTrackItr->pMax0,trackClustersClusIndexStartProtoTrackItr->pMax1, 
          					   trackClustersToMergeProtoTrackItr->pMin0,trackClustersToMergeProtoTrackItr->pMin1,
                                                   trackClustersToMergeProtoTrackItr->pMax0,trackClustersToMergeProtoTrackItr->pMax1);
        if(segmentDistance<fTrackClusterMergeCutoff) 
        {
          toMerge.push_back(trackClustersToMergeProtoTrackItr-trackClustersToMergeItr->clusterProtoTracks.begin());
          mergeSlope.push_back(trackClustersClusIndexStartProtoTrackItr->clusterSlope*xyScale);
       

          // Sum up number of protoTracks at the vertex
          //distance between two segments in the plane:
          //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
          //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
          double x11 = trackClustersToMergeProtoTrackItr->pMin0; 
          double y11 = trackClustersToMergeProtoTrackItr->pMin1; 
          double x12 = trackClustersToMergeProtoTrackItr->pMax0; 
          double y12 = trackClustersToMergeProtoTrackItr->pMax1; 
          double x21 = trackClustersClusIndexStartProtoTrackItr->pMin0; 
          double y21 = trackClustersClusIndexStartProtoTrackItr->pMin1; 
          double x22 = trackClustersClusIndexStartProtoTrackItr->pMax0; 
          double y22 = trackClustersClusIndexStartProtoTrackItr->pMax1; 

          // Compare toMergerItr min with clusIndexStart max
          double mergeRightClusIndexStartDist = std::sqrt(sqr(x11-x22) + sqr(y11-y22));
          // Compare toMergerItr max with clusIndexStart min
          double mergeLeftClusIndexStartDist = std::sqrt(sqr(x12-x21) + sqr(y12-y21));
          
          // Are we inside the vertex distance? This is smaller than the merge cutoff
          if(segmentDistance < fVertexLinesCutoff){ 
            if( mergeRightClusIndexStartDist > mergeLeftClusIndexStartDist )
              nInDistanceClusIndexStartLeft++;
            else
              nInDistanceClusIndexStartRight++;
          } 
        }

      }// End of loop over trackClustersToMergeItr->clusterProtoTracks.begin()


      mergeTheta.resize(toMerge.size());

      // Find the angle between the slopes
      for(auto mergeThetaItr = mergeTheta.begin(); mergeThetaItr != mergeTheta.end(); ++mergeThetaItr){
        double toMergeSlope = trackClustersToMergeItr->clusterProtoTracks[toMerge[mergeThetaItr-mergeTheta.begin()]].clusterSlope*xyScale;
        mergeTheta[mergeThetaItr-mergeTheta.begin()] = std::atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
      }


      // Perform the merge
      for(auto toMergeItr = toMerge.begin(); toMergeItr != toMerge.end(); toMergeItr++){

        // First check averages of charge and sigma charge for hits in lines closest to each other
        int closestToMerge=-1;
        int closestClusIndexStart=-1;
        double closestDistance=999999;
        for (auto toMergeHitItr = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin(); toMergeHitItr != trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.end(); toMergeHitItr++) {
          for (auto clusIndStHitItr = trackClustersClusIndexStartProtoTrackItr->hits.begin(); clusIndStHitItr != trackClustersClusIndexStartProtoTrackItr->hits.end(); clusIndStHitItr++) {
            //double distance = std::sqrt(sqr(clusIndStHitItr->first-(*toMergeHitItr).first)+
                      //sqr(clusIndStHitItr->second-toMergeHitItr->second));
            
            double distance = DistanceBetweenHits(*clusIndStHitItr,
                                                    *toMergeHitItr,
                                                    wire_dist,
                                                    tickToDist);
            if(distance < closestDistance){
              closestDistance = distance;
              closestToMerge=toMergeHitItr-trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin();
              closestClusIndexStart=clusIndStHitItr-trackClustersClusIndexStartProtoTrackItr->hits.begin();
            }
          }
        }

        // Find up to 9 more points closest to closestToMerge on the toMerge[i] line
        // check if it's closer, insert, delete
        using Pair_I_D = std::pair<int,double>;
        std::vector<Pair_I_D> closestToMergeDist;
        for (auto toMergeHitItr = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin(); toMergeHitItr != trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.end(); ++toMergeHitItr) {
          if(closestToMerge==toMergeHitItr-trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin())
            continue;
          double distance = DistanceBetweenHits(trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStart],
                                                  *toMergeHitItr,
                                                  wire_dist,
                                                  tickToDist);

          bool foundCloser = false;
          for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end(); closestToMergeDistItr++) {
            if(closestToMergeDistItr->second > distance){
                foundCloser = true;
                break;
              }
            }
            if(foundCloser 
                || closestToMergeDist.size() < trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.size()-1
                || closestToMergeDist.size() < 9){
              closestToMergeDist.emplace_back(toMergeHitItr-trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin(),distance);
              std::sort(closestToMergeDist.begin(), closestToMergeDist.end(),
                [](const Pair_I_D& a, const Pair_I_D& b){ return a.second < b.second; }
                );
            //  std::sort(closestToMergeDist.begin(), closestToMergeDist.end(), boost::bind(&std::pair<int,double>::second,_1) < boost::bind(&std::pair<int,double>::second,_2));
            }
            if(closestToMergeDist.size() > trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.size()-1 ||
              closestToMergeDist.size() > 9)
              closestToMergeDist.erase(closestToMergeDist.end());
        }
        //for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end();
          //closestToMergeDistItr++) 
          //std::cout << closestToMergeDistItr->first << " " << closestToMergeDistItr->second << std::endl;



        // Find up to 9 more points closest to closestToMerge on the clusIndexStart line
        std::vector<Pair_I_D> closestClusIndexStartDist;
        for (auto clusIndexStartHitItr = trackClustersClusIndexStartProtoTrackItr->hits.begin(); clusIndexStartHitItr != trackClustersClusIndexStartProtoTrackItr->hits.end(); ++clusIndexStartHitItr) {
          if(closestClusIndexStart==clusIndexStartHitItr-trackClustersClusIndexStartProtoTrackItr->hits.begin())
            continue;

          double distance = DistanceBetweenHits(*clusIndexStartHitItr,
                                                trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMerge],
                                                wire_dist,
                                                tickToDist);

          bool foundCloser = false;
          for(auto closestClusIndexStartDistItr = closestClusIndexStartDist.begin(); closestClusIndexStartDistItr != closestClusIndexStartDist.end(); ++closestClusIndexStartDistItr) {
            if(closestClusIndexStartDistItr->second > distance){
              foundCloser = true;
              break;
            }
          }
          if(foundCloser 
              || closestClusIndexStartDist.size() < trackClustersClusIndexStartProtoTrackItr->hits.size()-1
              || closestClusIndexStartDist.size() < 9){
            closestClusIndexStartDist.emplace_back(clusIndexStartHitItr-trackClustersClusIndexStartProtoTrackItr->hits.begin(),distance);
            std::sort(closestClusIndexStartDist.begin(), closestClusIndexStartDist.end(),
              [](const Pair_I_D& a, const Pair_I_D& b){ return a.second < b.second; }
              );
          }
          if(closestClusIndexStartDist.size() > trackClustersClusIndexStartProtoTrackItr->hits.size()-1 ||
            closestClusIndexStartDist.size() > 9)
            closestClusIndexStartDist.erase(closestClusIndexStartDist.end());
       }



        double toMergeAveCharge = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMerge]->Integral();
        double toMergeAveSigmaCharge = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMerge]->SigmaIntegral();
        for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end(); ++closestToMergeDistItr) {
          toMergeAveCharge+= trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMergeDistItr->first]->Integral();
          toMergeAveSigmaCharge+= trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMergeDistItr->first]->SigmaIntegral();
        }
        double clusIndexStartAveCharge = trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStart]->Integral();
        double clusIndexStartAveSigmaCharge = trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStart]->SigmaIntegral();
        for(auto closestClusIndexStartDistItr = closestClusIndexStartDist.begin(); closestClusIndexStartDistItr != closestClusIndexStartDist.end(); ++closestClusIndexStartDistItr) {
          clusIndexStartAveCharge+= trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStartDistItr->first]->Integral();
          clusIndexStartAveSigmaCharge+=trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStartDistItr->first]->SigmaIntegral();
        }



        double chargeAsymmetry = std::abs(toMergeAveCharge-clusIndexStartAveCharge)/(toMergeAveCharge+clusIndexStartAveCharge);
        double sigmaChargeAsymmetry = std::abs(toMergeAveSigmaCharge-clusIndexStartAveSigmaCharge)/(toMergeAveSigmaCharge+clusIndexStartAveSigmaCharge);
        double chargeAsymmetrySinAngle = chargeAsymmetry*std::abs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180));
        double sigmaChargeAsymmetrySinAngle = sigmaChargeAsymmetry*std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180));

        //std::cout << std::endl;
        //std::cout << chargeAsymmetry*std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)) << std::endl;
        //std::cout << sigmaChargeAsymmetry*std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)) << std::endl;
        

        if(chargeAsymmetrySinAngle > fChargeAsymAngleCut)
          continue;

        if(sigmaChargeAsymmetrySinAngle > fSigmaChargeAsymAngleCut)
          continue;


        // Veto the merge if the lines are not colinear 
      
        //distance between two segments in the plane:
        //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
        //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
        double x11 = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].pMin0; 
        double y11 = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].pMin1; 
        double x12 = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].pMax0; 
        double y12 = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].pMax1; 
        double x21 = trackClustersClusIndexStartProtoTrackItr->pMin0; 
        double y21 = trackClustersClusIndexStartProtoTrackItr->pMin1; 
        double x22 = trackClustersClusIndexStartProtoTrackItr->pMax0; 
        double y22 = trackClustersClusIndexStartProtoTrackItr->pMax1; 
        // BUG the double brace syntax is required to work around clang bug 21629
        // (https://bugs.llvm.org/show_bug.cgi?id=21629)
        std::array<double, 4> distances = {{
        // Compare toMergerItr min with clusIndexStart min, if this is the min distance, lines are not colinear, merge is vetoed
          ::norm(x11-x21, y11-y21),
          ::norm(x11-x22, y11-y22), // Compare toMergerItr min with clusIndexStart max
          ::norm(x12-x21, y12-y21), // Compare toMergerItr max with clusIndexStart min
        // Compare toMergerItr max with clusIndexStart max, if this is the min distance, lines are not colinear, merge is vetoed
          ::norm(x12-x22, y12-y22)
        }}; // distances

        double minDistance = 999999; 
        int minDistanceIndex = -1;
        for(unsigned int j = 0; j < distances.size(); j++){
          if (distances[j] < minDistance){
            minDistance = distances[j];
            minDistanceIndex = j;
          }
        }

        if(minDistanceIndex  == 0 || minDistanceIndex  == 3)
          continue;


        // How many lines do we have at the merging point? If too many, veto the merge
        //std::cout << nInDistanceClusIndexStartLeft << " " << nInDistanceClusIndexStartRight << std::endl;

        if(nInDistanceClusIndexStartLeft > fMaxVertexLines) {
          trackClustersClusIndexStartProtoTrackItr->mergedLeft = true;
          trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].mergedRight = true;
        }
        if(nInDistanceClusIndexStartRight > fMaxVertexLines) {
          trackClustersClusIndexStartProtoTrackItr->mergedRight = true;
          trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].mergedLeft = true;
        }




        // Check if we merged left or right already for clusIndexStart, we only do it once for each side
        if(trackClustersClusIndexStartProtoTrackItr->mergedLeft == true && minDistanceIndex == 2)
          continue;
        if(trackClustersClusIndexStartProtoTrackItr->mergedRight == true && minDistanceIndex == 1)
          continue;
        if(trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].mergedLeft == true && minDistanceIndex == 1)
          continue;
        if(trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].mergedRight == true && minDistanceIndex == 2)
          continue;

        //std::cout << "Potential merge" << std::endl;
        //std::cout << "main trackClustersItr slope: " << trackClusters->at(*toMergeItr).clusterSlope << " clusIndexStart slope: " << trackClusters->at(clusIndexStart).clusterSlope << std::endl;
        potentialBestMerge=true;     


        if(minDistance < bestToMergeTrackClusterProtoTrackDistance){
          bestToMergeTrackCluster=trackClustersToMergeItr-trackClusters->begin();
          bestToMergeTrackClusterProtoTrack=*toMergeItr;
          bestTrackClustersClusIndexStartProtoTrack=trackClustersClusIndexStartProtoTrackItr-trackClusters->at(clusIndexStart).clusterProtoTracks.begin();
          bestToMergeTrackClusterProtoTrackDistance=minDistance;

          // Did we merge left (0) or right (1)?
          if(minDistanceIndex == 1){
            bestToMergeRightLeft = 0;
            //bestClusIndexStartRightLeft = 1;
          }
          if(minDistanceIndex == 2){
            bestToMergeRightLeft = 1;
            //bestClusIndexStartRightLeft = 0;
          }

        }

      }// End of loop over toMerge
    }// End of loop over trackClusters->begin()+clusIndexStart+1
  }//End of loop over trackClusters->at(clusIndexStart).clusterProtoTracks

  if(potentialBestMerge){
    trackClusters->at(bestToMergeTrackCluster).clusterProtoTracks[bestToMergeTrackClusterProtoTrack].merged=true;
    trackClusters->at(clusIndexStart).clusterProtoTracks[bestTrackClustersClusIndexStartProtoTrack].merged=true;   
    if(bestToMergeRightLeft == 0){
      trackClusters->at(bestToMergeTrackCluster).clusterProtoTracks[bestToMergeTrackClusterProtoTrack].mergedLeft = true;
      trackClusters->at(clusIndexStart).clusterProtoTracks[bestTrackClustersClusIndexStartProtoTrack].mergedRight = true;   
      performedBestMerge=true;
    }
    if(bestToMergeRightLeft == 1){
      trackClusters->at(bestToMergeTrackCluster).clusterProtoTracks[bestToMergeTrackClusterProtoTrack].mergedRight = true;
      trackClusters->at(clusIndexStart).clusterProtoTracks[bestTrackClustersClusIndexStartProtoTrack].mergedLeft = true;   
      performedBestMerge=true;
    }
   
    if(performedBestMerge){ 
      trackClusters->at(clusIndexStart).addProtoTracks(trackClusters->at(bestToMergeTrackCluster).clusterProtoTracks);
      trackClusters->at(bestToMergeTrackCluster).clearProtoTracks();
      //std::cout << "Merged track-track" << std::endl;
    }

  }

  //lineMerged = true;
  //trackClustersClusIndexStartProtoTrackItr->merged = true;
  //trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].merged = true;

  //// For loop over all lines found to reassign lines to clusIndexStart that already belonged to toMerge 
  //// Need to delete trackClustersItr that gets merged, load protoTracks from one to the other 
  ////
  ////
  //for(auto trackClustersItr = trackClusters->begin(); trackClustersItr != trackClusters->end(); trackClustersItr++){
    //if((unsigned int)(*toMergeItr) == trackClustersItr-trackClusters->begin())
      //continue;

    //if(trackClustersItr->clusterNumber == trackClusters->at(*toMergeItr).clusterNumber){
      //trackClustersItr->clusterNumber = trackClusters->at(clusIndexStart).clusterNumber;
    //}
  //}
  //trackClusters->at(*toMergeItr).clusterNumber = trackClusters->at(clusIndexStart).clusterNumber;
  
  
  return performedBestMerge;

}



// Merges based on the distance between line segments
bool cluster::fuzzyClusterAlg::mergeShowerClusters(unsigned int clusIndexStart,
						     std::vector<showerCluster> *showerClusters,
						     double xyScale,
                                                     double wire_dist,
                                                     double tickToDist)
{

  // If we have zero or one Hough lines, move on 
  if(showerClusters->size() == 0 || showerClusters->size() == 1)
    return false;

  // If we reach the last Hough line, move on 
  if(showerClusters->size() == clusIndexStart+1)
    return false;

  std::vector<unsigned int> toMerge; 
  std::vector<double> mergeSlope;
  std::vector<double> mergeTheta;

  // toMerge trackCluster index, toMerge trackCluster proto track index
  //bool potentialBestMerge=false;
  bool performedBestMerge=false;
  unsigned int bestToMergeShowerCluster;
  unsigned int bestShowerClustersClusIndexStartProtoTrack;
  unsigned int bestToMergeShowerClusterProtoTrack;
  int bestToMergeRightLeft = -1;
  //int bestClusIndexStartRightLeft = -1;
  double bestToMergeShowerClusterProtoTrackDistance=999999;


  for(auto showerClustersClusIndexStartProtoTrackItr = showerClusters->at(clusIndexStart).clusterProtoTracks.begin();
           showerClustersClusIndexStartProtoTrackItr != showerClusters->at(clusIndexStart).clusterProtoTracks.end();
           ++showerClustersClusIndexStartProtoTrackItr){ 

    for(auto showerClustersToMergeItr = showerClusters->begin()+clusIndexStart+1; showerClustersToMergeItr != showerClusters->end(); ++showerClustersToMergeItr){
      //std::cout << "Made it here" << std::endl;

      toMerge.clear();
      mergeSlope.clear();
      mergeTheta.clear();

      for(auto showerClustersToMergeProtoTrackItr = showerClustersToMergeItr->clusterProtoTracks.begin();
               showerClustersToMergeProtoTrackItr != showerClustersToMergeItr->clusterProtoTracks.end();
               ++showerClustersToMergeProtoTrackItr){ 

        if(showerClustersToMergeProtoTrackItr->clusterNumber == showerClustersClusIndexStartProtoTrackItr->clusterNumber)
          continue;


        double segmentDistance = HoughLineDistance(showerClustersClusIndexStartProtoTrackItr->pMin0,showerClustersClusIndexStartProtoTrackItr->pMin1,
                                                   showerClustersClusIndexStartProtoTrackItr->pMax0,showerClustersClusIndexStartProtoTrackItr->pMax1, 
          					   showerClustersToMergeProtoTrackItr->pMin0,showerClustersToMergeProtoTrackItr->pMin1,
                                                   showerClustersToMergeProtoTrackItr->pMax0,showerClustersToMergeProtoTrackItr->pMax1);
        if(segmentDistance<fShowerClusterMergeCutoff) 
        {
          toMerge.push_back(showerClustersToMergeProtoTrackItr-showerClustersToMergeItr->clusterProtoTracks.begin());
          mergeSlope.push_back(showerClustersClusIndexStartProtoTrackItr->clusterSlope*xyScale);
        }

      }// End of loop over showerClustersToMergeItr->clusterProtoTracks.begin()


      mergeTheta.resize(toMerge.size());

      // Find the angle between the slopes
      for(auto mergeThetaItr = mergeTheta.begin(); mergeThetaItr != mergeTheta.end(); ++mergeThetaItr){
        double toMergeSlope = showerClustersToMergeItr->clusterProtoTracks[toMerge[mergeThetaItr-mergeTheta.begin()]].clusterSlope*xyScale;
        mergeTheta[mergeThetaItr-mergeTheta.begin()] = std::atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
      }


      // Perform the merge
      for(auto toMergeItr = toMerge.begin(); toMergeItr != toMerge.end(); ++toMergeItr){

        // Apply the angle cut
        if(mergeTheta[toMergeItr-toMerge.begin()] > fShowerClusterMergeAngle)
          continue;

        // Find the closest distance 
        //int closestToMerge=-1;
        //int closestClusIndexStart=-1;
        double closestDistance=999999;
        for (auto toMergeHitItr = showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin(); toMergeHitItr != showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.end(); ++toMergeHitItr) {
          for (auto clusIndStHitItr = showerClustersClusIndexStartProtoTrackItr->hits.begin(); clusIndStHitItr != showerClustersClusIndexStartProtoTrackItr->hits.end(); ++clusIndStHitItr) {
            //double distance = ::norm(clusIndStHitItr->first-(*toMergeHitItr).first+
                      //clusIndStHitItr->second-toMergeHitItr->second);
            
            double distance = DistanceBetweenHits(*clusIndStHitItr,
                                                    *toMergeHitItr,
                                                    wire_dist,
                                                    tickToDist);
            if(distance < closestDistance){
              closestDistance = distance;
              //closestToMerge=toMergeHitItr-showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin();
              //closestClusIndexStart=clusIndStHitItr-showerClustersClusIndexStartProtoTrackItr->hits.begin();
            }
          }
        }


        // Veto the merge if the lines are not colinear 
      
        //distance between two segments in the plane:
        //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
        //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
        double x11 = showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].pMin0; 
        double y11 = showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].pMin1; 
        double x12 = showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].pMax0; 
        double y12 = showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].pMax1; 
        double x21 = showerClustersClusIndexStartProtoTrackItr->pMin0; 
        double y21 = showerClustersClusIndexStartProtoTrackItr->pMin1; 
        double x22 = showerClustersClusIndexStartProtoTrackItr->pMax0; 
        double y22 = showerClustersClusIndexStartProtoTrackItr->pMax1; 
        // BUG the double brace syntax is required to work around clang bug 21629
        // (https://bugs.llvm.org/show_bug.cgi?id=21629)
        std::array<double, 4> distances = {{
        // Compare toMergerItr min with clusIndexStart min, if this is the min distance, lines are not colinear, merge is vetoed
          ::norm(x11-x21, y11-y21),
          ::norm(x11-x22, y11-y22), // Compare toMergerItr min with clusIndexStart max
          ::norm(x12-x21, y12-y21), // Compare toMergerItr max with clusIndexStart min
        // Compare toMergerItr max with clusIndexStart max, if this is the min distance, lines are not colinear, merge is vetoed
          ::norm(x12-x22, y12-y22)
        }}; // distances

        double minDistance = 999999; 
        int minDistanceIndex = -1;
        for(unsigned int j = 0; j < distances.size(); ++j){
          if (distances[j] < minDistance){
            minDistance = distances[j];
            minDistanceIndex = j;
          }
        }

        if(minDistanceIndex  == 0 || minDistanceIndex  == 3)
          continue;

          

        // Check if we merged left or right already for clusIndexStart, we only do it once for each side
        //if(showerClustersClusIndexStartProtoTrackItr->mergedLeft == true && minDistanceIndex == 2)
          //continue;
        //if(showerClustersClusIndexStartProtoTrackItr->mergedRight == true && minDistanceIndex == 1)
          //continue;
        //if(showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].mergedLeft == true && minDistanceIndex == 1)
          //continue;
        //if(showerClustersToMergeItr->clusterProtoTracks[*toMergeItr].mergedRight == true && minDistanceIndex == 2)
          //continue;

        //std::cout << "Potential merge" << std::endl;
        //std::cout << "main showerClustersItr slope: " << showerClusters->at(*toMergeItr).clusterSlope << " clusIndexStart slope: " << showerClusters->at(clusIndexStart).clusterSlope << std::endl;
        performedBestMerge=true;     


        if(closestDistance < bestToMergeShowerClusterProtoTrackDistance){
          bestToMergeShowerCluster=showerClustersToMergeItr-showerClusters->begin();
          bestToMergeShowerClusterProtoTrack=*toMergeItr;
          bestShowerClustersClusIndexStartProtoTrack=showerClustersClusIndexStartProtoTrackItr-showerClusters->at(clusIndexStart).clusterProtoTracks.begin();
          // Did we merge left (0) or right (1)?
          if(minDistanceIndex == 1){
            bestToMergeRightLeft = 0;
            //bestClusIndexStartRightLeft = 1;
          }
          if(minDistanceIndex == 2){
            bestToMergeRightLeft = 1;
            //bestClusIndexStartRightLeft = 0;
          }

        }

      }// End of loop over toMerge
    }// End of loop over showerClusters->begin()+clusIndexStart+1
  }//End of loop over showerClusters->at(clusIndexStart).clusterProtoTracks

  if(performedBestMerge){
    showerClusters->at(bestToMergeShowerCluster).clusterProtoTracks[bestToMergeShowerClusterProtoTrack].merged=true;
    showerClusters->at(clusIndexStart).clusterProtoTracks[bestShowerClustersClusIndexStartProtoTrack].merged=true;   
    if(bestToMergeRightLeft == 0){
      showerClusters->at(bestToMergeShowerCluster).clusterProtoTracks[bestToMergeShowerClusterProtoTrack].mergedLeft = true;
      showerClusters->at(clusIndexStart).clusterProtoTracks[bestShowerClustersClusIndexStartProtoTrack].mergedRight = true;   
    }
    if(bestToMergeRightLeft == 1){
      showerClusters->at(bestToMergeShowerCluster).clusterProtoTracks[bestToMergeShowerClusterProtoTrack].mergedRight = true;
      showerClusters->at(clusIndexStart).clusterProtoTracks[bestShowerClustersClusIndexStartProtoTrack].mergedLeft = true;   
    }
   
      showerClusters->at(clusIndexStart).addProtoTracks(showerClusters->at(bestToMergeShowerCluster).clusterProtoTracks);
      showerClusters->at(bestToMergeShowerCluster).clearProtoTracks();
      //std::cout << "Merged shower-shower" << std::endl;

  }

  return performedBestMerge;

}



//------------------------------------------------------------------------------
double cluster::fuzzyClusterAlg::DistanceBetweenHits(art::Ptr<recob::Hit> hit0,
                                                art::Ptr<recob::Hit> hit1,
                                                double wire_dist,
                                                double tickToDist)
{
  const double pHit0[2] = {
    (hit0->Channel())*wire_dist,
    hit0->PeakTime()*tickToDist
  };
  const double pHit1[2] = {
    (hit1->Channel())*wire_dist,
    hit1->PeakTime()*tickToDist
  };

  return ::norm( pHit0[0] - pHit1[0], pHit0[1] - pHit1[1]);

}



//------------------------------------------------------------------------------
double cluster::fuzzyClusterAlg::HoughLineDistance(double p0MinLine1, 
						double p1MinLine1, 
						double p0MaxLine1, 
						double p1MaxLine1, 
						double p0MinLine2, 
						double p1MinLine2, 
						double p0MaxLine2, 
						double p1MaxLine2)
{
  //distance between two segments in the plane:
  //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
  //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
  double x11 = p0MinLine1; 
  double y11 = p1MinLine1; 
  double x12 = p0MaxLine1; 
  double y12 = p1MaxLine1; 
  double x21 = p0MinLine2; 
  double y21 = p1MinLine2; 
  double x22 = p0MaxLine2; 
  double y22 = p1MaxLine2; 

  if(HoughLineIntersect(x11, y11, x12, y12, x21, y21, x22, y22)) return 0;
  // try each of the 4 vertices w/the other segment
  std::vector<double> distances;
  distances.push_back(PointSegmentDistance(x11, y11, x21, y21, x22, y22));
  distances.push_back(PointSegmentDistance(x12, y12, x21, y21, x22, y22));
  distances.push_back(PointSegmentDistance(x21, y21, x11, y11, x12, y12));
  distances.push_back(PointSegmentDistance(x22, y22, x11, y11, x12, y12));

  double minDistance = 999999; 
  for(unsigned int j = 0; j < distances.size(); j++){
    if (distances[j] < minDistance)
      minDistance = distances[j];
  }
  
  return minDistance;

}



//------------------------------------------------------------------------------
bool cluster::fuzzyClusterAlg::HoughLineIntersect(double x11,
					       double  y11,
					       double  x12,
					       double  y12,
					       double  x21,
					       double  y21,
					       double  x22,
					       double  y22)
{
  //whether two segments in the plane intersect:
  //one segment is (x11, y11) to (x12, y12)
  //the other is   (x21, y21) to (x22, y22)
  
  double dx1 = x12 - x11; // x2-x1
  double dy1 = y12 - y11; // y2-y1
  double dx2 = x22 - x21; // x4-x3
  double dy2 = y22 - y21; // y4-y3
  //double delta = dx2*dy1 - dy2*dx1; // (x4-x3)(y2-y1) - (y4-y3)(x2-x1)
  double delta = dy2*dx1 - dx2*dy1; // (y4-y3)(x2-x1) - (x4-x3)(y2-y1) 
  if (delta == 0) return false;  // parallel segments

  double t = (dx2*(y11 - y21) + dy2*(x21 - x11)) / delta; // ua
  double s = (dx1*(y11 - y21) + dy1*(x21 - x11)) / delta; // ub
  
  return (0 <= s && s <= 1 && 0 <= t && t <= 1);

}



//------------------------------------------------------------------------------
double cluster::fuzzyClusterAlg::PointSegmentDistance(double px,
						   double  py,
						   double  x1,
						   double  y1,
						   double  x2,
						   double  y2)
{
  double dx = x2 - x1;
  double dy = y2 - y1;
  if ( dx == 0 && dy == 0 )  // the segment's just a point
    return ::norm( px - x1, py - y1 );

  // Calculate the t that minimizes the distance.
  double t = ((px - x1)*dx + (py - y1)*dy) / ::sumsq(dx, dy);

  // See if this represents one of the segment's
  // end points or a point in the middle.
  if(t < 0){
    dx = px - x1;
    dy = py - y1;
  }
  else if(t > 1) {
    dx = px - x2;
    dy = py - y2;
  }
  else if(0 <= t && t <= 1) {
    double near_x = x1 + t * dx;
    double near_y = y1 + t * dy;
    dx = px - near_x;
    dy = py - near_y;
  }

  return ::norm(dx, dy);

}



