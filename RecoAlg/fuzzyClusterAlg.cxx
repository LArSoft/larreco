////////////////////////////////////////////////////////////////////////
//
// fuzzyClusterAlg.cxx
//
// Ben Carls, bcarls@fnal.gov
//
// This is an adpatation of the FLAME clustering algorithm, see the 
// notes below
//
/*===================================================================
  The standard implementation of FLAME data clustering algorithm.

  FLAME (FLAME clustering by Local Approximation of MEmberships)
  was first described in:
  "FLAME, a novel FLAME clustering method for the analysis of DNA
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
//
////////////////////////////////////////////////////////////////////////


#include <boost/bind.hpp>

//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//#include <TStopwatch.h>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "CLHEP/Random/RandFlat.h"
#include "Filters/ChannelFilter.h"

#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/fuzzyClusterAlg.h"
#include "RecoBase/Hit.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/AssociationUtil.h"

#include <time.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

// Define parameters that will tell us if we are doing a normal Hough line merge
// or a shower Hough line merge
static const int iMergeShower          = 0;
static const int iMergeNormal          = 1;
static const int iMergeShowerIntercept = 2;
static const int iMergeChargeAsymAngle = 3;


namespace cluster{
  const unsigned int kNO_CLUSTER    = UINT_MAX;
  const unsigned int kNOISE_CLUSTER = UINT_MAX-1;
}

//----------------------------------------------------------
// fuzzyClusterAlg stuff
//----------------------------------------------------------
cluster::fuzzyClusterAlg::fuzzyClusterAlg(fhicl::ParameterSet const& pset) 
   : fHBAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
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
  fIterations                     = p.get< int    >("Iterations");
  fEpsilon                        = p.get< double >("Epsilon");
  fKNN                            = p.get< double >("KNN");
  fDistanceMetric                 = p.get< int    >("DistanceMetric");
  fThreshold                      = p.get< double >("Threshold");

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
  fMaxVertexLines                 = p.get< int   >("MaxVertexLines"                 );
  fVertexLinesCutoff              = p.get< double >("VertexLinesCutoff"              );
  fHBAlg.reconfigure(p.get< fhicl::ParameterSet >("HoughBaseAlg"));
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

  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<util::DetectorProperties> detp;

  fWirePitch.push_back(fGeom->WirePitch(0,1,0));
  fWirePitch.push_back(fGeom->WirePitch(0,1,1));
  fWirePitch.push_back(fGeom->WirePitch(0,1,2));

  // Collect the hits in a useful form,
  // and take note of the maximum time width
  fMaxWidth=0.0;
  //fpsMat = TMatrixT<double>(allhits.size(),2);
  //fpsMembership = TMatrixT<double>(iNumClusters, allhits.size());
  fpsMat.ResizeTo(allhits.size(),2);
  
  NNumOfRows = allhits.size();
  MNumOfCols = 2;
  data = (double**) malloc(NNumOfRows * sizeof(double*) );


  double tickToDist = larp->DriftVelocity(larp->Efield(),larp->Temperature());
  tickToDist *= 1.e-3 * detp->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
  int dims = 3;//our point is defined by 3 elements:wire#,center of the hit, and the hit width
  std::vector<double> p(dims);
  for (auto allhitsItr = allhits.begin(); allhitsItr < allhits.end(); allhitsItr++){
        
    //    p[0] = ((*allhitsItr)->Channel())*fGeom->WirePitch(fGeom->View((*allhitsItr)->Channel()));
    //    std::cout << "plane " << (*allhitsItr)->WireID().Plane << "  pitch " << fWirePitch[(*allhitsItr)->WireID().Plane]
    //	      << std::endl;
    p[0] = ((*allhitsItr)->WireID().Wire)*fWirePitch[(*allhitsItr)->WireID().Plane];
    p[1] = (((*allhitsItr)->StartTime()+(*allhitsItr)->EndTime()  )/2.)*tickToDist;
    p[2] = (*allhitsItr)->Charge();   //width of a hit in cm

    // check on the maximum width condition
    if ( p[2] > fMaxWidth ) fMaxWidth = p[2];
    
    fps.push_back(p);

    // Store hits in the matrix needed for FLAME clustering
    fpsMat(allhitsItr-allhits.begin(),0) = p[0];
    fpsMat(allhitsItr-allhits.begin(),1) = p[1];

    data[allhitsItr-allhits.begin()] = (double*) malloc( MNumOfCols * sizeof(double*) );
    //data[allhitsItr-allhits.begin()][0] = ((*allhitsItr)->Channel())*fGeom->WirePitch(fGeom->View((*allhitsItr)->Channel()));
    data[allhitsItr-allhits.begin()][0] = ((*allhitsItr)->WireID().Wire)*fWirePitch[(*allhitsItr)->WireID().Plane];
    data[allhitsItr-allhits.begin()][1] = (((*allhitsItr)->StartTime()+(*allhitsItr)->EndTime()  )/2.)*tickToDist;
    //data[allhitsItr-allhits.begin()][3] = (*allhitsItr)->Charge();
  }

  mf::LogInfo("fuzzyCluster") << "Initfuzzy: hits vector size is " << fps.size();

  return;
}


//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
//
//  Ben Carls' implementation of fuzzyClusterAlg as much like examples as possible
void cluster::fuzzyClusterAlg::run_fuzzy_cluster(std::vector<art::Ptr<recob::Hit> >& allhits) {

  // Don't attempt to run the algorithm if we have 1 or fewer hits
  if(allhits.size() <= 1)
    return;

  Flame *flame;
  flame = Flame_New();

  // Setup data matrix
  Flame_SetDataMatrix(flame, data, NNumOfRows, MNumOfCols, fDistanceMetric);

  // Detecting Cluster Supporting Objects ...
  //Flame_DefineSupports( flame, 10, -2.0 );
  Flame_DefineSupports( flame, fKNN, fThreshold);

  // Propagating FLAME memberships ... 
  Flame_LocalApproximation( flame, fIterations, fEpsilon);

  // Defining clusters from FLAME memberships ... 
  Flame_MakeClusters( flame, -1.0 );



  fpointId_to_clusterId.resize(fps.size(), kNO_CLUSTER); // Not zero as before!
  fnoise.resize(fps.size(), false);
  fvisited.resize(fps.size(), false);

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;

  if(allhits.size()==0)
    return;

  //factor to make x and y scale the same units
//   uint32_t     channel = allhits[0]->Channel();
//   double wirePitch = geom->WirePitch(geom->View(channel));
    uint32_t     channel = allhits[0]->WireID().Wire;
//     double wirePitch[2];
//     wirePitch[0] = geom->WirePitch(0);
//     wirePitch[1] = geom->WirePitch(1);
//     wirePitch[2] = geom->WirePitch(2);

// saw this one

    double xyScale = .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
    xyScale*=detprop->SamplingRate();
    //  double wire_dist = wirePitch;

    double tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
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
    

    unsigned int cid = flame->cso_count+1;
    
  // Loop over clusters 
    for (int i = 0; i<= flame->cso_count; i++){
      // Loop over hits in cluster
      for (int j=0; j < flame->clusters[i].size; j++){
	if(i == flame->cso_count)
	  fpointId_to_clusterId[flame->clusters[i].array[j]] = kNOISE_CLUSTER;
	else
	  fpointId_to_clusterId[flame->clusters[i].array[j]] = i;
      }
    }
    
    unsigned int nClusters = cid;    
    unsigned int nClustersTemp = cid;
  
  // Loop over clusters with the Hough line finder to break the clusters up further
  // list of lines
  std::vector<protoTrack> protoTracksFound;
  if(nClustersTemp > 0)
    for (unsigned int i = 0; i <= (unsigned int)nClustersTemp-1; ++i){
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
      distance = (TMath::Abs((*hitsItr)->PeakTime()-protoTracksFoundItr->clusterSlope*(double)((*hitsItr)->WireID().Wire)-protoTracksFoundItr->clusterIntercept)/(std::sqrt(pow(xyScale/fWirePitch[(*hitsItr)->WireID().Plane]*protoTracksFoundItr->clusterSlope,2)+1)));
      /// Sum up background hits, use smart distance
      peakTimePerpMin=-(1/protoTracksFoundItr->clusterSlope)*(double)((*hitsItr)->WireID().Wire)+allhits[protoTracksFoundItr->iMinWire]->PeakTime()+(1/protoTracksFoundItr->clusterSlope)*(allhits[protoTracksFoundItr->iMinWire]->WireID().Wire);
      peakTimePerpMax=-(1/protoTracksFoundItr->clusterSlope)*(double)((*hitsItr)->WireID().Wire)+allhits[protoTracksFoundItr->iMaxWire]->PeakTime()+(1/protoTracksFoundItr->clusterSlope)*(allhits[protoTracksFoundItr->iMaxWire]->WireID().Wire);
      if(distance > 1*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling)
         && distance < 25*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling)){
        if((protoTracksFoundItr->clusterSlope < 0 && (*hitsItr)->PeakTime() < peakTimePerpMin && (*hitsItr)->PeakTime() > peakTimePerpMax)
            || (protoTracksFoundItr->clusterSlope > 0 && (*hitsItr)->PeakTime() > peakTimePerpMin && (*hitsItr)->PeakTime() < peakTimePerpMax)){
          totalBkgDistCharge+=distance/(*hitsItr)->Charge();
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
    if(!protoTracksFoundSizes.count(protoTracksFoundItr->clusterNumber))
      protoTracksFoundSizes[protoTracksFoundItr->clusterNumber] = std::sqrt( pow(protoTracksFoundItr->pMin0-protoTracksFoundItr->pMax0,2)+pow(protoTracksFoundItr->pMin1-protoTracksFoundItr->pMax1,2));
    else 
      protoTracksFoundSizes[protoTracksFoundItr->clusterNumber]+= std::sqrt( pow(protoTracksFoundItr->pMin0-protoTracksFoundItr->pMax0,2)+pow(protoTracksFoundItr->pMin1-protoTracksFoundItr->pMax1,2));
  }
 



  std::vector< art::Ptr<recob::Hit> > unclusteredhits;
  std::vector<unsigned int> unclusteredhitsToallhits;
  bool unclustered;
  double p0;
  double p1;
  double minDistanceShower;
  double minDistanceTrack;
  if(fDoFuzzyRemnantMerge){
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
      p0 = ((*allhitsItr)->Channel())*fWirePitch[ip];
      p1 = (((*allhitsItr)->StartTime()+(*allhitsItr)->EndTime())/2.)*tickToDist;

      // First try to group it with a shower-like cluster
      minDistanceShower = 999999;
      for(auto showerClustersItr = showerClusters.begin(); showerClustersItr != showerClusters.end(); ++showerClustersItr){
        for(auto protoTracksItr = showerClustersItr->clusterProtoTracks.begin(); protoTracksItr < showerClustersItr->clusterProtoTracks.end(); ++protoTracksItr){
          
          distance = PointSegmentDistance( p0, p1, protoTracksItr->pMin0, protoTracksItr->pMin1, protoTracksItr->pMax0, protoTracksItr->pMax1);

          if(distance > fFuzzyRemnantMergeCutoff)
            continue;

          distance/=pow(protoTracksFoundSizes[protoTracksItr->clusterNumber],1/4);
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
      minDistanceTrack = 999999;
      for(auto trackClustersItr = trackClusters.begin(); trackClustersItr != trackClusters.end(); ++trackClustersItr){
        for(auto protoTracksItr = trackClustersItr->clusterProtoTracks.begin(); protoTracksItr < trackClustersItr->clusterProtoTracks.end(); ++protoTracksItr){
          
          distance = PointSegmentDistance( p0, p1, protoTracksItr->pMin0, protoTracksItr->pMin1, protoTracksItr->pMax0, protoTracksItr->pMax1);

          if(distance > fFuzzyRemnantMergeCutoff)
            continue;

          //distance/=pow(protoTracksFoundSizes[protoTracksItr->clusterNumber],1/4);
          if(distance < minDistanceTrack){
            fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = protoTracksItr->clusterNumber;
            minDistanceTrack = distance;
            unclustered = false;
          }
        }
      }



      if(unclustered){
        unclusteredhitsToallhits.push_back(allhitsItr-allhits.begin());
        unclusteredhits.push_back(*allhitsItr);
      }
      
    }

  }













  cid = nClusters;
  
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
          double mergeRightClusIndexStartDist = std::sqrt(pow(x11-x22,2) + pow(y11-y22,2));
          // Compare toMergerItr max with clusIndexStart min
          double mergeLeftClusIndexStartDist = std::sqrt(pow(x12-x21,2) + pow(y12-y21,2));
         
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
        mergeTheta[mergeThetaItr-mergeTheta.begin()] = atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
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
            //double distance = std::sqrt(pow(clusIndStHitItr->first-(*toMergeHitItr).first,2)+
                      //pow(clusIndStHitItr->second-toMergeHitItr->second,2));
            
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
        std::vector<double> distances;

        // Compare toMergerItr min with clusIndexStart min, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt((x11-x21)*(x11-x21)+ (y11-y21)*(y11-y21)));
        // Compare toMergerItr min with clusIndexStart max
        distances.push_back(std::sqrt((x11-x22)*(x11-x22) + (y11-y22)*(y11-y22)));
        // Compare toMergerItr max with clusIndexStart min
        distances.push_back(std::sqrt((x12-x21)*(x12-x21) + (y12-y21)*(y12-y21)));
        // Compare toMergerItr max with clusIndexStart max, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt((x12-x22)*(x12-x22) + (y12-y22)*(y12-y22)));

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
          double mergeRightClusIndexStartDist = std::sqrt(pow(x11-x22,2) + pow(y11-y22,2));
          // Compare toMergerItr max with clusIndexStart min
          double mergeLeftClusIndexStartDist = std::sqrt(pow(x12-x21,2) + pow(y12-y21,2));
          
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
        mergeTheta[mergeThetaItr-mergeTheta.begin()] = atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
      }


      // Perform the merge
      for(auto toMergeItr = toMerge.begin(); toMergeItr != toMerge.end(); toMergeItr++){

        // First check averages of charge and sigma charge for hits in lines closest to each other
        int closestToMerge=-1;
        int closestClusIndexStart=-1;
        double closestDistance=999999;
        for (auto toMergeHitItr = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin(); toMergeHitItr != trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.end(); toMergeHitItr++) {
          for (auto clusIndStHitItr = trackClustersClusIndexStartProtoTrackItr->hits.begin(); clusIndStHitItr != trackClustersClusIndexStartProtoTrackItr->hits.end(); clusIndStHitItr++) {
            //double distance = std::sqrt(pow(clusIndStHitItr->first-(*toMergeHitItr).first,2)+
                      //pow(clusIndStHitItr->second-toMergeHitItr->second,2));
            
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
        std::vector<std::pair<int,double> > closestToMergeDist;
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
              closestToMergeDist.push_back(std::make_pair(toMergeHitItr-trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.begin(),distance));
              std::sort(closestToMergeDist.begin(), closestToMergeDist.end(), boost::bind(&std::pair<int,double>::second,_1) < boost::bind(&std::pair<int,double>::second,_2));
            }
            if(closestToMergeDist.size() > trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits.size()-1 ||
              closestToMergeDist.size() > 9)
              closestToMergeDist.erase(closestToMergeDist.end());
        }
        //for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end();
          //closestToMergeDistItr++) 
          //std::cout << closestToMergeDistItr->first << " " << closestToMergeDistItr->second << std::endl;



        // Find up to 9 more points closest to closestToMerge on the clusIndexStart line
        std::vector<std::pair<int,double> > closestClusIndexStartDist;
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
            closestClusIndexStartDist.push_back(std::make_pair(clusIndexStartHitItr-trackClustersClusIndexStartProtoTrackItr->hits.begin(),distance));
            std::sort(closestClusIndexStartDist.begin(), closestClusIndexStartDist.end(), boost::bind(&std::pair<int,double>::second,_1) < boost::bind(&std::pair<int,double>::second,_2));
          }
          if(closestClusIndexStartDist.size() > trackClustersClusIndexStartProtoTrackItr->hits.size()-1 ||
            closestClusIndexStartDist.size() > 9)
            closestClusIndexStartDist.erase(closestClusIndexStartDist.end());
       }



        double toMergeAveCharge = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMerge]->Charge();
        double toMergeAveSigmaCharge = trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMerge]->SigmaCharge();
        for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end(); ++closestToMergeDistItr) {
          toMergeAveCharge+= trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMergeDistItr->first]->Charge();
          toMergeAveSigmaCharge+= trackClustersToMergeItr->clusterProtoTracks[*toMergeItr].hits[closestToMergeDistItr->first]->SigmaCharge();
        }
        double clusIndexStartAveCharge = trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStart]->Charge();
        double clusIndexStartAveSigmaCharge = trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStart]->SigmaCharge();
        for(auto closestClusIndexStartDistItr = closestClusIndexStartDist.begin(); closestClusIndexStartDistItr != closestClusIndexStartDist.end(); ++closestClusIndexStartDistItr) {
          clusIndexStartAveCharge+= trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStartDistItr->first]->Charge();
          clusIndexStartAveSigmaCharge+=trackClustersClusIndexStartProtoTrackItr->hits[closestClusIndexStartDistItr->first]->SigmaCharge();
        }



        double chargeAsymmetry = std::abs(toMergeAveCharge-clusIndexStartAveCharge)/(toMergeAveCharge+clusIndexStartAveCharge);
        double sigmaChargeAsymmetry = std::abs(toMergeAveSigmaCharge-clusIndexStartAveSigmaCharge)/(toMergeAveSigmaCharge+clusIndexStartAveSigmaCharge);
        double chargeAsymmetrySinAngle = chargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1);
        double sigmaChargeAsymmetrySinAngle = sigmaChargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1);

        //std::cout << std::endl;
        //std::cout << chargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1) << std::endl;
        //std::cout << sigmaChargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1) << std::endl;
        

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
        std::vector<double> distances;

        // Compare toMergerItr min with clusIndexStart min, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt(pow(x11-x21,2) + pow(y11-y21,2)));
        // Compare toMergerItr min with clusIndexStart max
        distances.push_back(std::sqrt(pow(x11-x22,2) + pow(y11-y22,2)));
        // Compare toMergerItr max with clusIndexStart min
        distances.push_back(std::sqrt(pow(x12-x21,2) + pow(y12-y21,2)));
        // Compare toMergerItr max with clusIndexStart max, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt(pow(x12-x22,2) + pow(y12-y22,2)));

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
        mergeTheta[mergeThetaItr-mergeTheta.begin()] = atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
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
            //double distance = std::sqrt(pow(clusIndStHitItr->first-(*toMergeHitItr).first,2)+
                      //pow(clusIndStHitItr->second-toMergeHitItr->second,2));
            
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
        std::vector<double> distances;

        // Compare toMergerItr min with clusIndexStart min, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt(pow(x11-x21,2) + pow(y11-y21,2)));
        // Compare toMergerItr min with clusIndexStart max
        distances.push_back(std::sqrt(pow(x11-x22,2) + pow(y11-y22,2)));
        // Compare toMergerItr max with clusIndexStart min
        distances.push_back(std::sqrt(pow(x12-x21,2) + pow(y12-y21,2)));
        // Compare toMergerItr max with clusIndexStart max, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt(pow(x12-x22,2) + pow(y12-y22,2)));

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




























// Merges based on the distance between line segments
void cluster::fuzzyClusterAlg::mergeHoughLinesBySegment(unsigned int clusIndexStart,
						     std::vector<protoTrack> *tracksFound,
						     double xyScale,
                                                     int mergeStyle,
                                                     double wire_dist,
                                                     double tickToDist)
{



  // If we have zero or one Hough lines, move on 
  if(tracksFound->size() == 0 || tracksFound->size() == 1)
    return;

  // If a merge happened, trigger this function again to look at this cluster again!
  bool lineMerged = false;

  // If we reach the last Hough line, move on 
  if(tracksFound->size() == clusIndexStart+1)
    return;

  // Min to merge
  std::vector<unsigned int> toMerge; 
  std::vector<double> mergeSlope;
  std::vector<double> mergeTheta;

  // Check if segments are close enough
  for(auto tracksFoundToMergeItr = tracksFound->begin(); tracksFoundToMergeItr != tracksFound->end(); ++tracksFoundToMergeItr){
    if(tracksFound->at(clusIndexStart).clusterNumber == tracksFoundToMergeItr->clusterNumber)
      continue;
    double segmentDistance = HoughLineDistance(tracksFound->at(clusIndexStart).pMin0,tracksFound->at(clusIndexStart).pMin1,
                                               tracksFound->at(clusIndexStart).pMax0,tracksFound->at(clusIndexStart).pMax1, 
      					       tracksFoundToMergeItr->pMin0,tracksFoundToMergeItr->pMin1,
                                               tracksFoundToMergeItr->pMax0,tracksFoundToMergeItr->pMax1);
    if( (segmentDistance<fTrackClusterMergeCutoff && mergeStyle == iMergeNormal) 
        || (segmentDistance<fShowerClusterMergeCutoff && (mergeStyle == iMergeShower || mergeStyle == iMergeShowerIntercept))
        || (segmentDistance<fTrackClusterMergeCutoff && mergeStyle == iMergeChargeAsymAngle))
    {
      //std::cout << std::endl;
      //std::cout << tracksFoundClusIndStItr->minWire << " " << tracksFoundClusIndStItr->maxWire << std::endl;
      //std::cout << tracksFoundToMergeItr->minWire << " " << tracksFoundToMergeItr->maxWire << std::endl;
      //std::cout << segmentDistance << std::endl;
      toMerge.push_back(tracksFoundToMergeItr-tracksFound->begin());
      mergeSlope.push_back(tracksFound->at(clusIndexStart).clusterSlope*xyScale);
    }
  }

  mergeTheta.resize(toMerge.size());

  // Find the angle between the slopes
  for(auto mergeThetaItr = mergeTheta.begin(); mergeThetaItr != mergeTheta.end(); ++mergeThetaItr){
    double toMergeSlope = tracksFound->at(toMerge[mergeThetaItr-mergeTheta.begin()]).clusterSlope*xyScale;
    mergeTheta[mergeThetaItr-mergeTheta.begin()] = atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
    //std::cout << std::endl;
    //std::cout << "toMergeSlope: " << toMergeSlope/xyScale<< " mergeSlope[clusIndexStart]: " << mergeSlope[mergeThetaItr-mergeTheta.begin()]/tickToDist << std::endl;
    //std::cout << "mergeTheta: " << mergeTheta[mergeThetaItr-mergeTheta.begin()] << std::endl;
  }

  // Perform the merge
  for(auto toMergeItr = toMerge.begin(); toMergeItr != toMerge.end(); ++toMergeItr){
    if(  
        (mergeTheta[toMergeItr-toMerge.begin()] < fShowerClusterMergeAngle && mergeStyle == iMergeShower)
        || mergeStyle == iMergeShowerIntercept
        || mergeStyle == iMergeChargeAsymAngle){
     

      // First check averages of charge and sigma charge for hits in lines closest to each other
      int closestToMerge=-1;
      int closestClusIndexStart=-1;
      double closestDistance=999999;
      for (auto toMergeHitItr = tracksFound->at(*toMergeItr).hits.begin(); toMergeHitItr != tracksFound->at(*toMergeItr).hits.end(); ++toMergeHitItr) {
        for (auto clusIndStHitItr = tracksFound->at(clusIndexStart).hits.begin(); clusIndStHitItr != tracksFound->at(clusIndexStart).hits.end(); ++clusIndStHitItr) {
          //double distance = std::sqrt(pow(clusIndStHitItr->first-(*toMergeHitItr).first,2)+
                    //pow(clusIndStHitItr->second-toMergeHitItr->second,2));
          double distance = DistanceBetweenHits(*clusIndStHitItr,
                                                *toMergeHitItr,
                                                wire_dist,
                                                tickToDist);
          if(distance < closestDistance){
            closestDistance = distance;
            closestToMerge=toMergeHitItr-tracksFound->at(*toMergeItr).hits.begin();
            closestClusIndexStart=clusIndStHitItr-tracksFound->at(clusIndexStart).hits.begin();
          }
        }
      }

      // Find up to 9 more points closest to closestToMerge on the toMerge[i] line
      // check if it's closer, insert, delete
      std::vector<std::pair<int,double> > closestToMergeDist;
      for (auto toMergeHitItr = tracksFound->at(*toMergeItr).hits.begin(); toMergeHitItr != tracksFound->at(*toMergeItr).hits.end(); ++toMergeHitItr) {
        if(closestToMerge==toMergeHitItr-tracksFound->at(*toMergeItr).hits.begin())
          continue;
          double distance = DistanceBetweenHits(tracksFound->at(clusIndexStart).hits[closestClusIndexStart],
                                                *toMergeHitItr,
                                                wire_dist,
                                                tickToDist);

        bool foundCloser = false;
        for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end(); ++closestToMergeDistItr) {
          if(closestToMergeDistItr->second > distance){
            foundCloser = true;
            break;
          }
        }
        if(foundCloser 
            || closestToMergeDist.size() < tracksFound->at(*toMergeItr).hits.size()-1
            || closestToMergeDist.size() < 9){
          closestToMergeDist.push_back(std::make_pair(toMergeHitItr-tracksFound->at(*toMergeItr).hits.begin(),distance));
          std::sort(closestToMergeDist.begin(), closestToMergeDist.end(), boost::bind(&std::pair<int,double>::second,_1) < boost::bind(&std::pair<int,double>::second,_2));
        }
        if(closestToMergeDist.size() > tracksFound->at(*toMergeItr).hits.size()-1 ||
          closestToMergeDist.size() > 9)
          closestToMergeDist.erase(closestToMergeDist.end());
      }
      //for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end();
        //closestToMergeDistItr++) 
        //std::cout << closestToMergeDistItr->first << " " << closestToMergeDistItr->second << std::endl;



      // Find up to 9 more points closest to closestToMerge on the clusIndexStart line
      std::vector<std::pair<int,double> > closestClusIndexStartDist;
      for (auto clusIndexStartHitItr = tracksFound->at(clusIndexStart).hits.begin(); clusIndexStartHitItr != tracksFound->at(clusIndexStart).hits.end(); ++clusIndexStartHitItr) {
        if(closestClusIndexStart==clusIndexStartHitItr-tracksFound->at(clusIndexStart).hits.begin())
          continue;

        double distance = DistanceBetweenHits(*clusIndexStartHitItr,
                                              tracksFound->at(*toMergeItr).hits[closestToMerge],
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
            || closestClusIndexStartDist.size() < tracksFound->at(clusIndexStart).hits.size()-1
            || closestClusIndexStartDist.size() < 9){
          closestClusIndexStartDist.push_back(std::make_pair(clusIndexStartHitItr-tracksFound->at(clusIndexStart).hits.begin(),distance));
          std::sort(closestClusIndexStartDist.begin(), closestClusIndexStartDist.end(), boost::bind(&std::pair<int,double>::second,_1) < boost::bind(&std::pair<int,double>::second,_2));
        }
        if(closestClusIndexStartDist.size() > tracksFound->at(clusIndexStart).hits.size()-1 ||
          closestClusIndexStartDist.size() > 9)
          closestClusIndexStartDist.erase(closestClusIndexStartDist.end());
     }



      double toMergeAveCharge = tracksFound->at(*toMergeItr).hits[closestToMerge]->Charge();
      double toMergeAveSigmaCharge = tracksFound->at(*toMergeItr).hits[closestToMerge]->SigmaCharge();
      for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end(); ++closestToMergeDistItr) {
        toMergeAveCharge+=tracksFound->at(*toMergeItr).hits[closestToMergeDistItr->first]->Charge();
        toMergeAveSigmaCharge+=tracksFound->at(*toMergeItr).hits[closestToMergeDistItr->first]->SigmaCharge();
      }
      double clusIndexStartAveCharge = tracksFound->at(clusIndexStart).hits[closestClusIndexStart]->Charge();
      double clusIndexStartAveSigmaCharge = tracksFound->at(clusIndexStart).hits[closestClusIndexStart]->SigmaCharge();
      for(auto closestClusIndexStartDistItr = closestClusIndexStartDist.begin(); closestClusIndexStartDistItr != closestClusIndexStartDist.end(); ++closestClusIndexStartDistItr) {
        clusIndexStartAveCharge+=tracksFound->at(clusIndexStart).hits[closestClusIndexStartDistItr->first]->Charge();
        clusIndexStartAveSigmaCharge+=tracksFound->at(clusIndexStart).hits[closestClusIndexStartDistItr->first]->SigmaCharge();
      }



      double chargeAsymmetry = std::abs(toMergeAveCharge-clusIndexStartAveCharge)/(toMergeAveCharge+clusIndexStartAveCharge);
      double sigmaChargeAsymmetry = std::abs(toMergeAveSigmaCharge-clusIndexStartAveSigmaCharge)/(toMergeAveSigmaCharge+clusIndexStartAveSigmaCharge);
      double chargeAsymmetrySinAngle = chargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1);
      double sigmaChargeAsymmetrySinAngle = sigmaChargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1);

      //std::cout << std::endl;
      //std::cout << chargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1) << std::endl;
      //std::cout << sigmaChargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1) << std::endl;
      

      if(chargeAsymmetrySinAngle > fChargeAsymAngleCut &&
          mergeStyle == iMergeChargeAsymAngle)
        continue;

      if(sigmaChargeAsymmetrySinAngle > fSigmaChargeAsymAngleCut  &&
          mergeStyle == iMergeChargeAsymAngle)
        continue;



















      //double lineLengthAsymm = std::abs((double)tracksFound->at(clusIndexStart).hits.size() - (double)tracksFound->at(*toMergeItr).hits.size())/((double)tracksFound->at(clusIndexStart).hits.size() + (double)tracksFound->at(*toMergeItr).hits.size());

      //std::cout << "Length Asymm: " << lineLengthAsymm << std::endl;
      //std::cout << tracksFound->at(clusIndexStart).hits.size() << std::endl;

      //if(lineLengthAsymm > 0.75) 
        //continue; 



















      // Veto the merge if the lines are not colinear 
      if(mergeStyle == iMergeNormal || mergeStyle == iMergeChargeAsymAngle) {
        // Find where the lines are closest
  
        //distance between two segments in the plane:
        //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
        //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
        double x11 = tracksFound->at(*toMergeItr).pMin0; 
        double y11 = tracksFound->at(*toMergeItr).pMin1; 
        double x12 = tracksFound->at(*toMergeItr).pMax0; 
        double y12 = tracksFound->at(*toMergeItr).pMax1; 
        double x21 = tracksFound->at(clusIndexStart).pMin0; 
        double y21 = tracksFound->at(clusIndexStart).pMin1; 
        double x22 = tracksFound->at(clusIndexStart).pMax0; 
        double y22 = tracksFound->at(clusIndexStart).pMax1; 
        std::vector<double> distances;

        // Compare toMergerItr min with clusIndexStart min, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt(pow(x11-x21,2) + pow(y11-y21,2)));
        // Compare toMergerItr min with clusIndexStart max
        distances.push_back(std::sqrt(pow(x11-x22,2) + pow(y11-y22,2)));
        // Compare toMergerItr max with clusIndexStart min
        distances.push_back(std::sqrt(pow(x12-x21,2) + pow(y12-y21,2)));
        // Compare toMergerItr max with clusIndexStart max, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt(pow(x12-x22,2) + pow(y12-y22,2)));

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

      }
      //std::cout << tracksFound->at(*toMergeItr).clusterSlope << " " << tracksFound->at(clusIndexStart).clusterSlope << std::endl;



      // Check if both lines is in region that looks showerlike
      // Or merge if the distance between the lines is zero and one looks showerlike
      //double segmentDistance = HoughLineDistance(tracksFound->at(clusIndexStart).pMin0,tracksFound->at(clusIndexStart).pMin1,
                                                 //tracksFound->at(clusIndexStart).pMax0,tracksFound->at(clusIndexStart).pMax1, 
                                                 //tracksFound->at(*toMergeItr).pMin0,tracksFound->at(*toMergeItr).pMin1,
                                                 //tracksFound->at(*toMergeItr).pMax0,tracksFound->at(*toMergeItr).pMax1);
      //std::cout << "segmentDistance: " << segmentDistance << std::endl;
      //std::cout << tracksFound->at(*toMergeItr).showerLikeness << " " << tracksFound->at(clusIndexStart).showerLikeness << std::endl;

      
     
      // If doing a shower merge, only merge if the Hough lines look showerlike
      if(mergeStyle == iMergeShower){
        if(!(tracksFound->at(*toMergeItr).showerLikeness>fShowerLikenessCut) || !(tracksFound->at(clusIndexStart).showerLikeness>fShowerLikenessCut))
          continue;
      }
      
      
      //// If not doing a shower merge, only merge if one of the Hough lines doesn't look showerlike
      //if(mergeStyle == iMergeChargeAsymAngle){
        //if(tracksFound->at(*toMergeItr).showerLikeness>fShowerLikenessCut && tracksFound->at(clusIndexStart).showerLikeness>fShowerLikenessCut)
          //continue;
      //}

      
      if(mergeStyle == iMergeShowerIntercept){
        if((tracksFound->at(*toMergeItr).showerLikeness<fShowerLikenessCut) && (tracksFound->at(clusIndexStart).showerLikeness<fShowerLikenessCut))
          continue;
        double x11 = tracksFound->at(*toMergeItr).pMin0; 
        double y11 = tracksFound->at(*toMergeItr).pMin1; 
        double x12 = tracksFound->at(*toMergeItr).pMax0; 
        double y12 = tracksFound->at(*toMergeItr).pMax1; 
        double x21 = tracksFound->at(clusIndexStart).pMin0; 
        double y21 = tracksFound->at(clusIndexStart).pMin1; 
        double x22 = tracksFound->at(clusIndexStart).pMax0; 
        double y22 = tracksFound->at(clusIndexStart).pMax1; 
        if(HoughLineIntersect(x11, y11, x12, y12, x21, y21, x22, y22) == 0)
          continue;
      }

      //std::cout << "Merging" << std::endl;
      lineMerged = true;
      tracksFound->at(clusIndexStart).merged = true;
      tracksFound->at(*toMergeItr).merged = true;

      // For loop over all lines found to reassign lines to clusIndexStart that already belonged to toMerge 
      for(auto tracksFoundItr = tracksFound->begin(); tracksFoundItr != tracksFound->end(); ++tracksFoundItr){
        if((unsigned int)(*toMergeItr) == tracksFoundItr-tracksFound->begin())
          continue;

        if(tracksFoundItr->clusterNumber == tracksFound->at(*toMergeItr).clusterNumber){
          tracksFoundItr->clusterNumber = tracksFound->at(clusIndexStart).clusterNumber;
          //std::cout << "tracksFoundItr slope: " << tracksFoundItr->clusterSlope << " clusIndexStart slope: " << tracksFound->at(clusIndexStart).clusterSlope << std::endl;

        }
      }
      tracksFound->at(*toMergeItr).clusterNumber = tracksFound->at(clusIndexStart).clusterNumber;
      //std::cout << "main tracksFoundItr slope: " << tracksFound->at(*toMergeItr).clusterSlope << " clusIndexStart slope: " << tracksFound->at(clusIndexStart).clusterSlope << std::endl;
     
    }  
  }

                                                
  if(lineMerged)
    mergeHoughLinesBySegment(clusIndexStart,tracksFound,xyScale,mergeStyle,wire_dist,tickToDist);
  else
    mergeHoughLinesBySegment(clusIndexStart+1,tracksFound,xyScale,mergeStyle,wire_dist,tickToDist);
  
  return;

}


//------------------------------------------------------------------------------
double cluster::fuzzyClusterAlg::DistanceBetweenHits(art::Ptr<recob::Hit> hit0,
                                                art::Ptr<recob::Hit> hit1,
                                                double wire_dist,
                                                double tickToDist)
{
  double pHit0[2];
  pHit0[0] = (hit0->Channel())*wire_dist;
  pHit0[1] = ((hit0->StartTime()+hit0->EndTime())/2.)*tickToDist;
  double pHit1[2];
  pHit1[0] = (hit1->Channel())*wire_dist;
  pHit1[1] = ((hit1->StartTime()+hit1->EndTime())/2.)*tickToDist;

  return std::sqrt( pow(pHit0[0] - pHit1[0],2) + pow(pHit0[1] - pHit1[1],2));

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
    return std::sqrt( pow(px - x1,2) + pow(py - y1,2));

  // Calculate the t that minimizes the distance.
  double t = ((px - x1)*dx + (py - y1)*dy) / (dx*dx + dy*dy);

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

  return std::sqrt(dx*dx + dy*dy);

}





/* Quick Sort.
 * Adam Drozdek: Data Structures and Algorithms in C++, 2nd Edition.
 */
void cluster::fuzzyClusterAlg::PartialQuickSort( Indexdouble *data, int first, int last, int part )
{
	int lower=first+1, upper=last;
	double pivot;
	Indexdouble val;
	if( first >= last ) return;
	val = data[first];
	data[first] = data[ (first+last)/2 ];
	data[ (first+last)/2 ] = val;
	pivot = data[ first ].value;

	while( lower <= upper ){
		while( lower <= last && data[lower].value < pivot ) lower ++;
		while( pivot < data[upper].value ) upper --;
		if( lower < upper ){
			val = data[lower];
			data[lower] = data[upper];
			data[upper] = val;
			upper --;
		}
		lower ++;
	}
	val = data[first];
	data[first] = data[upper];
	data[upper] = val;
	if( first < upper-1 ) PartialQuickSort( data, first, upper-1, part );
	if( upper >= part ) return;
	if( upper+1 < last ) PartialQuickSort( data, upper+1, last, part );
}

DistFunction basicDistFunctions[] =
{
        cluster::fuzzyClusterAlg::Flame_Euclidean ,
        cluster::fuzzyClusterAlg::Flame_CosineDist ,
        cluster::fuzzyClusterAlg::Flame_PearsonDist ,
        cluster::fuzzyClusterAlg::Flame_UCPearsonDist ,
        cluster::fuzzyClusterAlg::Flame_SQPearsonDist ,
        cluster::fuzzyClusterAlg::Flame_DotProductDist ,
        cluster::fuzzyClusterAlg::Flame_CovarianceDist ,
        cluster::fuzzyClusterAlg::Flame_Manhattan 
};

double cluster::fuzzyClusterAlg::Flame_Euclidean( double *x, double *y, int m )
{
	double d = 0;
	int i;
	for(i=0; i<m; i++ ) d += ( x[i] - y[i] ) * ( x[i] - y[i] );
	return sqrt( d );
}
double cluster::fuzzyClusterAlg::Flame_Cosine( double *x, double *y, int m )
{
	double r =0, x2 =0, y2 =0;
	int i;
	for(i=0; i<m; i++ ){
		r += x[i] * y[i];
		x2 += x[i] * x[i];
		y2 += y[i] * y[i];
	}
	return r / ( sqrt( x2 * y2 ) + EPSILON );
}
double cluster::fuzzyClusterAlg::Flame_Pearson( double *x, double *y, int m )
{
	int i;
	double r, x2, y2, xavg, yavg;
	if( m ==0 ) return 0;
	xavg = yavg = 0;
	r = x2 = y2 = 0;
	for( i=0; i<m; i++ ){
		xavg += x[i];
		yavg += y[i];
	}
	xavg = xavg/m;
	yavg = yavg/m;
	for( i=0; i<m; i++ ){
		r  += ( x[i] - xavg ) * ( y[i] - yavg );
		x2 += ( x[i] - xavg ) * ( x[i] - xavg );
		y2 += ( y[i] - yavg ) * ( y[i] - yavg );
	}
	return r / ( sqrt( x2 * y2 ) + EPSILON );
}
double cluster::fuzzyClusterAlg::Flame_UCPearson( double *x, double *y, int m )
{
	int i;
	double r, x2, y2, xavg, yavg;
	if( m ==0 ) return 0;
	xavg = yavg = 0;
	r = x2 = y2 = 0;
	for( i=0; i<m; i++ ){
		xavg += x[i];
		yavg += y[i];
	}
	xavg = xavg/m;
	yavg = yavg/m;
	for( i=0; i<m; i++ ){
		r  += x[i] * y[i];
		x2 += ( x[i] - xavg ) * ( x[i] - xavg );
		y2 += ( y[i] - yavg ) * ( y[i] - yavg );
	}
	return r / ( sqrt( x2 * y2 ) + EPSILON );
}
double cluster::fuzzyClusterAlg::Flame_SQPearson( double *x, double *y, int m )
{
	int i;
	double r, x2, y2, xavg, yavg;
	if( m ==0 ) return 0;
	xavg = yavg = 0;
	r = x2 = y2 = 0;
	for( i=0; i<m; i++ ){
		xavg += x[i];
		yavg += y[i];
	}
	xavg = xavg/m;
	yavg = yavg/m;
	for( i=0; i<m; i++ ){
		r  += ( x[i] - xavg ) * ( y[i] - yavg );
		x2 += ( x[i] - xavg ) * ( x[i] - xavg );
		y2 += ( y[i] - yavg ) * ( y[i] - yavg );
	}
	return r*r / ( x2 * y2 + EPSILON );
}
double cluster::fuzzyClusterAlg::Flame_DotProduct( double *x, double *y, int m )
{
	int i;
	double r = 0;
	for(i=0; i<m; i++ ) r += x[i] * y[i];
	if( m == 0 ) return 0;
	return r / m;
}
double cluster::fuzzyClusterAlg::Flame_Covariance( double *x, double *y, int m )
{
	int i;
	double r, x2, y2, xavg, yavg;
	if( m ==0 ) return 0;
	xavg = yavg = 0;
	r = x2 = y2 = 0;
	for( i=0; i<m; i++ ){
		xavg += x[i];
		yavg += y[i];
	}
	xavg = xavg/m;
	yavg = yavg/m;
	for( i=0; i<m; i++ ) r += ( x[i] - xavg ) * ( y[i] - yavg );
	if( m <= 1 ) return 0;
	return r / (m-1);
}
double cluster::fuzzyClusterAlg::Flame_Manhattan( double *x, double *y, int m )
{
	double d = 0;
	int i;
	for(i=0; i<m; i++ ) d += fabs( x[i] - y[i] );
	return d;
}
double cluster::fuzzyClusterAlg::Flame_CosineDist( double *x, double *y, int m )
{
	return 1-Flame_Cosine( x, y, m );
}
double cluster::fuzzyClusterAlg::Flame_PearsonDist( double *x, double *y, int m )
{
	return 1-Flame_Pearson( x, y, m );
}
double cluster::fuzzyClusterAlg::Flame_UCPearsonDist( double *x, double *y, int m )
{
	return 1-Flame_UCPearson( x, y, m );
}
double cluster::fuzzyClusterAlg::Flame_SQPearsonDist( double *x, double *y, int m )
{
	return 1-Flame_SQPearson( x, y, m );
}
double cluster::fuzzyClusterAlg::Flame_DotProductDist( double *x, double *y, int m )
{
	return 1-Flame_DotProduct( x, y, m );
}
double cluster::fuzzyClusterAlg::Flame_CovarianceDist( double *x, double *y, int m )
{
	return 1-Flame_Covariance( x, y, m );
}

cluster::fuzzyClusterAlg::Flame* cluster::fuzzyClusterAlg::Flame_New()
{
	Flame *self = (Flame*) malloc( sizeof(Flame) );
	memset( self, 0, sizeof(Flame) );
	return self;
}
void cluster::fuzzyClusterAlg::Flame_Clear( Flame *self )
{
	int i;
	for(i=0; i<self->N; i++){
		free( self->graph[i] );
		free( self->dists[i] );
		free( self->weights[i] );
		free( self->fuzzyships[i] );
	}
	if( self->clusters ){
		for(i=0; i<=self->cso_count; i++){
			if( self->clusters[i].array ) free( self->clusters[i].array );
		}
		free( self->clusters );
		self->clusters = NULL;
	}
	if( self->graph ) free( self->graph );
	if( self->dists ) free( self->dists );
	if( self->nncounts ) free( self->nncounts );
	if( self->weights ) free( self->weights );
	if( self->fuzzyships ) free( self->fuzzyships );
	if( self->obtypes ) free( self->obtypes );
	self->graph = NULL;
	self->dists = NULL;
	self->nncounts = NULL;
	self->weights = NULL;
	self->obtypes = NULL;
	self->fuzzyships = NULL;
	self->N = self->K = self->KMAX = self->cso_count = 0;
}

/* If m==0, data is distance matrix. */
void cluster::fuzzyClusterAlg::Flame_SetMatrix( Flame *self, double *data[], int n, int m )
{
	int i, j;
	int MAX = sqrt( n ) + 10;
	Indexdouble *vals = (Indexdouble*) calloc( n, sizeof(Indexdouble) );
	if( MAX >= n ) MAX = n - 1;
		
	Flame_Clear( self );
	self->N = n;
	self->KMAX = MAX;
	
	self->graph = (int**) calloc( n, sizeof(int*) );
	self->dists = (double**) calloc( n, sizeof(double*) );
	self->weights = (double**) calloc( n, sizeof(double*) );
	self->nncounts = (int*) calloc( n, sizeof(int) );
	self->obtypes = (char*) calloc( n, sizeof(char) );
	self->fuzzyships = (double**) calloc( n, sizeof(double*) );

	for(i=0; i<n; i++){
		self->graph[i] = (int*) calloc( MAX, sizeof(int) );
		self->dists[i] = (double*) calloc( MAX, sizeof(double) );
		self->weights[i] = (double*) calloc( MAX, sizeof(double) );
		if( m ==0 ){
			/* data is distance matrix. */
			for(j=0; j<n; j++){
				vals[j].index = j;
				vals[j].value = data[i][j];
			}
		}else{
			/* data is raw data matrix. */
			for(j=0; j<n; j++){
				vals[j].index = j;
				vals[j].value = self->distfunc( data[i], data[j], m );
			}
		}
		PartialQuickSort( vals, 0, n-1, MAX+1 );
		/* Store MAX number of nearest neighbors. */
		for(j=0; j<MAX; j++){
			self->graph[i][j] = vals[j+1].index;
			self->dists[i][j] = vals[j+1].value;
		}
	}
	free( vals );
}
void cluster::fuzzyClusterAlg::Flame_SetDataMatrix( Flame *self, double *data[], int n, int m, int dt )
{
	self->simtype = dt;
	if( dt >0 && dt < DST_NULL ) self->distfunc = basicDistFunctions[ dt-1 ];
	if( self->distfunc == NULL ) self->distfunc = basicDistFunctions[0];
	Flame_SetMatrix( self, data, n, m );
}
void cluster::fuzzyClusterAlg::Flame_SetDistMatrix( Flame *self, double *data[], int n )
{
	cluster::fuzzyClusterAlg::Flame_SetMatrix( self, data, n, DST_USER );
}
void cluster::fuzzyClusterAlg::Flame_DefineSupports( Flame *self, int knn, double thd )
{
	int i, j, k;
	int n = self->N;
	int kmax = self->KMAX;
	double **dists = self->dists;
	double *density = (double*) calloc( n, sizeof(double) );
	double d, sum, sum2, fmin, fmax = 0.0;
	
	if( knn > kmax ) knn = kmax;
	self->K = knn;
	for(i=0; i<n; i++){
		/* To include all the neighbors that have distances equal to the
		 * distance of the most distant one of the K-Nearest Neighbors */
		k = knn;
		d = dists[i][knn-1];
		for(j=knn; j<kmax; j++) if( dists[i][j] == d ) k ++; else break;
		self->nncounts[i] = k;

		/* The definition of weights in this implementation is 
		 * different from the previous implementations where distances 
		 * or similarities often have to be transformed in some way.
		 *
		 * But in this definition, the weights are only dependent on 
		 * the ranking of distances of the neighbors, so it is more 
		 * robust against distance transformations. */
		sum = 0.5*k*(k+1.0);
		for(j=0; j<k; j++) self->weights[i][j] = (k-j) / sum;
		
		sum = 0.0;
		for(j=0; j<k; j++) sum += dists[i][j];
		density[i] = 1.0 / (sum + EPSILON);
	}
	sum = 0.0;
	sum2 = 0.0;
	for(i=0; i<n; i++){
		sum += density[i];
		sum2 += density[i] * density[i];
	}
	sum = sum / n;
	/* Density threshold for possible outliers. */
	thd = sum + thd * sqrt( sum2 / n - sum * sum );

	memset( self->obtypes, 0, n*sizeof(char) );
	self->cso_count = 0;
	for(i=0; i<n; i++){
		k = self->nncounts[i];
		fmax = 0.0;
		fmin = density[i] / density[ self->graph[i][0] ];
		for(j=1; j<k; j++){
			d = density[i] / density[ self->graph[i][j] ];
			if( d > fmax ) fmax = d;
			if( d < fmin ) fmin = d;
			/* To avoid defining neighboring objects or objects close 
			 * to an outlier as CSOs.  */
			if( self->obtypes[ self->graph[i][j] ] ) fmin = 0.0;
		}
		if( fmin >= 1.0 ){
			self->cso_count ++;
			self->obtypes[i] = OBT_SUPPORT;
		}else if( fmax <= 1.0 && density[i] < thd ){
			self->obtypes[i] = OBT_OUTLIER;
		}
	}
	free( density );
}
void cluster::fuzzyClusterAlg::Flame_LocalApproximation( Flame *self, int steps, double epsilon )
{
	int i;
        int j;
        int k; 
        int t;
	int n = self->N;
        int m = self->cso_count;
	double **fuzzyships = self->fuzzyships;
	double **fuzzyships2 = (double**)calloc( n, sizeof(double*) );
	char *obtypes = self->obtypes;
	char even = 0;
	double dev;

	k = 0;
	for(i=0; i<n; i++){
		fuzzyships[i] = (double*) realloc( fuzzyships[i], (m+1)*sizeof(double) );
		fuzzyships2[i] = (double*) calloc( m+1, sizeof(double) );
		memset( fuzzyships[i], 0, (m+1)*sizeof(double) );
		if( obtypes[i] == OBT_SUPPORT ){
			/* Full membership to the cluster represented by itself. */
			fuzzyships[i][k] = 1.0;
			fuzzyships2[i][k] = 1.0;
			k ++;
		}else if( obtypes[i] == OBT_OUTLIER ){
			/* Full membership to the outlier group. */
			fuzzyships[i][m] = 1.0;
			fuzzyships2[i][m] = 1.0;
		}else{
			/* Equal memberships to all clusters and the outlier group.
			 * Random initialization does not change the results. */
			for(j=0; j<=m; j++)
				fuzzyships[i][j] = fuzzyships2[i][j] = 1.0/(m+1);
		}
	}
	for(t=0; t<steps; t++){
		dev = 0;
		for(i=0; i<n; i++){
			int knn = self->nncounts[i];
			int *ids = self->graph[i];
			double *wt = self->weights[i];
			double *fuzzy = fuzzyships[i];
			double **fuzzy2 = fuzzyships2;
			double sum = 0.0;
			if( self->obtypes[i] != OBT_NORMAL ) continue;
			if( even ){
				fuzzy = fuzzyships2[i];
				fuzzy2 = fuzzyships;
			}
			/* Update membership of an object by a linear combination of 
			 * the memberships of its nearest neighbors. */
			for(j=0; j<=m; j++){
				fuzzy[j] = 0.0;
				for(k=0; k<knn; k++) fuzzy[j] += wt[k] * fuzzy2[ ids[k] ][j];
				dev += (fuzzy[j] - fuzzy2[i][j]) * (fuzzy[j] - fuzzy2[i][j]);
				sum += fuzzy[j];
			}
			for(j=0; j<=m; j++) fuzzy[j] = fuzzy[j] / sum;
		}
		even = ! even;
		if( dev < epsilon ) break;
	}
	/* update the membership of all objects to remove clusters 
	 * that contains only the CSO. */
	for(i=0; i<n; i++){
		int knn = self->nncounts[i];
		int *ids = self->graph[i];
		double *wt = self->weights[i];
		double *fuzzy = fuzzyships[i];
		double **fuzzy2 = fuzzyships2;
		for(j=0; j<=m; j++){
			fuzzy[j] = 0.0;
			for(k=0; k<knn; k++) fuzzy[j] += wt[k] * fuzzy2[ ids[k] ][j];
			dev += (fuzzy[j] - fuzzy2[i][j]) * (fuzzy[j] - fuzzy2[i][j]);
		}
	}
	for(i=0; i<n; i++) free( fuzzyships2[i] );
	free( fuzzyships2 );
}

void cluster::fuzzyClusterAlg::IntArray_Push( IntArray *self, int value )
{
	if( self->size >= self->bufsize ){
		self->bufsize += self->bufsize /10 + 10;
		self->array = (int*)realloc( self->array, self->bufsize*sizeof(int));
	}
	self->array[ self->size ] = value;
	self->size ++;
}
void cluster::fuzzyClusterAlg::Flame_MakeClusters( Flame *self, double thd )
{
	int i, j, imax;
	int N = self->N;
	int C = self->cso_count+1;
	double fmax;
	double **fuzzyships = self->fuzzyships;
	//IntArray *clust;
	Indexdouble *vals = (Indexdouble*) calloc( N, sizeof(Indexdouble) );
	
	/* Sort objects based on the "entropy" of fuzzy memberships. */
	for(i=0; i<N; i++){
		vals[i].index = i;
		vals[i].value = 0.0;
		for(j=0; j<C; j++){
			double fs = fuzzyships[i][j];
			if( fs > EPSILON ) vals[i].value -= fs * log( fs );
		}
	}
	PartialQuickSort( vals, 0, N-1, N );

	if( self->clusters ){
		for(i=0; i<C; i++)
			if( self->clusters[i].array ) free( self->clusters[i].array );
		free( self->clusters );
	}
	self->clusters = (IntArray*) calloc( C, sizeof(IntArray) );
	if( thd <0 || thd > 1.0 ){
		/* Assign each object to the cluster 
		 * in which it has the highest membership. */
		for(i=0; i<N; i++){
			int id = vals[i].index;
			fmax = 0;
			imax = -1;
			for(j=0; j<C; j++){
				if( fuzzyships[id][j] > fmax ){
					imax = j;
					fmax = fuzzyships[id][j];
				}
			}
			IntArray_Push( self->clusters + imax, id );
		}
	}else{
		/* Assign each object to all the clusters
		 * in which it has membership higher than thd,
		 * otherwise, assign it to the outlier group.*/
		for(i=0; i<N; i++){
			int id = vals[i].index;
			imax = -1;
			for(j=0; j<C; j++){
				if( fuzzyships[id][j] > thd || ( j == C-1 && imax <0 ) ){
					imax = j;
					//clust = self->clusters + j;
					IntArray_Push( self->clusters + j, id );
				}
			}
		}
	}
	/* removing empty clusters */
	C = 0;
	for(i=0; i<self->cso_count; i++){
		if( self->clusters[i].size >0 ){
			self->clusters[C] = self->clusters[i];
			C ++;
		}
	}
	/* keep the outlier group, even if its empty */
	self->clusters[C] = self->clusters[self->cso_count];
	C ++;
	for(i=C; i<self->cso_count+1; i++) memset( self->clusters+i, 0, sizeof(IntArray) );
	self->count = C;
	free( vals );
}






