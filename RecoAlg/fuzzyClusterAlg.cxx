////////////////////////////////////////////////////////////////////////
//
// fuzzyClusterAlg.cxx
//
// Ben Carls, bcarls@fnal.gov
//
// This code looks for clusters using a Gustafson-Kessel variant on fuzzy c-means algorithm. The
// clusters are then examined by the HoughBaseAlg to identify Hough lines
// which can then be split off into their own clusters. See the webpage below
// for more information on the fuzzy clustering algorithm.
//
//http://homes.di.unimi.it/~valenti/SlideCorsi/Bioinformatica05/Fuzzy-Clustering-lecture-Babuska.pdf
//http://biosoft.kaist.ac.kr/BISL_homepage/publication/p20050006.pdf
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
  fFuzzifier                      = p.get< double  >("Fuzzifier");
  fMaxNumClusters                 = p.get< int    >("MaxNumClusters");
  nIterations                     = p.get< int    >("Iterations");
  fDoFuzzyRemnantMerge            = p.get< int    >("DoFuzzyRemnantMerge"            );
  fFuzzyRemnantMergeCutoff        = p.get< double >("FuzzyRemnantMergeCutoff"        );
  fDoTrackClusterMerge            = p.get< int    >("DoTrackClusterMerge"            );
  fTrackClusterMergeCutoff        = p.get< double >("TrackClusterMergeCutoff"        );
  fChargeAsymAngleCut             = p.get< double >("ChargeAsymAngleCut"             );
  fSigmaChargeAsymAngleCut        = p.get< double >("SigmaChargeAsymAngleCut"        );
  fDoShowerClusterMerge           = p.get< int    >("DoShowerClusterMerge"           );
  fDoShowerTrackClusterMerge      = p.get< int    >("DoShowerTrackClusterMerge"      );
  fShowerClusterMergeCutoff       = p.get< double >("ShowerClusterMergeCutoff"       );
  fShowerClusterMergeAngle        = p.get< double >("ShowerClusterMergeAngle"        );
  fShowerTrackClusterMergeCutoff  = p.get< double >("ShowerTrackClusterMergeCutoff"  );
  fShowerTrackClusterMergeAngle   = p.get< double >("ShowerTrackClusterMergeAngle"   );
  fShowerLikenessCut              = p.get< double >("ShowerLikenessCut"              );
  fMaxVertexLines                 = p.get< int   >("MaxVertexLines"                 );
  fVertexLinesCutoff              = p.get< double >("VertexLinesCutoff"              );
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

  // Clear the matrix that stores the points needed for fuzzy clustering
  fpsMat.Clear();
  fpsMembership.Clear();
  fpsNewMembership.Clear();
  fpsCentroids.Clear(); 

  //------------------------------------------------------------------
  // Determine spacing between wires (different for each detector)
  ///get 2 first wires and find their spacing (wire_dist)

  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<util::DetectorProperties> detp;
 
  fLargestSijLast=1;

  // Collect the hits in a useful form,
  // and take note of the maximum time width
  fMaxWidth=0.0;
  //fpsMat = TMatrixT<double>(allhits.size(),2);
  //fpsMembership = TMatrixT<double>(iNumClusters, allhits.size());
  fpsMat.ResizeTo(allhits.size(),2);
  double tickToDist = larp->DriftVelocity(larp->Efield(),larp->Temperature());
  tickToDist *= 1.e-3 * detp->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
  int dims = 3;//our point is defined by 3 elements:wire#,center of the hit, and the hit width
  std::vector<double> p(dims);
  for (auto allhitsItr = allhits.begin(); allhitsItr < allhits.end(); allhitsItr++){
        
    p[0] = ((*allhitsItr)->Channel())*fGeom->WirePitch(fGeom->View((*allhitsItr)->Channel()));
    p[1] = (((*allhitsItr)->StartTime()+(*allhitsItr)->EndTime()  )/2.)*tickToDist;
    p[2] = ((*allhitsItr)->EndTime()  -(*allhitsItr)->StartTime())*tickToDist;   //width of a hit in cm

    // check on the maximum width condition
    if ( p[2] > fMaxWidth ) fMaxWidth = p[2];
    
    fps.push_back(p);

    // Store hits in the matrix needed for fuzzy clustering
    fpsMat(allhitsItr-allhits.begin(),0) = p[0];
    fpsMat(allhitsItr-allhits.begin(),1) = p[1];

  }

  mf::LogInfo("fuzzyCluster") << "InitFuzzy: hits vector size is " << fps.size();
  
  return;
}


//----------------------------------------------------------
inline void cluster::fuzzyClusterAlg::computeCentroids(int k)
{
  // fpsCentroids are the weighted centers of the clusters.
  // We multiply fpsMembership by fpsMat to find fpsCentroids, then
  // divide by the normalization (sum of the weights).

  int iNumClusters = k;

  //fpsMembership.Print();
  //fpsCentroids.Print();
  //fpsMat.Print();

  fpsCentroids.ResizeTo(k,2);

  // Zero the centroid matrix
  for ( int j = 0; j < fpsCentroids.GetNrows(); j++)
    for ( int i = 0; i < fpsCentroids.GetNcols(); i++)
      fpsCentroids(j,i) = 0;

  // Determine the elements of u^m_ij
  TMatrixT<double> Uji_m(iNumClusters, fpsMat.GetNrows());
  // For each cluster
  for ( int j = 0; j < iNumClusters; ++j)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); ++i)
      // Determine Uji_m
      Uji_m(j,i) = pow(fpsMembership(j,i), fFuzzifier); 


  // Now find sum^N_i=1 u^m_ij*x_i
  // For each cluster
  for ( int j = 0; j < iNumClusters; ++j)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); ++i)
      // For each dimension
      for ( int f = 0; f < 2; f++)
        fpsCentroids(j,f) = fpsCentroids(j,f) + Uji_m(j,i)*fpsMat(i,f);


  // Divide centroids by the normalization (sum of weights)
  // For each cluster
  for ( int j = 0; j < iNumClusters; ++j){
    double normalizationFactor = 0;
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); ++i)
      normalizationFactor += Uji_m(j,i);
    // For each dimension
    for ( int f = 0; f < 2; f++)
      fpsCentroids(j,f) /= normalizationFactor;
  }

}
//----------------------------------------------------------
inline bool cluster::fuzzyClusterAlg::updateMembership(int *k)
{
  // We are updating the membership of the data points based on the centroids
  // determined. This is determined using
 
  //std::cout << "New updateMembership" << std::endl;

  //fpsMat.Print();
  //fpsCentroids.Print();
  //fpsMembership.Print();
  
  std::vector<TMatrixT<double>> clusterCovarianceMats(*k);
  std::vector<double> clusterRadii(*k);

  // Determine the elements of u^m_ij
  TMatrixT<double> Uji_m(*k, fpsMat.GetNrows());
  // For each cluster
  for ( int j = 0; j < *k; ++j)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); ++i)
      // Determine Uji_m
      Uji_m(j,i) = pow(fpsMembership(j,i), fFuzzifier); 



  TVectorT<double> fpsMat_row(2);
  TVectorT<double> fpsCentroids_row(2);
  TVectorT<double> fpsMat_col(2);
  TVectorT<double> fpsCentroids_col(2);
  TVectorT<double> fpsMat_row_t(2);
  TVectorT<double> fpsCentroids_row_t(2);
  TMatrixT<double> fpsMatMinusCent_col(2,1);
  TMatrixT<double> fpsMatMinusCent_row(1,2);
  TMatrixT<double> fpsDistances(*k,fpsMat.GetNrows());
  TMatrixT<double> clusCovarianceMat(2,2);
  // Calculate the covariance matrix
  // For each cluster
  double Uji_m_sum;
  for ( int j = 0; j < *k; ++j){
    fpsCentroids_row = TMatrixDRow(fpsCentroids,j);
    Uji_m_sum = 0;
    //For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++){
      fpsMat_row = TMatrixDRow(fpsMat,i);
      fpsMatMinusCent_col(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_col(1,0)=fpsMat_row(1)-fpsCentroids_row(1);
      fpsMatMinusCent_row(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_row(0,1)=fpsMat_row(1)-fpsCentroids_row(1);
      //std::cout << "fpsMat" << std::endl;
      //fpsMat_row.Print();
      //std::cout << "fpsCentroids" << std::endl;
      //fpsCentroids_row.Print();
      clusCovarianceMat += Uji_m(j,i)*fpsMatMinusCent_col*fpsMatMinusCent_row;
      Uji_m_sum+=Uji_m(j,i);
    }
    //std::cout << "Uji_m_sum: " << Uji_m_sum << std::endl;
    clusCovarianceMat=clusCovarianceMat*(1/Uji_m_sum);
    //clusCovarianceMat.Print();
    clusterCovarianceMats[j].ResizeTo(clusCovarianceMat);
    clusterCovarianceMats[j] = clusCovarianceMat;

    // Is the cluster covariance matrix singular?
    double clusCovarianceMatDeterminant = clusCovarianceMat(0,0)*clusCovarianceMat(1,1)-clusCovarianceMat(0,1)*clusCovarianceMat(1,0);

    // Check if the covariance matrix determinant is zero or nan
    if(clusCovarianceMatDeterminant > 0 
        && std::isfinite(clusCovarianceMatDeterminant)
        && std::isnormal(clusCovarianceMatDeterminant)){
      clusterRadii[j] = fBeta*pow(clusCovarianceMatDeterminant,0.25)/((double)*k);
    }
    else{
      mf::LogVerbatim("fuzzyCluster") << "updateMembership: Covariance matrix is singular 1";
      clusterRadii[j]=0;
      continue;
    }
    //std::cout << "Made it past exception 1" << std::endl;


    //try
    //{
      //clusterRadii[j] = fBeta*pow(clusCovarianceMat.Determinant(),0.25)/((double)*k);
    //}
    //catch(...){
      //mf::LogVerbatim("fuzzyCluster") << "updateMembership: Covariance matrix is singular 1";
      //continue;
    //}
  }

  
  for ( int j = 0; j < *k; ++j){
    bool clusterCovarianceSingular = false;

    //TMatrixT<double>& clusCovarianceMatInv = clusterCovarianceMats[j];
    
    TDecompSVD clusCovarianceMatInvSVD(clusterCovarianceMats[j]);
    clusCovarianceMatInvSVD.SetTol(1e-20);
    TMatrixT<double> clusCovarianceMatInv = clusCovarianceMatInvSVD.Invert();
    

    //std::cout << "cov. matrix: " << clusCovarianceMatInv.Determinant() << std::endl;
    double clusCovarianceMatInvDeterminant = clusCovarianceMatInv(0,0)*clusCovarianceMatInv(1,1)-clusCovarianceMatInv(0,1)*clusCovarianceMatInv(1,0);



    //std::cout << "inverse cov. matrix: " << clusCovarianceMatInv.Determinant() << " " << std::endl;
    fpsCentroids_row = TMatrixDRow(fpsCentroids,j);
    //TMatrixT<double> fpsMatMinusCent_row(1,2);
    //TMatrixT<double> fpsMatMinusCent_col(2,1);
    double clusCovarianceMatInvDet = clusCovarianceMatInvDeterminant;
    for ( int i = 0; i < fpsMat.GetNrows(); ++i){
      fpsMat_row = TMatrixDRow(fpsMat,i);
      fpsMatMinusCent_row(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_row(0,1)=fpsMat_row(1)-fpsCentroids_row(1);
      fpsMatMinusCent_col(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_col(1,0)=fpsMat_row(1)-fpsCentroids_row(1);
      TMatrixT<double> tempDistanceSquared = (fpsMatMinusCent_row*(clusCovarianceMatInv*fpsMatMinusCent_col));
      if(clusCovarianceMatInvDeterminant != 0 
        && std::isfinite(clusCovarianceMatInvDeterminant)
        && std::isnormal(clusCovarianceMatInvDeterminant)){
          //fpsDistances(j,i) = std::sqrt(std::max((double)0,tempDistanceSquared(0,0)/sqrt(clusterCovarianceMats[j].Determinant()) - pow(clusterRadii[j],2)));
          fpsDistances(j,i) = std::sqrt(std::max((double)0,tempDistanceSquared(0,0)/std::sqrt(clusCovarianceMatInvDet) - (clusterRadii[j])*(clusterRadii[j])));
      }
      else{
        clusterCovarianceSingular = true;
        fpsDistances(j,i) = 999999;
      }
      //std::cout << "fpsDistances(j,i): " << fpsDistances(j,i) << std::endl;
    }
    if(clusterCovarianceSingular)
      mf::LogVerbatim("fuzzyCluster") << "updateMembership: Covariance matrix is singular 3";
  }

  //Determine the new elements of u_ij
  //std::cout << "fpsDistances: " << std::endl;
  //fpsDistances.Print();
  fpsNewMembership.ResizeTo(fpsMembership);
  double fCoeff;
  int clusWithZeroDistance = 0;
  // For each hit
  for ( int i = 0; i < fpsMat.GetNrows(); ++i){
    // Does the hit have nonzero distance to all clusters?
    clusWithZeroDistance = 0;
    for ( int j = 0; j < *k; ++j){
      if(fpsDistances(j,i) == 0)
        clusWithZeroDistance++;
    }
    // For each cluster
    for ( int j = 0; j < *k; ++j){
      if(clusWithZeroDistance == 0){ 
        fCoeff = 0;
        // For each cluster
        for ( int l = 0; l < *k; ++l){
          //std::cout << fpsDistances(j,i) << " " << fpsDistances(k,i) << std::endl;
          fCoeff += pow((fpsDistances(j,i)/fpsDistances(l,i)),2/(fFuzzifier - 1));
        }
        //std::cout << "fCoeff: " << fCoeff << std::endl;
        fpsNewMembership(j,i) = 1/fCoeff;
      }
      if(clusWithZeroDistance>0){
        if(fpsDistances(j,i) > 0){
          fpsNewMembership(j,i) = 0;
        }
        if(fpsDistances(j,i) == 0){
          fpsNewMembership(j,i) = 1/(double)clusWithZeroDistance;
          //fpsNewMembership(j,i) = 1/(double)iNumClusters;
        }
      }
    }
  }
  //fpsNewMembership.Print();





  if(!canStop()){
    fpsMembership = fpsNewMembership;
    return false;
  }
  return true;

}  




//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
//
//  Ben Carls' implementation of fuzzyClusterAlg as much like examples as possible
void cluster::fuzzyClusterAlg::run_fuzzy_cluster(std::vector<art::Ptr<recob::Hit> >& allhits) {

  int nMaxClusters = fMaxNumClusters;
  if(fpsMat.GetNrows() < nMaxClusters)
    nMaxClusters = fpsMat.GetNrows();
  
  std::vector<TMatrixD> fpsMembershipStore;
  fpsMembershipStore.clear();
  fpsMembershipStore.resize(nMaxClusters);

  fpointId_to_clusterId.resize(fps.size(), kNO_CLUSTER); // Not zero as before!
  fnoise.resize(fps.size(), false);
  fvisited.resize(fps.size(), false);

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;

  if(allhits.size()==0)
    return;
  
  //factor to make x and y scale the same units
  uint32_t     channel = allhits[0]->Wire()->RawDigit()->Channel();
  double wirePitch = geom->WirePitch(geom->View(channel));
  double xyScale  = .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
  xyScale        *= detprop->SamplingRate()/wirePitch;
  double wire_dist = wirePitch;
  double tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
  tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
  
  double indcolscaling = 0.;       //a parameter to account for the different 
        			   ////characteristic hit width of induction and collection plane
  /// \todo: the collection plane's characteristic hit width's are, 
  /// \todo: on average, about 5 time samples wider than the induction plane's. 
  /// \todo: this is hard-coded for now.
  geo::SigType_t sigt = geom->SignalType(channel);
  if(sigt == geo::kInduction)
    indcolscaling = 0.;
  else
    indcolscaling = 1.;
  
  //fpsMat.Print();


  int k = nMaxClusters;
  if (k > fpsMat.GetNrows() || k <= 0)
    return;

  fpsMembership.ResizeTo(k, fps.size());
  fpsNewMembership.ResizeTo(k, fps.size());
  fpsMembershipStore[k-1].ResizeTo(k, fps.size());
  fpsCentroids.ResizeTo(k,2);


  /// Get the random number generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine & engine = rng -> getEngine();
  CLHEP::RandFlat flat1(engine);

  //Randomize membership for each hit for fuzzy
  double normalizationFactor;
  for( int i = 0; i < fpsMat.GetNrows(); ++i){
    normalizationFactor = 0;
    for( int j = 0; j < k; j++){
      fpsMembership(j,i) = flat1.fire(); //(rand() / (RAND_MAX + 0.0));
      normalizationFactor += fpsMembership(j,i);
    }
    for( int j = 0; j < k; j++)
      fpsMembership(j,i) /= normalizationFactor;
  }

  // Compute initial centroids
  computeCentroids(k);


  //for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid ;
    //for (int l = 0; l < k; l++)
      //mf::LogInfo("fuzzyCluster") << l  << fpsMembership(l,pid) ;
  //} 

  // Run iterations of the clustering algorithm
  fBeta = 1;
  int i = 0;
  while(!updateMembership(&k)){
    //std::cout << "k: " << k << std::endl;
    if(k == 1)
      break;
    if(mergeClusters()){
      k--;
      //fpsNewMembership.ResizeTo(fpsMembership);
      //fpsNewMembership = fpsMembership;
    }
    computeCentroids(k);
    if(i + 1 == nIterations)
      break;
    i++;
  }
  
  //for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid ;
    //for (int l = 0; l < k; l++)
      //mf::LogInfo("fuzzyCluster") << l  << fpsMembership(l,pid) ;
  //} 





  fpsMembershipFinal.ResizeTo(k, fps.size());
  fpsMembershipFinal = fpsMembership; 











 

  int nClusters = 0;
  if(k > 0) nClusters = fpsMembershipFinal.GetNrows();
  unsigned int cid = fpsMembershipFinal.GetNrows();
  //int nClusters = iMinXBClusterNum;
  //unsigned int cid = iMinXBClusterNum;
  //mf::LogInfo("fuzzyCluster") << "Number of clusters found after merging: " << nClusters   ;
  //int nClusters = 2;
  //unsigned int cid = 2;

  //mf::LogInfo("fuzzyCluster") << iMinXBClusterNum  << nClusters ;
  //std::cout << "nClusters: " << nClusters  << std::endl;
  int iCluster;
  double maxClusMembership;
  for (size_t pid = 0; pid < fps.size(); ++pid){
    //mf::LogInfo("fuzzyCluster") << pid ;
    iCluster = kNO_CLUSTER;
    maxClusMembership = -1;
    // not already visited
    if (fpointId_to_clusterId[pid] == kNO_CLUSTER) {
      for (int i = 0; i <= nClusters-1; ++i){
        //mf::LogInfo("fuzzyCluster") << i  << fpsMembershipStore[nClusters-1](i,pid) ;
        if ( fpsMembershipFinal(i,pid) > maxClusMembership ) {
          maxClusMembership = fpsMembershipFinal(i,pid); 
          iCluster = i;
        }
      }
    } // if (!visited
    fpointId_to_clusterId[pid] = iCluster;
  } // for

  // Run EndPointAlg over hits to see if 
  //mf::LogVerbatim("fuzzyCluster") << "New plane: " ;
  std::vector<unsigned int> corners;
  corners.clear();
  // nClustersTEmp is how many fuzzy clusters we originally found, it is not how many hough lines we found
  unsigned int nClustersTemp = nClusters;

 
  




  
  // Loop over clusters with the Hough line finder to break the clusters up further
  // list of lines
  //std::cout << "Starting Hough" << std::endl;
  std::vector<protoTrack> protoTracksFound;
  //TStopwatch w;
  //double timeTotal = 0;
  if(nClustersTemp > 0)
    for (unsigned int i = 0; i <= (unsigned int)nClustersTemp-1; ++i){
      //w.Start();
      fHBAlg.Transform(allhits, &fpointId_to_clusterId, i, &nClusters, &protoTracksFound);
      //w.Stop();
      //timeTotal+=w.CpuTime();
    }
  //std::cout << "Hough over, took for cpu time: " << timeTotal << std::endl;

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
      distance = (TMath::Abs((*hitsItr)->PeakTime()-protoTracksFoundItr->clusterSlope*(double)((*hitsItr)->WireID().Wire)-protoTracksFoundItr->clusterIntercept)/(std::sqrt(pow(xyScale*protoTracksFoundItr->clusterSlope,2)+1)));
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
      trackMerged = mergeTrackClusters(i,&trackClusters,xyScale,wire_dist,tickToDist);
      if(trackMerged)
        continue;
      else
        i++;
    } 
  }

  if(fDoShowerClusterMerge && showerClusters.size() > 1){
    unsigned int i = 0;
    while(i < showerClusters.size()-1){ 
      showerMerged = mergeShowerClusters(i,&showerClusters,xyScale,wire_dist,tickToDist);
      if(showerMerged)
        continue;
      else
        i++;
    } 
  }

  if(fDoShowerTrackClusterMerge && showerClusters.size() > 0 && trackClusters.size() >0){
    unsigned int i = 0;
    while(i < showerClusters.size()){ 
      unsigned int j = 0;
      while(j < trackClusters.size()){
        showerTrackMerged = mergeShowerTrackClusters(&showerClusters[i],&trackClusters[j],xyScale,wire_dist,tickToDist);
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
      showerMerged = mergeShowerClusters(i,&showerClusters,xyScale,wire_dist,tickToDist);
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
  int nDBClusters = 0;
  bool unclustered;
  double p0;
  double p1;
  double minDistance;
  if(fDoFuzzyRemnantMerge){
    for(auto allhitsItr = allhits.cbegin(); allhitsItr != allhits.cend(); ++allhitsItr){
      unclustered = true;
      // nClusters is the number of fuzzy clusters we found, we only assign hits to lines here
      // if they are not already part of hough lines
      if(fpointId_to_clusterId.at(allhitsItr-allhits.begin()) >= (unsigned int) nClustersTemp){
        unclustered = false;
        continue;
      }
      p0 = ((*allhitsItr)->Wire()->RawDigit()->Channel())*wire_dist;
      p1 = (((*allhitsItr)->StartTime()+(*allhitsItr)->EndTime())/2.)*tickToDist;
      minDistance = 999999;

      for(auto showerClustersItr = showerClusters.begin(); showerClustersItr != showerClusters.end(); ++showerClustersItr){
        for(auto protoTracksItr = showerClustersItr->clusterProtoTracks.begin(); protoTracksItr < showerClustersItr->clusterProtoTracks.end(); ++protoTracksItr){
          
          distance = PointSegmentDistance( p0, p1, protoTracksItr->pMin0, protoTracksItr->pMin1, protoTracksItr->pMax0, protoTracksItr->pMax1);
          // Is the point behind or ahead of the line?
          //if(protoTracksFoundItr->pMin0 > p0){
             //double protoTracksFoundSlope = (protoTracksFoundItr->pMax1 - protoTracksFoundItr->pMin1)/(protoTracksFoundItr->pMax0 - protoTracksFoundItr->pMin0);
             //double pMinHitSlope = (p1 - protoTracksFoundItr->pMin1)/(p0 - protoTracksFoundItr->pMin0);
             //double slopeAngle = atan(std::abs((protoTracksFoundSlope - pMinHitSlope)/(1 + protoTracksFoundSlope*pMinHitSlope)))*(180/TMath::Pi());
             //if(distance < 10 && slopeAngle < 10){
               //fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = protoTracksFoundItr->clusterNumber;
               //unclustered = false;
               //break;
             //}
          //}
          //if (protoTracksFoundItr->pMax0 < p0){
             //double protoTracksFoundSlope = (protoTracksFoundItr->pMax1 - protoTracksFoundItr->pMin1)/(protoTracksFoundItr->pMax0 - protoTracksFoundItr->pMin0);
             //double pMaxHitSlope = (protoTracksFoundItr->pMax1-p1)/(protoTracksFoundItr->pMax0-p0);
             //double slopeAngle = atan(std::abs((protoTracksFoundSlope - pMaxHitSlope)/(1 + protoTracksFoundSlope*pMaxHitSlope)))*(180/TMath::Pi());
             //if(distance < 10 && slopeAngle < 10){
               //fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = protoTracksFoundItr->clusterNumber;
               //unclustered = false;
               //break;
             //}
          //}
          
          // If the line does not look showerlike, skip it
          if(protoTracksItr->showerLikeness<fShowerLikenessCut)
            continue;

          if(distance > fFuzzyRemnantMergeCutoff)
            continue;

          distance/=pow(protoTracksFoundSizes[protoTracksItr->clusterNumber],1/4);
          if(distance < minDistance){
            fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = protoTracksItr->clusterNumber;
            minDistance = distance;
            unclustered = false;
          }
        }
      }


      if(unclustered){
        unclusteredhitsToallhits.push_back(allhitsItr-allhits.begin());
        unclusteredhits.push_back(*allhitsItr);
      }
      
    }

    // Setup DBSCAN for noise and extra hits
    // Start by getting the ChannelFilter
    filter::ChannelFilter chanFilt;
    fDBScan.InitScan(unclusteredhits, chanFilt.SetOfBadChannels());
    fDBScan.run_cluster();

    nDBClusters = fDBScan.fclusters.size();
    for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){	  
      if (fDBScan.fpointId_to_clusterId[j]== kNO_CLUSTER || fDBScan.fpointId_to_clusterId[j]==kNOISE_CLUSTER) {
        fpointId_to_clusterId.at(unclusteredhitsToallhits[j]) = kNOISE_CLUSTER;
      } 
      else {
        fpointId_to_clusterId.at(unclusteredhitsToallhits[j]) = fDBScan.fpointId_to_clusterId[j] + nClusters;
      }
    }
  }













  cid = nClusters + nDBClusters;
  
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















//----------------------------------------------------------
bool cluster::fuzzyClusterAlg::mergeClusters()
{


  // If we only have one cluster, move on 
  if(fpsMembership.GetNrows() == 1)
    return false;

  //for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid  << fpsMembershipCrisp[pid]  << fpsMat(pid,0)  << fpsMat(pid,1) ;
  //}
  
  double largestSij=0;
  double largesti=-1;
  double largestj=-1;
  double sumMinUikUjk = 0;
  double sumUik = 0;
  double sumUjk = 0;
  for( int i = 0; i < fpsMembership.GetNrows()-1; ++i){
    for( int j = i+1; j < fpsMembership.GetNrows(); ++j){
      //std::cout << i << " " << j << std::endl;
      sumMinUikUjk = 0;
      sumUik = 0;
      sumUjk = 0;
      for( int k = 0; k < fpsMembership.GetNcols(); ++k){
        sumMinUikUjk+=std::min(fpsMembership(i,k),fpsMembership(j,k));
        sumUik += fpsMembership(i,k);
        sumUjk += fpsMembership(j,k);
      }
      double Sij = sumMinUikUjk/std::min(sumUik,sumUjk);
      if(Sij > largestSij){
        largestSij = Sij;
        largesti=i;
        largestj=j;
      }
    }
  }
 

  if(std::abs(fLargestSijLast-largestSij) > 0.1){
    fLargestSijLast=largestSij;
    return false;
  }
  
  //std::cout << "largest Sij: " << largestSij << std::endl;
  if( largestSij > 1/(fpsMembership.GetNrows()-1)){
    mf::LogVerbatim("fuzzyCluster") << "You're going to Merge!";
  }
  else {
    fBeta = std::min((double)fpsMembership.GetNrows(),fBeta+1);
    return false;
  } 

  if(largesti > -1 && largestj > -1){
    TMatrixD fpsMembershipTemp(fpsMembership.GetNrows()-1, fpsMat.GetNrows());
    TMatrixDRow(fpsMembership,largesti) += TMatrixDRow(fpsMembership,largestj);
    for( int j = 0; j < fpsMembership.GetNrows()-1; ++j){
      //std::cout << j << " " << largestj << " " << largesti << std::endl;
      if(j < largestj)  TMatrixDRow(fpsMembershipTemp,j) = TMatrixDRow(fpsMembership,j);
      if(j >= largestj) TMatrixDRow(fpsMembershipTemp,j) = TMatrixDRow(fpsMembership,j+1);
    }
    fpsMembership.ResizeTo(fpsMembershipTemp);
    fpsMembership=fpsMembershipTemp;
  }
  
  return true;

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
  pHit0[0] = (hit0->Wire()->RawDigit()->Channel())*wire_dist;
  pHit0[1] = ((hit0->StartTime()+hit0->EndTime())/2.)*tickToDist;
  double pHit1[2];
  pHit1[0] = (hit1->Wire()->RawDigit()->Channel())*wire_dist;
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



