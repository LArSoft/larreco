////////////////////////////////////////////////////////////////////////
//
// ClusterPCA class
//
// This algorithm is designed to find the principal axis and diffuseness 
//  of clusters
//  
////////////////////////////////////////////////////////////////////////

#include <string>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

//Framework includes:
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft includes:
#include "Geometry/Geometry.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"

// Root includes
#include "TPrincipal.h"
#include "TTree.h"



namespace cluster {
   
  class ClusterPCA : public art::EDAnalyzer {
    
  public:
    
    explicit ClusterPCA(fhicl::ParameterSet const& pset); 
    ~ClusterPCA();

    void PerformClusterPCA( std::vector<art::Ptr<recob::Hit> >& HitsThisCluster, double* PrincDirectionWT, double& PrincValue, double& TotalCharge, bool NormPC);
    
    void analyze(art::Event const& evt);
    void beginJob();
    
  private:
    
    std::string     fClusterModuleLabel;
    bool            fNormPC;
    
    TTree *         fTree;
    
    Int_t fView;
    Float_t fPrincDirW;
    Float_t fPrincDirT;
    Float_t fPrincValue;
    Float_t fTotalCharge;
    Float_t fNHits;
    
  }; // class ClusterPCA

}

//#endif 



namespace cluster{

  //-------------------------------------------------
  ClusterPCA::ClusterPCA(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset)
    , fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel"))
    , fNormPC            (pset.get<bool>("NormPC"))
  {
  }

  //-------------------------------------------------
  ClusterPCA::~ClusterPCA()
  {
  }

  //-------------------------------------------------
  // Set up analysis tree
  void ClusterPCA::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("PCATree","PCATree");
    fTree->Branch("View",      &fView,      "View/I");
    fTree->Branch("PrincDirW", &fPrincDirW, "PrincDirW/F");
    fTree->Branch("PrincDirT", &fPrincDirT, "PrincDirT/F");
    fTree->Branch("PrincValue",&fPrincValue, "PrincValue/F");
    fTree->Branch("TotalCharge",&fTotalCharge, "TotalCharge/F");
    fTree->Branch("NHits",      &fNHits,     "fNHits/F");
  }
    
  //------------------------------------------------------------------------------------//
  void ClusterPCA::analyze(art::Event const& evt)
  { 
    // Get a Handle for the input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    evt.getByLabel(fClusterModuleLabel,clusterVecHandle);

    art::ServiceHandle<geo::Geometry> geo;
    int nplanes = geo->Nplanes();

    //one PtrVector for each plane in the geometry
    std::vector< art::PtrVector<recob::Cluster> > Cls(nplanes);

    // loop over the input Clusters
    for(size_t i = 0; i < clusterVecHandle->size(); ++i){
      
      //get a art::Ptr to each Cluster
      art::Ptr<recob::Cluster> cl(clusterVecHandle, i);
      
      switch(cl->View()){
      case geo::kU :
	fView=0;
	break;
      case geo::kV :
	fView=1;
	break;
      case geo::kZ :
	fView=2;
	break;
      default :
	break;
      }// end switch on view

      Cls[fView].push_back(cl);
    }// end loop over input clusters


    for(fView = 0; fView < nplanes; ++fView){
      
      art::FindManyP<recob::Hit> fmh(Cls[fView], evt, fClusterModuleLabel);
      
      for(size_t c = 0; c < Cls[fView].size(); ++c){
	
	
	// find the hits associated with the current cluster
	std::vector< art::Ptr<recob::Hit> > ptrvs = fmh.at(c);

	double PrincDir[2], PrincValue=0;
	double TotalCharge=0;

	PerformClusterPCA( ptrvs, PrincDir, PrincValue, TotalCharge, fNormPC);
	
	fPrincDirW   = PrincDir[0];
	fPrincDirT   = PrincDir[1];
	fPrincValue  = PrincValue;
	fTotalCharge = TotalCharge;
	fNHits       = ptrvs.size();

	fTree->Fill();
	
      }// end loop over first cluster iterator
     }// end loop over planes
     
    return;

  }








// Perform PCA analysis of the hits in a cluster.
//
// This method is supplied with the vector of ptrs to hits in the cluster.
//
// It returns 2 things:
//  1. the principal direction, in [w,t] coordinates.
//  2. the principal eigenvalue, which is a measure of colinearity
//
// There is also a flag which can be specified to determine whether to normalize the PCA.
// For NormPC=true direction finding is not properly handled yet (easy to fix later)
// It is not clear yet whether true or false is more effective for cluster showeriness studies.

  void ClusterPCA::PerformClusterPCA( std::vector<art::Ptr<recob::Hit> >& HitsThisCluster, double* PrincipalDirection, double& PrincipalEigenvalue, double& TotalCharge, bool NormPC)
  {
  
  double Center[2] = {0,0};
  TotalCharge = 0;
  
  for(auto itHit = HitsThisCluster.begin(); itHit!=HitsThisCluster.end(); ++itHit)
    {
      Center[0] += (*itHit)->WireID().Wire;
      Center[1] += (*itHit)->PeakTime();
      TotalCharge += (*itHit)->Charge();
    }

  Center[0] /= float(HitsThisCluster.size());
  Center[1] /= float(HitsThisCluster.size());
      

  double WireTime[2];
  
  std::string OptionString;
  
  if(NormPC==false) 
    OptionString = "D";
  else              
    OptionString = "ND";
 
  TPrincipal  pc(2, OptionString.c_str());


  for(auto itHit = HitsThisCluster.begin(); itHit!=HitsThisCluster.end(); ++itHit)
    {
      WireTime[0] = (*itHit)->WireID().Wire - Center[0];
      WireTime[1] = (*itHit)->PeakTime()    - Center[1];
      
      pc.AddRow(WireTime);
     
    }
  
  pc.MakePrincipals();

  PrincipalEigenvalue = (*pc.GetEigenValues())[0];
  
  for(size_t n=0; n!=2; ++n)
    {
      PrincipalDirection[n]=        (*pc.GetEigenVectors())[0][n];
    }
  

  // Comment this out if you want to shut it up
  pc.Print("MSEV");
  
  pc.Clear();
  return;
}


} // end namespace

namespace cluster{

  DEFINE_ART_MODULE(ClusterPCA)
  
} // end namespace 

