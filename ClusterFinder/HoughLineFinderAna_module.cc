////////////////////////////////////////////////////////////////////////
//
// HoughLineFinder class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find lines (Houghclusters) from clusters found by DBSCAN after deconvolution and hit finding.
//  The algorithm is based on: 
//  Queisser, A. "Computing the Hough Transform", C/C++ Users Journal 21, 12 (Dec. 2003).
//  Niblack, W. and Petkovic, D. On Improving the Accuracy of the Hough Transform", Machine Vision and Applications 3, 87 (1990)  
////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include <vector>
#include <string>



extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
// ROOT includes
#include <TH1D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <iostream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/AssociationUtil.h"




class TH1F;
class TTree;
namespace cluster {
   
  class HoughLineFinderAna : public art::EDAnalyzer {
    
  public:
    
    explicit HoughLineFinderAna(fhicl::ParameterSet const& pset); 
    ~HoughLineFinderAna();
         
    void analyze(const art::Event&);
    void beginJob();
     
  private:

    std::string fHoughModuleLabel;
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fDBScanModuleLabel;       
    TTree* ftree;
    int fm_run;          // Run number
    unsigned long int fm_run_timestamp;          // Run number
    int fm_event;        // Event number
    int fm_plane;        // Plane number
    int fm_dbsize;
    int fm_clusterid;    // Cluster ID
    int fm_wirespan;    // Wire spanned by track
    int fm_sizeClusterZ;  //Number of clusters
    int fm_sizeHitZ;      //Number of Hits
    float fm_clusterslope;
    float fm_clusterintercept;
    int *fm_wireZ;
    int *fm_hitidZ;
    float *fm_mipZ;
    float *fm_drifttimeZ;
    float *fm_widthZ;
    float *fm_upadcZ;
      
  };
  
  
} // end namespace cluster

//#endif // HoughLineFinderAna_H


namespace cluster {

  HoughLineFinderAna::HoughLineFinderAna(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset) 
    , fHoughModuleLabel (pset.get< std::string >("HoughModuleLabel"))
    , fDigitModuleLabel (pset.get< std::string >("DigitModuleLabel"))
    , fHitsModuleLabel  (pset.get< std::string >("HitsModuleLabel"))
    , fDBScanModuleLabel(pset.get< std::string >("DBScanModuleLabel"))
    , fm_run(0) 
    , fm_event(0) 
    , fm_plane(0)
    , fm_dbsize(0)
    , fm_clusterid(0) 
    , fm_wirespan(0) 
    , fm_sizeClusterZ(10000) 
    , fm_sizeHitZ(10000)
    , fm_clusterslope(0) 
    , fm_clusterintercept(0)
  {
  }
  
  //-------------------------------------------------
  HoughLineFinderAna::~HoughLineFinderAna()
  {
    delete fm_hitidZ;
    delete fm_mipZ;
    delete fm_drifttimeZ;
    delete fm_widthZ;
    delete fm_upadcZ;
    delete fm_wireZ;
  }
  
  //-------------------------------------------------
  void HoughLineFinderAna::beginJob()
  {
  
  
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    ftree= tfs->make<TTree>("HoughTree","HoughTree");
    fm_hitidZ = new int[fm_sizeHitZ];
    fm_mipZ = new float[fm_sizeHitZ];
    fm_drifttimeZ = new float[fm_sizeHitZ];
    fm_widthZ = new float[fm_sizeHitZ];
    fm_upadcZ = new float[fm_sizeHitZ];
    fm_wireZ = new int[fm_sizeHitZ];
    ftree->Branch("run", &fm_run, "run/I");
    ftree->Branch("run_timestamp", &fm_run_timestamp, "run_timestamp/l"); //l is for ULong64_t
    ftree->Branch("event", &fm_event, "event/I");
    ftree->Branch("plane",&fm_plane,"plane/I");
    ftree->Branch("dbsize",&fm_dbsize,"dbsize/I");
    ftree->Branch("clusterid",&fm_clusterid,"clusterid/I");
    ftree->Branch("clusterslope",&fm_clusterslope,"clusterslope/F");
    ftree->Branch("clusterintercept",&fm_clusterintercept,"clusterintecept/F");
    ftree->Branch("wirespan",&fm_wirespan,"wirespan/I");
    ftree->Branch("numberHits",&fm_sizeHitZ,"numberHits/I");
    ftree->Branch("numberClusters",&fm_sizeClusterZ,"numberClusters/I");
    ftree->Branch("hitidZ",fm_hitidZ,"hitidZ[numberHits]/I");
    ftree->Branch("wireZ",fm_wireZ,"wireZ[numberHits]/I");
    ftree->Branch("mipZ",fm_mipZ,"mipZ[numberHits]/F");
    ftree->Branch("drifttimeZ",fm_drifttimeZ,"drifttitmeZ[numberHits]/F");
    ftree->Branch("widthZ",fm_widthZ,"widthZ[numberHits]/F");
  }
  
  
  void HoughLineFinderAna::analyze(const art::Event& evt)
  {
  
    art::Handle< std::vector<recob::Cluster> > hlfListHandle;
    evt.getByLabel(fHoughModuleLabel,hlfListHandle);
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    evt.getByLabel(fHitsModuleLabel,hitListHandle);
    art::Handle< std::vector<recob::Cluster> > dbscanListHandle;
    evt.getByLabel(fDBScanModuleLabel,dbscanListHandle);
    
    art::FindManyP<recob::Hit> fmh(dbscanListHandle, evt, fDBScanModuleLabel);
    art::FindManyP<recob::Hit> fmhhl(hlfListHandle, evt, fHoughModuleLabel);
  
    art::PtrVector<recob::Cluster> clusters;  
    art::PtrVector<recob::Cluster> dbclusters;
    //   art::PtrVector<recob::Hit> hits;// unused, as yet. EC, 5-Oct-2010.
      
    for (size_t ii = 0; ii <  hlfListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> cluster(hlfListHandle,ii);
      clusters.push_back(cluster);
    }
    
    for (size_t ii = 0; ii <  dbscanListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> dbcluster(dbscanListHandle,ii);
      dbclusters.push_back(dbcluster);
    }
    
    LOG_VERBATIM("HoughLineFinderAna") << "run    : " << evt.id().run();
    //std::cout << "subrun : " << evt.subRun() << std::endl;
    LOG_VERBATIM("HoughLineFinderAna") << "event  : " << evt.id().event();
    fm_run=evt.id().run();
    fm_event=evt.id().event();
    fm_run_timestamp=evt.time().value(); // won't cast, EC, 7-Oct-2010.
    unsigned int firstwire=0;
    unsigned int lastwire=0;
    fm_sizeClusterZ=0;
    fm_sizeHitZ=0;
    fm_dbsize=0;  
    art::ServiceHandle<geo::Geometry> geo;

    for(auto view : geo->Views()){

      fm_dbsize       = 0;
      fm_sizeClusterZ = clusters.size();
      
      for(size_t j = 0; j < dbclusters.size(); ++j) {
	if(dbclusters[j]->View() == view){
	  std::vector< art::Ptr<recob::Hit> > _dbhits = fmh.at(j);
	  fm_dbsize += _dbhits.size();
	  if(_dbhits.size() > 0) fm_plane   = _dbhits.at(0)->WireID().Plane;
	} 
      }
      
      for(size_t j = 0; j < clusters.size(); ++j) {
	if(clusters[j]->View() == view){
	  fm_clusterid=clusters[j]->ID();
	  std::vector< art::Ptr<recob::Hit> > _hits = fmhhl.at(j);
	  fm_clusterslope=(double)clusters[j]->dTdW();
	  fm_clusterintercept=(double)clusters[j]->StartPos()[1];
	  if(_hits.size()!=0){
	    fm_plane   = _hits.at(0)->WireID().Plane;
	    firstwire = _hits[0]->WireID().Wire;
	    lastwire  = _hits[_hits.size()-1]->WireID().Wire;
	    fm_wirespan = lastwire-firstwire;
	    fm_sizeHitZ = _hits.size();
  	    
	    for(unsigned int i = 0; i < _hits.size(); ++i){	     
	      
	      fm_hitidZ[i]     = i;         
	      fm_wireZ[i]      = _hits[i]->WireID().Wire;
	      fm_mipZ[i]       = (double)_hits[i]->Charge();
	      fm_drifttimeZ[i] = (double)_hits[i]->PeakTime();
	      fm_widthZ[i]     = (double)_hits[i]->EndTime()-_hits[i]->StartTime();
	      fm_upadcZ[i]     = (double)_hits[i]->Charge();
	    } 
	    
	    ftree->Fill();  
	  }
	}//end if in the correct view
      }// end loop over clusters
    }// end loop over views
  
  }
  
}// end namespace


namespace cluster{

  DEFINE_ART_MODULE(HoughLineFinderAna)
  
} // end namespace caldata

