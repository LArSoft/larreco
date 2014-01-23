////////////////////////////////////////////////////////////////////////
/// \file  ShowerSelectorFilter_module.cc
/// \brief Extract showers based on cluster Parameters information
///
/// \version $Id: ShowerSelectorFilter_module.cc,v  $
/// \author andrzej.szelc@yale.edu
////////////////////////////////////////////////////////////////////////
// This class solves the following problem:
//
// Create a cluster saving its slope angle

// Framework includes

  // Framework includes
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" //#include "ClusterFinder/ShowerSelector.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "CLHEP/Random/JamesRandom.h"


#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoAlg/ClusterParamsAlg.h" 
#include "Utilities/AssociationUtil.h"
#include "Utilities/SeedCreator.h"



namespace cluster {

  
class ShowerSelectorFilter : public art::EDFilter {
   public:
     explicit ShowerSelectorFilter(fhicl::ParameterSet const& pset);
    virtual ~ShowerSelectorFilter() { }
    virtual bool filter(art::Event& e);
    void    reconfigure(fhicl::ParameterSet const& pset);

  private:

    // Control parameter: 1 to select odd numbered events; 
    //                    else select even numbered event.
    int fkeepOddOrEven;
    ClusterParamsAlg fCParAlg;
    std::string fClusterModuleLabel;
  };

   ShowerSelectorFilter::ShowerSelectorFilter(fhicl::ParameterSet const& pset)
   :fCParAlg(pset.get< fhicl::ParameterSet >("ClusterParamsAlg"),pset.get< std::string >("module_type"))
  {
    this->reconfigure(pset);
    
       // Create random number engine needed for PPHT
  createEngine(SeedCreator::CreateRandomNumberSeed(),"HepJamesRandom");
  }

  void ShowerSelectorFilter::reconfigure(fhicl::ParameterSet const& pset)
  {
   // fkeepOddOrEven = pset.get<int>("keepOddOrEven");
    fClusterModuleLabel= pset.get<std::string>("ClusterModuleLabel");
    fCParAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterParamsAlg"));
  }

  bool ShowerSelectorFilter::filter(art::Event& e)
  {

    
    double wire_start=0,wire_end=0,time_start=0,time_end=0;
    double lineslope=0,lineintercept=0,goodness=0;
    
    std::cout << " cluster module label: " << fClusterModuleLabel << std::endl;
    
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    e.getByLabel(fClusterModuleLabel,clusterListHandle);
    art::Handle< std::vector<art::PtrVector < recob::Cluster> > > clusterAssociationHandle;
    art::FindManyP<recob::Hit> fmh(clusterListHandle, e, fClusterModuleLabel);
    
   int nofshowerclusters=0; 
    
   for(unsigned int iClust = 0; iClust < clusterListHandle->size(); iClust++){

    art::Ptr<recob::Cluster> cl(clusterListHandle, iClust);
    std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iClust);
   
    std::cout << std::endl << "+++++ iClust: " << iClust << ", plane: " << hitlist[0]->WireID().Plane << std::endl;
    std::cout << " hitlist size: " << hitlist.size() << std::endl;
    if(hitlist.size()<10)
    {
     std::cout << "continuing due to small number of hits" << std::endl;
      continue;    // can't do much with less
      
    }
    if(fCParAlg.Find2DAxisRough(lineslope,lineintercept,goodness,hitlist)==-1)
      continue;
    std::cout << " lineslope, intercept " << lineslope << " "<< lineintercept << std::endl;

    
  //  fCParAlg.Find2DAxisRoughHighCharge(lineslope,lineintercept,goodness,hitlist);
  //  std::cout << "High Charge parameters lineslope, intercept " << lineslope << " "<< lineintercept << std::endl;

    
    if(fCParAlg.Find2DStartPointsHighCharge( hitlist,wire_start,time_start,wire_end,time_end)==-1)
    {
      std::cout << " exiting on StartPointsHighCharge " << std::endl;
      continue;
      
    }
    std::cout << " high charge basic start points: (" << wire_start<<","<<time_start<<"), end: ( "<<wire_end << ","<<time_end <<")" <<  std::endl;

    double wstn=0,tstn=0,wendn=0,tendn=0;
    fCParAlg.FindTrunk(hitlist,wstn,tstn,wendn,tendn,lineslope,lineintercept);
    std::cout << " trunk start points: (" << wstn<<","<<tstn<<"), end: ( "<<wendn << ","<<tendn <<")" <<  std::endl;
    
    double HiBin,LowBin,invHiBin,invLowBin;
    fCParAlg.FindDirectionWeights(lineslope,wstn,tstn,wendn,tendn,hitlist,HiBin,LowBin,invHiBin,invLowBin); 
    std::cout << " Direction weights:  norm: " << HiBin << " " << LowBin << " Inv: " << invHiBin << " " << invLowBin << std::endl;
    
    int locDirection=1;
    fCParAlg.RefineStartPointsHough(hitlist, wstn,tstn,wendn,tendn,locDirection); 
    
       
    if(fCParAlg.isShower(lineslope,wstn,tstn,wendn,tendn,hitlist))
      nofshowerclusters++;
    
//     int fDirection=fCParAlg.DecideClusterDirection(hitlist,lineslope,wstn,tstn,wendn,tendn);
//      std::cout << " direction start points: (" << wstn<<","<<tstn<<"), end: ( "<<wendn << ","<<tendn <<")" << "Direction: " << fDirection << std::endl;
//     wire_start=wstn;
//     time_start=tstn;
//     wire_end=wendn;
//     time_end=tendn;
//     fCParAlg.RefineStartPointsHough(hitlist, wire_start,time_start,wire_end,time_end,fDirection); 
//      std::cout << " Hough line refine start points: (" << wire_start<<","<<time_start<<"), end: ( "<<wire_end << ","<<time_end <<")" << "Direction: " << fDirection << std::endl; 
    
   }
   
   //divide the shower weight by length!!!
   
    // Get event number from the event.
    int event = e.id().event();

    std::cout << "+++ inside of filter: event " << event  << " n of shower clusters: " << nofshowerclusters << std::endl;
    // Always keep event 3.
    if(nofshowerclusters>=2) 
      return true;
    else
      return false;

  }

 
  
  DEFINE_ART_MODULE(ShowerSelectorFilter)

}
