//############################################################################
//### Name:        ShowerGenericTool                                       ###
//### Author:      You                                                     ###
//### Date:        13.05.19                                                ###
//### Description: Generic form of the shower tools                        ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"

namespace ShowerRecoTools {


  class ShowerGenericTool: public IShowerTool {

    public:

      ShowerGenericTool(const fhicl::ParameterSet& pset);

      ~ShowerGenericTool();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      // Function to initialise the producer i.e produces<std::vector<recob::Vertex> >();
      // commands go here.
      void InitialiseProducers() override;

      //Function to add the assoctions
      int AddAssociations(art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder) override;
  };


  ShowerGenericTool::ShowerGenericTool(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
  {
  }

  ShowerGenericTool::~ShowerGenericTool()
  {
  }

  void ShowerGenericTool::InitialiseProducers(){
  }

  int ShowerGenericTool::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){
    return 0;
  }

  int ShowerGenericTool::AddAssociations(art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){
    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerGenericTool)

