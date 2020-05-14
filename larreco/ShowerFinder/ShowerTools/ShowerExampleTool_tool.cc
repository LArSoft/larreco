//############################################################################
//### Name:        ShowerExampleTool                                       ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk          ###
//### Date:        26.06.19                                                ###
//### Description: Example form of the shower tools                        ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/RecoAlg/TRACSAlg.h"

namespace ShowerRecoTools {

  class ShowerExampleTool : public IShowerTool {

  public:
    ShowerExampleTool(const fhicl::ParameterSet& pset);

    ~ShowerExampleTool();

    //Example Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //Function to initialise the producer i.e produces<std::vector<recob::Vertex> >(); commands go here.
    void InitialiseProducers() override;

    //Function to add the assoctions
    int AddAssociations(art::Event& Event,
                        reco::shower::ShowerElementHolder& ShowerEleHolder) override;

    //prehaps you want a fcl parameter.
    art::InputTag fPFParticleModuleLabel;
  };

  ShowerExampleTool::ShowerExampleTool(const fhicl::ParameterSet& pset)
    : //Setup the algs and others here
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel", ""))
  {}

  ShowerExampleTool::~ShowerExampleTool() {}

  void
  ShowerExampleTool::InitialiseProducers()
  {

    //Do you create something and you want to save it the event. Initialsie here. For every event with have a vector of showers so each one has a vertex. This is what we are saving. Make sure to use the name "myvertex" later down the line.
    InitialiseProduct<std::vector<recob::Vertex>>("myvertex");

    //We can also do associations
    InitialiseProduct<art::Assns<recob::Shower, recob::Vertex>>("myvertexassan");
  }

  int
  ShowerExampleTool::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                      art::Event& Event,
                                      reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //In here calculate a shower or element (or multiple). It can be something used to create the recob::shower i.e. the direction. These have specific names so be careful to make these correctly. Alternative you can create something completely new e.g. recob::Vertex and add it the shower element holder

    //Now we are calculating the property of the shower like pfparticle. You have access to everything in the event. Maybe you want the vertex.
    art::Handle<std::vector<recob::Vertex>> vtxHandle;
    std::vector<art::Ptr<recob::Vertex>> vertices;
    if (Event.getByLabel(fPFParticleModuleLabel, vtxHandle))
      art::fill_ptr_vector(vertices, vtxHandle);
    else {
      throw cet::exception("ShowerExampleTool")
        << "Could not get the pandora vertices. Something is not configured correctly. Please give "
           "the correct pandora module label. Stopping";
      return 1;
    }

    //Remember the module goes through the tools and if you want to (fcl param) it will loop over them twice. You can check to see if a element has been set with a specific name:
    bool shower_direction_set = ShowerEleHolder.CheckElement("ShowerDirection");

    TVector3 ShowerDirection = {-999, -999, -999};

    //Then you can go and get that element if you want to use it and fill it in for you.
    if (shower_direction_set) { ShowerEleHolder.GetElement("ShowerDirection", ShowerDirection); }

    //Do some crazy physics - Some legacy code in here for ease.
    art::Ptr<recob::Vertex> proposed_vertex = vertices[0];
    double xyz[3] = {-999, -999, -999};
    proposed_vertex->XYZ(xyz);

    if (ShowerDirection.X() < 0) {
      xyz[0] = -xyz[0];
      xyz[1] = -xyz[1];
      xyz[2] = -xyz[2];
    }
    recob::Vertex new_vertex = recob::Vertex(xyz);
    TVector3 recobshower_vertex = {xyz[0], xyz[1], xyz[2]};
    TVector3 recobshower_err = {xyz[0] * 0.1, xyz[1] * 0.1, xyz[2] * 0.1};
    //You can set elements of the recob::shower just choose the right name (you can acess them later). You can give the property an error anf this must be done the for standard recob::shower properties; The standard is to access the name via a fcl file.
    ShowerEleHolder.SetElement(recobshower_vertex, recobshower_err, "ShowerStartPosition");

    //You can also set the same element with a different name so that you can compare downstream two tools.
    //The standard is to actually define the name in fcl.
    ShowerEleHolder.SetElement(
      recobshower_vertex, recobshower_err, "ShowerExampleTool_ShowerStartPosition");

    //Or you can set one of the save elements
    ShowerEleHolder.SetElement(new_vertex, "myvertex");

    //Or a new unsave one.
    std::vector<double> xyz_vec = {xyz[0], xyz[1], xyz[2]};
    ShowerEleHolder.SetElement(xyz_vec, "xyz");

    //If you want to check if your element was actually made before the shower is made you can set a bool. If partial showers is turned off then the shower will not be made if this element is not filled. Properties i.e. elements with errors i.e. ShowerStartPosition  will not be checked. There is no way to store properties in the Event, only products are stored. You can make your own class which holds the error. The defualt is not to check the element. The recob::shower properties are checked however.
    ShowerEleHolder.SetElement(xyz_vec, "xyz", true);

    //You can see if an element will be checked before the shower is save with
    bool will_be_checked = ShowerEleHolder.CheckElementTag("xyz");

    if (will_be_checked) { std::cout << "Element checked at save time" << std::endl; }

    //You can also changed the tag.
    ShowerEleHolder.SetElementTag("xyz", false);

    //Note: Elements that are actually saved because you defined them in InitialiseProducers will be checked regardless. We don't want you saving nothign now.

    //You can also get the shower number that you are current one (the first shower number is 0).
    int showernum = ShowerEleHolder.GetShowerNumber();
    std::cout << "You on are shower: " << showernum << std::endl;

    //You can also read out what ptr are set and what elements are set:.
    PrintPtrs();
    PrintPtr("myvertex");
    ShowerEleHolder.PrintElements();
    ShowerEleHolder.PrintElement("myvertex");

    //Remember to add make a new fcl parmas list for your new tool. For examles see showertools.fcl. And remember to add it the the list in the module fcl params list.

    return 0;
  }

  int
  ShowerExampleTool::AddAssociations(art::Event& Event,
                                     reco::shower::ShowerElementHolder& ShowerEleHolder)
  {
    //Here you add elements to associations defined. You can get the art::Ptrs by  GetProducedElementPtr<T>. Then you can add single like a usally association using AddSingle<assn<T>. Assn below.

    //First check the element has been set
    if (!ShowerEleHolder.CheckElement("myvertex")) {
      mf::LogError("ShowerExampleTooAddAssn") << "vertex not set." << std::endl;
      return 1;
    }

    //Then you can get the size of the vector which the unique ptr hold so that you can do associations. If you are comfortable in the fact that your element will always be made when a shower is made you don't need to to do this you can just get the art ptr as:      const art::Ptr<recob::Vertex> vertexptr = GetProducedElementPtr<recob::Vertex>("myvertex", ShowerEleHolder);. Note doing this when you allow partial showers to be set can screw up the assocation for the partial shower.
    int ptrsize = GetVectorPtrSize("myvertex");

    const art::Ptr<recob::Vertex> vertexptr =
      GetProducedElementPtr<recob::Vertex>("myvertex", ShowerEleHolder, ptrsize);
    const art::Ptr<recob::Shower> showerptr =
      GetProducedElementPtr<recob::Shower>("shower", ShowerEleHolder);
    AddSingle<art::Assns<recob::Shower, recob::Vertex>>(showerptr, vertexptr, "myvertexassan");

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerExampleTool)
