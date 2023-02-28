////////////////////////////////////////////////////////////////////////
//
// ClusterPCA class
//
// This algorithm is designed to find the principal axis and diffuseness
//  of clusters
//
////////////////////////////////////////////////////////////////////////

#include <array>
#include <string>

//Framework includes:
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft includes:
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

// Root includes
#include "TPrincipal.h"
#include "TTree.h"

namespace cluster {

  class ClusterPCA : public art::EDAnalyzer {

  public:
    explicit ClusterPCA(fhicl::ParameterSet const& pset);
    ~ClusterPCA();

  private:
    void PerformClusterPCA(const std::vector<art::Ptr<recob::Hit>>& HitsThisCluster,
                           double* PrincDirectionWT,
                           double& PrincValue,
                           double& TotalCharge,
                           bool NormPC);

    void analyze(art::Event const& evt);
    void beginJob();

    std::string fClusterModuleLabel;
    bool fNormPC;

    TTree* fTree;

    Int_t fView;
    Float_t fPrincDirW;
    Float_t fPrincDirT;
    Float_t fPrincValue;
    Float_t fTotalCharge;
    Float_t fNHits;

  }; // class ClusterPCA

}

//#endif

namespace cluster {

  //-------------------------------------------------
  ClusterPCA::ClusterPCA(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel"))
    , fNormPC(pset.get<bool>("NormPC"))
  {}

  //-------------------------------------------------
  ClusterPCA::~ClusterPCA() {}

  //-------------------------------------------------
  // Set up analysis tree
  void ClusterPCA::beginJob()
  {
    art::ServiceHandle<art::TFileService const> tfs;
    fTree = tfs->make<TTree>("PCATree", "PCATree");
    fTree->Branch("View", &fView, "View/I");
    fTree->Branch("PrincDirW", &fPrincDirW, "PrincDirW/F");
    fTree->Branch("PrincDirT", &fPrincDirT, "PrincDirT/F");
    fTree->Branch("PrincValue", &fPrincValue, "PrincValue/F");
    fTree->Branch("TotalCharge", &fTotalCharge, "TotalCharge/F");
    fTree->Branch("NHits", &fNHits, "fNHits/F");
  }

  //------------------------------------------------------------------------------------//
  void ClusterPCA::analyze(art::Event const& evt)
  {
    // Get a Handle for the input Cluster object(s).
    art::Handle<std::vector<recob::Cluster>> clusterVecHandle;
    evt.getByLabel(fClusterModuleLabel, clusterVecHandle);

    constexpr size_t nViews = 3;

    // one vector of cluster index in the original collection, per view;
    // to support all present views, replace with a dynamically growing array
    std::array<std::vector<size_t>, nViews> ClsIndex;

    // loop over the input Clusters
    for (size_t i = 0; i < clusterVecHandle->size(); ++i) {

      //get a art::Ptr to each Cluster
      art::Ptr<recob::Cluster> cl(clusterVecHandle, i);

      switch (cl->View()) {
      case geo::kU: fView = 0; break;
      case geo::kV: fView = 1; break;
      case geo::kZ: fView = 2; break;
      default:
        mf::LogError("ClusterPCA") << "Hit belongs to an unsupported view (#" << cl->View() << ")";
        break;
      } // end switch on view

      ClsIndex[fView].push_back(i);
    } // end loop over input clusters

    art::FindManyP<recob::Hit> fm(clusterVecHandle, evt, fClusterModuleLabel);
    for (fView = 0; fView < (int)nViews; ++fView) {

      const std::vector<size_t>& ClsIndices = ClsIndex[fView];

      for (size_t c = 0; c < ClsIndices.size(); ++c) {

        // find the hits associated with the current cluster
        const std::vector<art::Ptr<recob::Hit>>& ptrvs = fm.at(ClsIndices[c]);

        double PrincDir[2], PrincValue = 0;
        double TotalCharge = 0;

        PerformClusterPCA(ptrvs, PrincDir, PrincValue, TotalCharge, fNormPC);

        fPrincDirW = PrincDir[0];
        fPrincDirT = PrincDir[1];
        fPrincValue = PrincValue;
        fTotalCharge = TotalCharge;
        fNHits = ptrvs.size();

        fTree->Fill();

      } // end loop over first cluster iterator
    }   // end loop over planes

    return;

  } // ClusterPCA::analyze()

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

  void ClusterPCA::PerformClusterPCA(const std::vector<art::Ptr<recob::Hit>>& HitsThisCluster,
                                     double* PrincipalDirection,
                                     double& PrincipalEigenvalue,
                                     double& TotalCharge,
                                     bool NormPC)
  {

    double Center[2] = {0, 0};
    TotalCharge = 0;

    for (auto itHit = HitsThisCluster.begin(); itHit != HitsThisCluster.end(); ++itHit) {
      Center[0] += (*itHit)->WireID().Wire;
      Center[1] += (*itHit)->PeakTime();
      TotalCharge += (*itHit)->Integral();
    }

    Center[0] /= float(HitsThisCluster.size());
    Center[1] /= float(HitsThisCluster.size());

    double WireTime[2];

    std::string OptionString;

    if (NormPC == false)
      OptionString = "D";
    else
      OptionString = "ND";

    TPrincipal pc(2, OptionString.c_str());

    for (auto itHit = HitsThisCluster.begin(); itHit != HitsThisCluster.end(); ++itHit) {
      WireTime[0] = (*itHit)->WireID().Wire - Center[0];
      WireTime[1] = (*itHit)->PeakTime() - Center[1];

      pc.AddRow(WireTime);
    }

    pc.MakePrincipals();

    PrincipalEigenvalue = (*pc.GetEigenValues())[0];

    for (size_t n = 0; n != 2; ++n) {
      PrincipalDirection[n] = (*pc.GetEigenVectors())[0][n];
    }

    // Comment this out if you want to shut it up
    pc.Print("MSEV");

    pc.Clear();
    return;
  }

} // end namespace

namespace cluster {

  DEFINE_ART_MODULE(ClusterPCA)

} // end namespace
