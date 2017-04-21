////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdTrainingNuevent
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan (Dorota.Stefan@cern.ch), May 2016
//
// Training data for PointIdAlg
//
//      We use this to dump deconv. ADC for preparation of various classifiers.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PointIdTrainingNuevent_Module
#define PointIdTrainingNuevent_Module

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

// Framework includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

// C++ Includes
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

namespace nnet	 {

  struct NUVTX
  {
  	int interaction;
  	int nupdg;
  	int cryo;
  	int tpc;
  	TVector2 position[3];
  };

  class PointIdTrainingNuevent : public art::EDAnalyzer
  {
  public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<nnet::TrainingDataAlg::Config> TrainingDataAlg {
			Name("TrainingDataAlg")
		};

		fhicl::Atom<art::InputTag> GenModuleLabel {
			Name("GenModuleLabel"),
			Comment("...")
		};

		fhicl::Atom<std::string> OutTextFilePath {
			Name("OutTextFilePath"),
			Comment("...")
		};

		fhicl::Atom<double> FidVolCut {
			Name("FidVolCut"),
			Comment("...")
		};

		fhicl::Sequence<int> SelectedView {
			Name("SelectedView"),
			Comment("...")
		};
		
		fhicl::Atom<bool> Crop {
		Name("Crop"),
		Comment("...")
		};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit PointIdTrainingNuevent(Parameters const& config);

    virtual void analyze(const art::Event& event) override;

  private:
  
  	bool PrepareEv(const art::Event& event);
  	bool InsideFidVol(TVector3 const & pvtx) const;
  	void CorrOffset(TVector3& vec, const simb::MCParticle& particle);
  	
	TVector2 GetProjVtx(TVector3 const & vtx3d, const size_t cryo, const size_t tpc, const size_t plane) const;

    nnet::TrainingDataAlg fTrainingDataAlg;

	art::InputTag fGenieGenLabel;
	std::string fOutTextFilePath;

	std::vector<int> fSelectedView;

	int fEvent;     ///< number of the event being processed
	int fRun;       ///< number of the run being processed
	int fSubRun;    ///< number of the sub-run being processed
	
	bool fCrop;
		
	double fFidVolCut;
	
	NUVTX fPointid;

	geo::GeometryCore const* fGeometry;
  };

  //-----------------------------------------------------------------------
  PointIdTrainingNuevent::PointIdTrainingNuevent(PointIdTrainingNuevent::Parameters const& config) : art::EDAnalyzer(config),
	fTrainingDataAlg(config().TrainingDataAlg()),
	fGenieGenLabel(config().GenModuleLabel()),
	fOutTextFilePath(config().OutTextFilePath()),
	fSelectedView(config().SelectedView()),
	fCrop(config().Crop()),
	fFidVolCut(config().FidVolCut())
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());
  }
  
  //-----------------------------------------------------------------------
  bool PointIdTrainingNuevent::PrepareEv(const art::Event& event)
  {
  	art::ValidHandle< std::vector<simb::MCTruth> > mctruthHandle 
  		= event.getValidHandle< std::vector<simb::MCTruth> >(fGenieGenLabel);
	
	for (auto const & mc : (*mctruthHandle))
	{
		if (mc.Origin() == simb::kBeamNeutrino)
		{	
			const TLorentzVector& pvtx = mc.GetNeutrino().Nu().Position();
			TVector3 vtx(pvtx.X(), pvtx.Y(), pvtx.Z());	
				
			CorrOffset(vtx, mc.GetNeutrino().Nu());		
		
			if (InsideFidVol(vtx)) 
			{	
				double v[3] = {vtx.X(), vtx.Y(), vtx.Z()}; 
				unsigned int cryo = fGeometry->FindCryostatAtPosition(v);
				unsigned int tpc = fGeometry->FindTPCAtPosition(v).TPC;
				
				if (fSelectedView.empty())
				{
					if (fGeometry->TPC(tpc, cryo).HasPlane(geo::kU))
						fSelectedView.push_back((int)geo::kU);
					if (fGeometry->TPC(tpc, cryo).HasPlane(geo::kV))
						fSelectedView.push_back((int)geo::kV);
					if (fGeometry->TPC(tpc, cryo).HasPlane(geo::kZ))
						fSelectedView.push_back((int)geo::kZ);
				}
				
				for (size_t i = 0; i < 3; ++i)
				{
					if (fGeometry->TPC(tpc, cryo).HasPlane(i))
					{fPointid.position[i] = GetProjVtx(vtx, cryo, tpc, i);}
					else
					{fPointid.position[i] = TVector2(0, 0);}
				}
				
				fPointid.nupdg = mc.GetNeutrino().Nu().PdgCode();
				fPointid.interaction = mc.GetNeutrino().CCNC();
				fPointid.cryo = cryo;
				fPointid.tpc = tpc;
				
				return true;	
			}
			else {return false;}
		}
	}	
	
	return false;
  }  

  //-----------------------------------------------------------------------
  void PointIdTrainingNuevent::analyze(const art::Event& event) 
  {
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    
    if (PrepareEv(event))
    {
	std::ostringstream os;
	os << "event_" << fEvent << "_run_" << fRun << "_subrun_" << fSubRun;

	std::ofstream fout_raw, fout_deposit, fout_pdg, fout_nuin;

	for (size_t v = 0; v < fSelectedView.size(); ++v)
	{
		std::ostringstream ss1;
		ss1 << fOutTextFilePath << "/raw_" << os.str()
			<< "_tpc_" << fPointid.tpc
			<< "_view_" << fSelectedView[v];

		fout_raw.open(ss1.str() + ".raw");
		fout_deposit.open(ss1.str() + ".deposit");
		fout_pdg.open(ss1.str() + ".pdg");
		fout_nuin.open(ss1.str() + ".nuin");

		fTrainingDataAlg.setEventData(event, fSelectedView[v], fPointid.tpc, fPointid.cryo);
		
		unsigned int w0, w1, d0, d1;
		if (fCrop)
		{
			if (fTrainingDataAlg.findCrop(0.004F, w0, w1, d0, d1))
      {
      	std::cout << "   crop: " << w0 << " " << w1 << " " << d0 << " " << d1 << std::endl;
      }
      else
      {
       	std::cout << "   skip empty tpc:" << fPointid.tpc << " / view:" << fSelectedView[v] << std::endl;
        continue;
      }
		}
		else
		{
			w0 = 0;
      w1 = fTrainingDataAlg.NWires();
      d0 = 0;
      d1 = fTrainingDataAlg.NScaledDrifts();
		}

		for (size_t w = w0; w < w1; ++w)
		{
			auto const & raw = fTrainingDataAlg.wireData(w);
			for (size_t d = d0; d < d1; ++d)
			{
				fout_raw << raw[d] << " ";
			}
			fout_raw << std::endl;

			auto const & edep = fTrainingDataAlg.wireEdep(w);
			for (size_t d = d0; d < d1; ++d)
			{
				fout_deposit << edep[d] << " ";
			}
			fout_deposit << std::endl;

			auto const & pdg = fTrainingDataAlg.wirePdg(w);
			for (size_t d = d0; d < d1; ++d)
			{
				fout_pdg << pdg[d] << " ";
			}
			fout_pdg << std::endl;
		}

		fout_nuin << fPointid.interaction << " " << fPointid.nupdg 
			<< " " << fPointid.position[v].X() << " " << fPointid.position[v].Y() << std::endl;

		fout_raw.close();
		fout_deposit.close();
		fout_pdg.close();
		fout_nuin.close();
	}
	
     }

  } // Raw2DRegionID::analyze()
  
  //-----------------------------------------------------------------------
  void PointIdTrainingNuevent::CorrOffset(TVector3& vec, const simb::MCParticle& particle)
  {
	auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
	float corrt0x = particle.T() * 1.e-3 * detprop->DriftVelocity();

	float px = vec.X();
	if (px > 0) px += corrt0x;
	else px -= corrt0x;
	vec.SetX(px);
  }

  //-----------------------------------------------------------------------
  bool PointIdTrainingNuevent::InsideFidVol(TVector3 const & pvtx) const
  {
	double vtx[3] = {pvtx.X(), pvtx.Y(), pvtx.Z()};
	bool inside = false;

	if (!fGeometry->FindTPCAtPosition(vtx).isValid) return false;

	geo::TPCID idtpc = fGeometry->FindTPCAtPosition(vtx);

	if (fGeometry->HasTPC(idtpc))
	{		
		const geo::TPCGeo& tpcgeo = fGeometry->GetElement(idtpc);
		double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
		double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
		double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

		//x
		double dista = fabs(minx - pvtx.X());
		double distb = fabs(pvtx.X() - maxx); 

		if ((pvtx.X() > minx) && (pvtx.X() < maxx) &&
		 	(dista > fFidVolCut) && (distb > fFidVolCut))
		{ 
			inside = true;
		}
		else { inside = false; }

		//y
		dista = fabs(maxy - pvtx.Y());
		distb = fabs(pvtx.Y() - miny);
		if (inside && (pvtx.Y() > miny) && (pvtx.Y() < maxy) &&
		 	(dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
		else inside = false;

		//z
		dista = fabs(maxz - pvtx.Z());
		distb = fabs(pvtx.Z() - minz);
		if (inside && (pvtx.Z() > minz) && (pvtx.Z() < maxz) &&
		 	(dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
		else inside = false;
	}
		
	return inside;
  }
  
  //-----------------------------------------------------------------------
  TVector2 PointIdTrainingNuevent::GetProjVtx(TVector3 const & vtx3d, const size_t cryo, const size_t tpc, const size_t plane) const
  {

	TVector2 vtx2d = pma::GetProjectionToPlane(vtx3d, plane, tpc, cryo);
	TVector2 vtxwd = pma::CmToWireDrift(vtx2d.X(), vtx2d.Y(), plane, tpc, cryo);
	
	return vtxwd;
  }


  DEFINE_ART_MODULE(PointIdTrainingNuevent)
}

#endif // PointIdTrainingNuevent_Module

