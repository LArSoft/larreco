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
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

// Framework includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

// C++ Includes
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#include "TTree.h"
#include "TH2C.h" // ADC map
#include "TH2I.h" // PDG+vertex info map
#include "TH2F.h" // deposit map

namespace nnet	 {

  struct NUVTX
  {
  	int interaction;
  	int nupdg;
  	int cryo;
  	int tpc;
  	TVector2 position[3];
  };

  class EventImageData // full image for one plane
  {
  public:
    EventImageData(size_t weff, size_t nw, size_t nt, size_t nx, size_t ny, size_t nz, bool saveDep) :
        fVtxX(-9999), fVtxY(-9999),
        fLayerOffset(0), fTpcEffW(weff), fTpcSizeW(nw), fTpcSizeD(nt),
        fNTpcX(nx), fNTpcY(ny), fNTpcZ(nz),
        fSaveDep(saveDep)
    {
        fGeometry = &*(art::ServiceHandle<geo::Geometry>());

        size_t ntotw = (fNTpcZ - 1) * fTpcEffW + fTpcSizeW;
        fAdc.resize(ntotw, std::vector<float>(fNTpcX * fTpcSizeD, 0));
        if (saveDep) { fDeposit.resize(ntotw, std::vector<float>(fNTpcX * fTpcSizeD, 0)); }
        fPdg.resize(ntotw, std::vector<unsigned int>(fNTpcX * fTpcSizeD, 0));
    }

    void addTpc(const TrainingDataAlg & dataAlg);
    bool findCrop(size_t max_area_cut, unsigned int & w0, unsigned int & w1, unsigned int & d0, unsigned int & d1) const;

    const std::vector< std::vector<float> > & adcData(void) const { return fAdc; }
    const std::vector<float> & wireAdc(size_t widx) const { return fAdc[widx]; }

    const std::vector< std::vector<float> > & depData(void) const { return fDeposit; }
    const std::vector<float> & wireDep(size_t widx) const { return fDeposit[widx]; }

    const std::vector< std::vector<unsigned int> > & pdgData(void) const { return fPdg; }
    const std::vector<unsigned int> & wirePdg(size_t widx) const { return fPdg[widx]; }

  private:
    std::vector< std::vector<float> > fAdc, fDeposit;
    std::vector< std::vector<unsigned int> > fPdg;
    int fVtxX, fVtxY;
    size_t fLayerOffset, fTpcEffW, fTpcSizeW, fTpcSizeD;
    size_t fNTpcX, fNTpcY, fNTpcZ;
    bool fSaveDep;

    geo::GeometryCore const* fGeometry;
  };

  void EventImageData::addTpc(const TrainingDataAlg & dataAlg)
  {
    size_t tpc_z = dataAlg.TPC() / (fNTpcX * fNTpcY);
    size_t tpc_y = (dataAlg.TPC() / fNTpcX) % fNTpcY;
    size_t tpc_x = dataAlg.TPC() % fNTpcX;
    bool flipd = (fGeometry->TPC(dataAlg.TPC(), dataAlg.Cryo()).DetectDriftDirection() > 0);
    bool flipw = (dataAlg.Plane() == 1);
    float zero = dataAlg.ZeroLevel();
    
    size_t gw = tpc_z * fTpcEffW + tpc_y * fLayerOffset;
    size_t gd = tpc_x * fTpcSizeD;

    std::cout << "      flipw:" << flipw << " flipd:" << flipd << " nwires:" << dataAlg.NWires() << std::endl;
    for (size_t w = 0; w < dataAlg.NWires(); ++w)
    {
        std::vector<float> & dst = fAdc[gw + w];
        const float* src = 0;
        if (flipw)
        {
            src = dataAlg.wireData(dataAlg.NWires() - w - 1).data();
        }
        else
        {
            src = dataAlg.wireData(w).data();
        }
        if (flipd)
        {
            for (size_t d = 0; d < fTpcSizeD; ++d)
            {
                dst[gd + d] += src[fTpcSizeD - d - 1] - zero;
            }
        }
        else
        {
            for (size_t d = 0; d < fTpcSizeD; ++d)
            {
                dst[gd + d] += src[d] - zero;
            }
        }
    }
  }

  bool EventImageData::findCrop(size_t max_area_cut, unsigned int & w0, unsigned int & w1, unsigned int & d0, unsigned int & d1) const
  {
    size_t max_cut = max_area_cut / 4;
    float adcThr = 10;

    w0 = 0;
    size_t cut = 0;
    while (w0 < fAdc.size())
    {
        for (auto const d : fAdc[w0]) { if (d > adcThr) cut++; }
        if (cut < max_cut) w0++;
        else break;
    }
    w1 = fAdc.size() - 1;
    cut = 0;
    while (w1 > w0)
    {
        for (auto const d : fAdc[w1]) { if (d > adcThr) cut++; }
        if (cut < max_cut) w1--;
        else break;
    }
    w1++;

    d0 = 0;
    cut = 0;
    while (d0 < fAdc.front().size())
    {
        for (size_t i = w0; i < w1; ++i) { if (fAdc[i][d0] > adcThr) cut++; }
        if (cut < max_cut) d0++;
        else break;
    }
    d1 = fAdc.front().size() - 1;
    cut = 0;
    while (d1 > d0)
    {
        for (size_t i = w0; i < w1; ++i) { if (fAdc[i][d1] > adcThr) cut++; }
        if (cut < max_cut) d1--;
        else break;
    }
    d1++;

    unsigned int margin = 32;
    if ((w1 - w0 > 8) && (d1 - d0 > 8))
    {
        if (w0 < margin) w0 = 0;
        else w0 -= margin;

        if (w1 > fAdc.size() - margin) w1 = fAdc.size();
        else w1 += margin;
        
        if (d0 < margin) d0 = 0;
        else d0 -= margin;
        
        if (d1 > fAdc.front().size() - margin) d1 = fAdc.front().size();
        else d1 += margin;
        
        return true;
    }
    else return false;
  }

  class PointIdTrainingNuevent : public art::EDAnalyzer
  {
  public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<nnet::TrainingDataAlg::Config> TrainingDataAlg { Name("TrainingDataAlg") };
		fhicl::Atom<art::InputTag> GenModuleLabel { Name("GenModuleLabel"), Comment("Neutrino generator label.") };
		fhicl::Atom<double> FidVolCut { Name("FidVolCut"), Comment("Take events with vertex inside this cut on volume.") };
		fhicl::Sequence<int> SelectedView { Name("SelectedView"), Comment("Use selected views only, or all views if empty list.") };
		fhicl::Atom<bool> SaveDepositMap { Name("SaveDepositMap"), Comment("Save projections of the true energy depositions.") };
		fhicl::Atom<bool> SavePdgMap { Name("SavePdgMap"), Comment("Save vertex info and PDG codes map.") };
		fhicl::Atom<bool> Crop { Name("Crop"), Comment("Crop the projection to the event region plus margin.") };
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit PointIdTrainingNuevent(Parameters const& config);
    
    void beginJob() override;

    void analyze(const art::Event& event) override;

  private:
  	void ResetVars();
  
  	bool prepareEv(const art::Event& event);
  	bool InsideFidVol(TVector3 const & pvtx) const;
  	void CorrOffset(TVector3& vec, const simb::MCParticle& particle);


	TVector2 GetProjVtx(TVector3 const & vtx3d, const size_t cryo, const size_t tpc, const size_t plane) const;

    nnet::TrainingDataAlg fTrainingDataAlg;
    art::InputTag fGenieGenLabel;
	
	std::vector<int> fSelectedPlane;

	TTree *fTree;
	TTree *fTree2D;
	int fEvent;     ///< number of the event being processed
	int fRun;       ///< number of the run being processed
	int fSubRun;    ///< number of the sub-run being processed
	
	bool fCrop, fSaveDepositMap, fSavePdgMap;
		
	double fFidVolCut;
	
	NUVTX fPointid;
	int fCryo, fTpc, fPlane;
	int fPdg, fInteraction;
	float fPosX, fPosY;

	geo::GeometryCore const* fGeometry;
	detinfo::DetectorProperties const* fDetProp;
  };

  //-----------------------------------------------------------------------
  PointIdTrainingNuevent::PointIdTrainingNuevent(PointIdTrainingNuevent::Parameters const& config) : art::EDAnalyzer(config),
	fTrainingDataAlg(config().TrainingDataAlg()),
	fGenieGenLabel(config().GenModuleLabel()),
	fSelectedPlane(config().SelectedView()),
	fCrop(config().Crop()),
	fSaveDepositMap(config().SaveDepositMap()),
	fSavePdgMap(config().SavePdgMap()),
	fFidVolCut(config().FidVolCut())
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());
    fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }
  
  //-----------------------------------------------------------------------
  void PointIdTrainingNuevent::beginJob()
  {
		art::ServiceHandle<art::TFileService> tfs;

		fTree = tfs->make<TTree>("nu vertex","nu vertex tree");
		fTree->Branch("fRun", &fRun, "fRun/I");
		fTree->Branch("fEvent", &fEvent, "fEvent/I");
		fTree->Branch("fCryo", &fCryo, "fCryo/I");
		fTree->Branch("fTpc", &fTpc, "fTpc/I");
		fTree->Branch("fPdg", &fPdg, "fPdg/I");
		fTree->Branch("fInteraction", &fInteraction, "fInteraction/I");
		
		fTree2D = tfs->make<TTree>("nu vertex 2d","nu vertex 2d tree");
		fTree2D->Branch("fPlane", &fPlane, "fPlane/I");
		fTree2D->Branch("fPosX", &fPosX, "fPosX/F");
		fTree2D->Branch("fPosY", &fPosY, "fPosY/F");
  }
  
  //-----------------------------------------------------------------------
  bool PointIdTrainingNuevent::prepareEv(const art::Event& event)
  {
	art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
	if (!event.getByLabel(fGenieGenLabel, mctruthHandle)) { return false; }

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
				
				if (fSelectedPlane.empty())
				{
					if (fGeometry->TPC(tpc, cryo).HasPlane(geo::kZ))
						fSelectedPlane.push_back((int)geo::kZ);
					if (fGeometry->TPC(tpc, cryo).HasPlane(geo::kV))
						fSelectedPlane.push_back((int)geo::kV);
					if (fGeometry->TPC(tpc, cryo).HasPlane(geo::kU))
						fSelectedPlane.push_back((int)geo::kU);
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
				
				fCryo = cryo;
				fTpc = tpc;
				
				return true;	
			}
			else { return false; }
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
    ResetVars();
    
    prepareEv(event);

	std::ostringstream os;
	os << "event_" << fEvent << "_run_" << fRun << "_subrun_" << fSubRun;
	
	std::cout << "analyze " << os.str() << std::endl;

    size_t cryo = 0;
    size_t weff[3] = { 400, 400, 480 };
    size_t tpcs[3] = { 2, 2, 6 };   // FD workspace
    //size_t tpcs[3] = { 4, 1, 3 }; // ProtoDUNE

    bool goodEvent = false;
    unsigned int gd0 = 0, gd1 = 0;
	for (size_t p : fSelectedPlane)
	{
	    std::cout << "Plane: " << p << ", wires: " << fGeometry->Nwires(p, 0, cryo) << ", eff: " << weff[p] << ", zero:" << fTrainingDataAlg.ZeroLevel() << std::endl;
	    EventImageData fullimg(weff[p],
	        fGeometry->Nwires(p, 0, cryo), fDetProp->NumberTimeSamples() / fTrainingDataAlg.DriftWindow(),
	        tpcs[0], tpcs[1], tpcs[2],
	        fSaveDepositMap);

	    for (size_t t = 0; t < fGeometry->NTPC(cryo); ++t)
	    {
		    fTrainingDataAlg.setEventData(event, p, t, cryo);
		    std::cout << "   TPC: " << t << " wires: " << fTrainingDataAlg.NWires() << std::endl;
		    std::cout << "   ADC sum: " << fTrainingDataAlg.getAdcSum() << ", area: " << fTrainingDataAlg.getAdcArea() << std::endl;
		    if (fTrainingDataAlg.getAdcArea() > 150)
		    {
		        std::cout << "   add data..." << std::endl;
		        fullimg.addTpc(fTrainingDataAlg);
		    }
		    else { std::cout << "   ---> too low signal, skip tpc" << std::endl; }
		}

        std::cout << "   find crop..." << std::endl;
        unsigned int w0, w1, d0, d1;
		if (fullimg.findCrop(40, w0, w1, d0, d1))
   		{
   		    if (goodEvent)
   		    {
   		        d0 = gd0; d1 = gd1;
   		    }
   		    else
   		    {
   		        goodEvent = true;
   		        gd0 = d0; gd1 = d1;
   		    }
   			std::cout << "   crop: " << w0 << " " << w1 << " " << d0 << " " << d1 << std::endl;
   		}
	    else { std::cout << "   skip empty event" << std::endl; break; }

		std::ostringstream ss1;
   		ss1 << os.str() << "_plane_" << p; // TH2's name

		art::ServiceHandle<art::TFileService> tfs;
   		TH2C* rawHist = tfs->make<TH2C>((ss1.str() + "_raw").c_str(), "ADC", w1 - w0, w0, w1, d1 - d0, d0, d1);
   		
   		float zero = fTrainingDataAlg.ZeroLevel();
        for (size_t w = w0; w < w1; ++w)
        {
            auto const & raw = fullimg.wireAdc(w);
            for (size_t d = d0; d < d1; ++d)
            {
                rawHist->Fill(w, d, (char)(raw[d] + zero));
            }
   		}
/*
   		if (fSaveDepositMap)
   		{
       		TH2F* depHist = tfs->make<TH2F>((ss1.str() + "_deposit").c_str(), "Deposit", w1 - w0, w0, w1, d1 - d0, d0, d1);
            for (size_t w = w0; w < w1; ++w)
            {
                auto const & edep = fullimg.wireEdep(w);
                for (size_t d = d0; d < d1; ++d) { depHist->Fill(w, d, edep[d]); }
       		}
       	}

       	if (fSavePdgMap)
       	{
   			TH2I* pdgHist = tfs->make<TH2I>((ss1.str() + "_pdg").c_str(), "PDG", w1 - w0, w0, w1, d1 - d0, d0, d1);
            for (size_t w = w0; w < w1; ++w)
            {
                auto const & pdg = fullimg.wirePdg(w);
                for (size_t d = d0; d < d1; ++d) { pdgHist->Fill(w, d, pdg[d]); }
       		}
   		}
*/
   		fPdg = fPointid.nupdg;
		fInteraction = fPointid.interaction;

		fPlane = p;
			  
  		fPosX = fPointid.position[p].X();
		fPosY = fPointid.position[p].Y();

		fTree2D->Fill();
	}
	if (goodEvent) { fTree->Fill(); }

} // Raw2DRegionID::analyze()
  
  //-----------------------------------------------------------------------
  void PointIdTrainingNuevent::CorrOffset(TVector3& vec, const simb::MCParticle& particle)
  {
  	double vtx[3] = {vec.X(), vec.Y(), vec.Z()};
  	geo::TPCID tpcid = fGeometry->FindTPCAtPosition(vtx);
  	
		auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
	
		if (tpcid.isValid)
	  {
	  	float corrt0x = particle.T() * 1.e-3 * detprop->DriftVelocity();
	  	if (fGeometry->TPC(tpcid).DetectDriftDirection() == 1) { corrt0x = corrt0x*(-1); }
	  	
	  	vtx[0] = vec.X() + corrt0x;
	  }

		vec.SetX(vtx[0]);
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

	void PointIdTrainingNuevent::ResetVars()
	{
		fCryo = 0;
		fTpc = 0;
		fPlane = 0;
		fPdg = 0;
		fInteraction = 0;
		fPosX = 0.0;
		fPosY = 0.0;
	}

  DEFINE_ART_MODULE(PointIdTrainingNuevent)
}

#endif // PointIdTrainingNuevent_Module

