////////////////////////////////////////////////////////////////////////
// Class:       MultiEMShowers
// Module Type: analyzer
// File:        MultiEMShowers_module.cc
// Author: dorota.stefan@cern.ch robert.sulej@cern.ch
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"

#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Shower.h"
#include "MCCheater/BackTracker.h"

#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include <memory>

#include "DirOfGamma/DirOfGamma.h"

// ROOT includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMathBase.h"

namespace ems {
	class MCinfo;
  	class MultiEMShowers;
}

class ems::MCinfo
{
	public:
	MCinfo(const art::Event& evt);
	void Info(const art::Event& evt);

	int GetNgammas(void) const { return fNgammas; }

	double GetMompi0(void) const { return fMompi0; }
	double GetMomGamma1(void) const { return fGammamom1; }
	double GetMomGamma2(void) const { return fGammamom2; }

	double GetCosine(void) { return fCosine; }

	TVector3 GetPrimary(void) const & { return fPrimary; }
	TVector3 GetPospi0(void) const & { return fPi0pos; }
	TVector3 GetPosgamma1(void) const & { return fConvgamma1; }
	TVector3 GetPosgamma2(void) const & { return fConvgamma2; }

	TVector3 GetDirgamma1(void) const & { return fDirgamma1; }
	TVector3 GetDirgamma2(void) const & { return fDirgamma2; }

	bool IsInside1(void) const & { return fInside1; }
	bool IsInside2(void) const & { return fInside2; }

	bool IsCompton(void) const & { return fCompton; }

	private:
	bool insideFidVol(const TLorentzVector& pvtx);
	
	double fFidVolCut;

	int fNgammas;

	double fMompi0;
	double fGammamom1;
	bool fInside1;
	double fGammamom2;
	bool fInside2;
	
	double fCosine;

	bool fCompton;

	TVector3 fPrimary;
	TVector3 fPi0pos;
	TVector3 fConvgamma1;
	TVector3 fConvgamma2;
	TVector3 fDirgamma1;
	TVector3 fDirgamma2;
};

ems::MCinfo::MCinfo(const art::Event& evt) :
fFidVolCut(2.0)
{
	Info(evt);
}

void ems::MCinfo::Info(const art::Event& evt)
{
	fMompi0 = 0.0; fPi0pos.SetXYZ(0,0,0); 
	fNgammas = 0;
	fCosine = 0.0;
	fInside1 = false; fInside2 = false;
	fCompton = false;

	fGammamom1 = 0.0; fGammamom2 = 0.0;
	fConvgamma1.SetXYZ(0,0,0); fConvgamma2.SetXYZ(0,0,0); 
	fDirgamma1.SetXYZ(0,0,0); fDirgamma2.SetXYZ(0,0,0);

	art::ServiceHandle< cheat::BackTracker > bt;
	const sim::ParticleList& plist = bt->ParticleList();
	for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar)
	{
		simb::MCParticle* particle = ipar->second;
		if (particle->Process() == "primary") 
		{
			TLorentzVector posvec = particle->Position();
			TVector3 pose(posvec.X(), posvec.Y(), posvec.Z());
			fPrimary = pose;
		}

		if ((particle->Process() == "primary") && (particle->PdgCode() == 111))
		{
			fMompi0 = particle->P();
			
			TLorentzVector posvec3 = particle->Position();
			TVector3 pospi0(posvec3.X(), posvec3.Y(), posvec3.Z());
			fPi0pos =  pospi0;

			if ((particle->NumberDaughters() == 2) &&
			    (bt->TrackIDToParticle(particle->Daughter(0))->PdgCode() == 22) &&
			    (bt->TrackIDToParticle(particle->Daughter(1))->PdgCode() == 22)) // pi0
			{
				fNgammas = particle->NumberDaughters();
				TLorentzVector mom1 = bt->TrackIDToParticle(particle->Daughter(0))->Momentum();
				TLorentzVector mom2 = bt->TrackIDToParticle(particle->Daughter(1))->Momentum();

				const simb::MCParticle* daughter1 = bt->TrackIDToParticle(particle->Daughter(0));
				// compton process
				if (daughter1->EndProcess() == "phot") fCompton = true;
				
				const simb::MCParticle* daughter2 = bt->TrackIDToParticle(particle->Daughter(1));
				if (daughter2->EndProcess() == "phot") fCompton = true;

				TVector3 mom1vec3(mom1.Px(), mom1.Py(), mom1.Pz());
				fGammamom1 = bt->TrackIDToParticle(particle->Daughter(0))->P();
				TVector3 mom2vec3(mom2.Px(), mom2.Py(), mom2.Pz());
				fGammamom2 = bt->TrackIDToParticle(particle->Daughter(1))->P();

				TLorentzVector pos1 = bt->TrackIDToParticle(particle->Daughter(0))->EndPosition();
				TLorentzVector pos2 = bt->TrackIDToParticle(particle->Daughter(1))->EndPosition();
				
				if (insideFidVol(pos1)) fInside1 = true; 
				if (insideFidVol(pos2)) fInside2 = true;
				
				TVector3 pos1vec3(pos1.X(), pos1.Y(), pos1.Z());
				fConvgamma1 = pos1vec3;
				TVector3 pos2vec3(pos2.X(), pos2.Y(), pos2.Z());
				fConvgamma2 = pos2vec3;

				TVector3 vecnorm1 = mom1vec3 * (1.0 / mom1vec3.Mag());
				fDirgamma1 = vecnorm1;
				TVector3 vecnorm2 = mom2vec3 * (1.0 / mom2vec3.Mag());
				fDirgamma2 = vecnorm2;
		
				fCosine = fDirgamma1 * fDirgamma2;
				break;
				
			}
			else
			{
				fNgammas = particle->NumberDaughters();
			}
		}
		else 
		{
			
		}	
	}
}

bool ems::MCinfo::insideFidVol(const TLorentzVector& pvtx) 
{
	art::ServiceHandle<geo::Geometry> geom;
	double vtx[3] = {pvtx.X(), pvtx.Y(), pvtx.Z()};
	bool inside = false;

	geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);
	if (geom->HasTPC(idtpc))
	{
		
		const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
		double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
		double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
		double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

		for (size_t c = 0; c < geom->Ncryostats(); c++)
		{
			const geo::CryostatGeo& cryostat = geom->Cryostat(c);
			for (size_t t = 0; t < cryostat.NTPC(); t++)
			{
				const geo::TPCGeo& tpcg = cryostat.TPC(t);
				if (tpcg.MinX() < minx) minx = tpcg.MinX();
				if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX(); 
				if (tpcg.MinY() < miny) miny = tpcg.MinY();
				if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
				if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
				if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
			}
		}	

		//x
		double dista = fabs(minx - pvtx.X());
		double distb = fabs(pvtx.X() - maxx); 
		if ((pvtx.X() > minx) && (pvtx.X() < maxx) &&
		 	(dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
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

class ems::MultiEMShowers : public art::EDAnalyzer {
public:
 	explicit MultiEMShowers(fhicl::ParameterSet const & p);

  	MultiEMShowers(MultiEMShowers const &) = delete;
 	MultiEMShowers(MultiEMShowers &&) = delete;
  	MultiEMShowers & operator = (MultiEMShowers const &) = delete;
  	MultiEMShowers & operator = (MultiEMShowers &&) = delete;

	void beginJob() override;
	void endJob() override;
  
  	void analyze(art::Event const & e) override;

	void reconfigure(fhicl::ParameterSet const& p);

private:
	bool convCluster(art::Event const & evt);
	double getMinDist(std::vector< art::Ptr<recob::Hit> > const & v, 
								TVector3 const & convmc, 
								size_t view, size_t tpc, size_t cryo);
	int fConvGood;
	int fConvWrong;
	int fConvBothGood;
	int fGammasInside;

	// ROOT
	TTree* fEvTree; 
	int fEvNumber;
	int fNGroups;

	// mc 
	double fPi0mom;
	double fGmom1;
	double fGmom2;
	double fMcth;
	int fNgammas;
	int fEvFidVol;
	int fEvComp;
	int fEvGMomCut;	
	int fEvInput;
	TVector3 fGdir1;
	TVector3 fGdir2;
	TVector3 fPrimary;

	//reco
	int fEvReco;
	int fEv2Groups;
	int fEv2Good;
	int fCountph;
	int fCountreco;
	
	TTree* fShTree; 	
	TTree* fRecoTree;
	double fStartX; double fStartY; double fStartZ;
	double fDedxZ; double fDedxV; double fDedxU;
	double fMCrecovtx; double fMCrecoTh;
	double fMCrecovtxgood; double fMCrecoThgood;
	double fRecth; double fRecthgood;
	double fDistConvrecomc1; double fDistConvrecomc2;
	double fGdirmcreco1; double fGdirmcreco2;
	double fGdirmcreco1good; double fGdirmcreco2good;

  	std::string fHitsModuleLabel;
	std::string fCluModuleLabel;
	std::string fTrk3DModuleLabel;
	std::string fVtxModuleLabel;
	std::string fShsModuleLabel;
};


ems::MultiEMShowers::MultiEMShowers(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  
{
	fConvGood = 0;
	fConvWrong = 0;
	fConvBothGood = 0;
	fEvFidVol = 0;
	fEvComp = 0;
	fEvGMomCut = 0;
	fEvReco = 0;
	fEvInput = 0;
	fEv2Groups = 0;
	fEv2Good = 0;
	reconfigure(p);
}

void ems::MultiEMShowers::reconfigure(fhicl::ParameterSet const& p)
{
  fHitsModuleLabel = p.get< std::string >("HitsModuleLabel");
  fCluModuleLabel = p.get< std::string >("ClustersModuleLabel");
  fTrk3DModuleLabel = p.get< std::string >("Trk3DModuleLabel");
  fVtxModuleLabel = p.get< std::string >("VtxModuleLabel");
  fShsModuleLabel = p.get< std::string >("ShsModuleLabel");
  
  return;
}

void ems::MultiEMShowers::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;

	fEvTree = tfs->make<TTree>("MultiShowers", "showers3d");
	fEvTree->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fEvTree->Branch("fNGroups", &fNGroups, "fNGroups/I");
	fEvTree->Branch("fPi0mom", &fPi0mom, "fPi0mom/D");
	fEvTree->Branch("fNgammas", &fNgammas, "fNgammas/I");
	fEvTree->Branch("fGmom1", &fGmom1, "fGmom1/D");
	fEvTree->Branch("fGmom2", &fGmom2, "fGmom2/D");
	fEvTree->Branch("fMcth", &fMcth, "fMcth/D");
	fEvTree->Branch("fRecth", &fRecth, "fRecth/D");
	fEvTree->Branch("fMCrecovtx", &fMCrecovtx, "fMCrecovtx/D");
	fEvTree->Branch("fMCrecoTh", &fMCrecoTh, "fMCrecoTh/D");
	fEvTree->Branch("fConvGood", &fConvGood, "fConvGood/I");
	fEvTree->Branch("fConvWrong", &fConvWrong, "fConvWrong/I");
	fEvTree->Branch("fConvBothGood", &fConvBothGood, "fConvBothGood/I");
	fEvTree->Branch("fGammasInside", &fGammasInside, "fGammasInside/I");
	fEvTree->Branch("fCountph", &fCountph, "fCountph/I");
	fEvTree->Branch("fCountreco", &fCountreco, "fCountreco/I");

	fRecoTree = tfs->make<TTree>("Cascades", "conv points");
	fRecoTree->Branch("fRecth", &fRecth, "fRecth/D");
	fRecoTree->Branch("fMCrecovtx", &fMCrecovtx, "fMCrecovtx/D");
	fRecoTree->Branch("fMCrecoTh", &fMCrecoTh, "fMCrecoTh/D");
	fRecoTree->Branch("fRecthgood", &fRecthgood, "fRecthgood/D");
	fRecoTree->Branch("fMCrecovtxgood", &fMCrecovtxgood, "fMCrecovtxgood/D");
	fRecoTree->Branch("fMCrecoThgood", &fMCrecoThgood, "fMCrecoThgood/D");
	fRecoTree->Branch("fGdirmcreco1", &fGdirmcreco1, "fGdirmcreco1/D");
	fRecoTree->Branch("fGdirmcreco2", &fGdirmcreco2, "fGdirmcreco2/D");
	fRecoTree->Branch("fGdirmcreco1good", &fGdirmcreco1good, "fGdirmcreco1good/D");
	fRecoTree->Branch("fGdirmcreco2good", &fGdirmcreco2good, "fGdirmcreco2good/D");

	fShTree = tfs->make<TTree>("Shower", "conv point");
	
	fShTree->Branch("fStartX", &fStartX, "fStartX/D");
	fShTree->Branch("fStartY", &fStartY, "fStartY/D");
	fShTree->Branch("fStartZ", &fStartZ, "fStartZ/D");
	fShTree->Branch("fDedxZ", &fDedxZ, "fDedxZ/D");
	fShTree->Branch("fDedxV", &fDedxV, "fDedxV/D");
	fShTree->Branch("fDedxU", &fDedxU, "fDedxU/D");
}

void ems::MultiEMShowers::endJob()
{
	std::cout << "******************** fEvFidVol =  " << fEvFidVol << std::endl;
	std::cout << "******************** fEvGMomCut = " << fEvGMomCut << std::endl;
	std::cout << "******************** fEvComp =    " << fEvComp << std::endl;
	std::cout << "******************** fEvReco =    " << fEvReco << std::endl;
	std::cout << "******************** fEvInput =   " << fEvInput << std::endl;
	std::cout << "******************** fEv2Groups = " << fEv2Groups << std::endl;
	std::cout << "******************** fEv2Good =   " << fEv2Good << std::endl;
	if (fEvInput)
	std::cout << "******************** reco %  =    " << double(fEvReco)/double(fEvInput) << std::endl; 
}

void ems::MultiEMShowers::analyze(art::Event const & e)
{
	fEvNumber = e.id().event();
	fNGroups = 0; 
	fStartX = 0.0; fStartY = 0.0; fStartZ = 0.0;
	fPi0mom = 0.0; fNgammas = 0;
	fDistConvrecomc1 = 0.0; fDistConvrecomc2 = 0.0;
	fMCrecovtx = -400.0;	
	fMCrecovtxgood = -400.0;
	fRecth = -400.0;
	fRecthgood = -400.0;
	fMCrecoTh = -400.0;
	fMCrecoThgood = -400.0;
	fGammasInside = 0;
	fCountph = 0;
	fCountreco = 0;
	fGdirmcreco1 = 0.0;
	fGdirmcreco2 = 0.0;
	fGdirmcreco1good = 0.0;
	fGdirmcreco2good = 0.0;
	fDedxZ = 0.0; fDedxV = 0.0; fDedxU = 0.0;

	ems::MCinfo mc(e);
	fPrimary = mc.GetPrimary();
	fPi0mom = mc.GetMompi0();
	fGmom1 = mc.GetMomGamma1();
	fGmom2 = mc.GetMomGamma2();
	fGdir1 = mc.GetDirgamma1();
	fGdir2 = mc.GetDirgamma2();
	fNgammas = mc.GetNgammas();
	TVector3 pospi0 = mc.GetPospi0();

	double cosinemc = mc.GetCosine();
	fMcth = 180.0F * (std::acos(cosinemc)) / TMath::Pi();
	TVector3 convp[2]; 
	convp[0] = mc.GetPosgamma1();
	convp[1] = mc.GetPosgamma2();
	const double maxdist = 2.0; //cm

	// check whether two photons are inside fid vol
	if (mc.IsInside1() && mc.IsInside2()) 
	{
		fGammasInside = 1;
		fEvFidVol++;
	}

	if ((fGmom1 > 0.1) && (fGmom2 > 0.1)) fEvGMomCut++;

	if (mc.IsCompton()) fEvComp++;

	art::Handle< std::vector<recob::Shower> > shsListHandle;
 	art::Handle< std::vector<recob::Track> > trkListHandle;
	art::Handle< std::vector<recob::Vertex> > vtxListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	art::Handle< std::vector<recob::Hit> > hitListHandle;
	if (e.getByLabel(fShsModuleLabel, shsListHandle) &&
		e.getByLabel(fTrk3DModuleLabel, trkListHandle) &&
	    e.getByLabel(fVtxModuleLabel, vtxListHandle) &&
	    e.getByLabel(fCluModuleLabel, cluListHandle) &&
	    e.getByLabel(fHitsModuleLabel, hitListHandle))
	{
		art::FindManyP< recob::Cluster > cluFromShs(shsListHandle, e, fShsModuleLabel);
		art::FindManyP< recob::Cluster > cluFromTrk(trkListHandle, e, fTrk3DModuleLabel);
		art::FindManyP< recob::Vertex > vtxFromTrk(trkListHandle, e, fVtxModuleLabel);
		art::FindManyP< recob::Hit > hitFromClu(cluListHandle, e, fCluModuleLabel);

		fNGroups = shsListHandle->size();	

		fCountph = 0;
		if (fNgammas == 2) // pi0
		{
			int idph = -1; 
			for (size_t s = 0; s < shsListHandle->size(); ++s)
			{
				const recob::Shower& sh = (*shsListHandle)[s];
				double mindist = maxdist; bool found = false; 
			
				for (int i = 0; i < fNgammas; ++i)
				{
					double dist = sqrt(pma::Dist2(sh.ShowerStart(), convp[i]));
					if ((dist < mindist) && (idph != i)) 
					{ mindist =  dist; idph = i; found = true; }
				}
				if (found) { fConvGood++; fCountph++; }
				else { fConvWrong++; }
			}
			if (fCountph == 2) fConvBothGood++;
		
			// plot a few variables if there are 2 showers
			if (fCountph == 2)
				for (size_t s = 0; s < shsListHandle->size(); ++s)	
				{
					const recob::Shower& sh = (*shsListHandle)[s];
					TVector3 pos = sh.ShowerStart(); 
					fStartX = pos.X(); fStartY = pos.Y(); fStartZ = pos.Z();
					std::vector<double> vecdedx = sh.dEdx();
					
					if (vecdedx.size() == 3)
					{
						fDedxZ = vecdedx[0]; fDedxV = vecdedx[1]; fDedxU = vecdedx[2];
					}
					
					fShTree->Fill();	
				}
		}
		else  // other than pi0
		{
			for (size_t s = 0; s < shsListHandle->size(); ++s)
			{
				const recob::Shower& sh = (*shsListHandle)[s];
				double mindist = maxdist; 
			
				double dist = sqrt(pma::Dist2(sh.ShowerStart(), fPrimary));
				if (dist < mindist) 
				{ 
					TVector3 pos = sh.ShowerStart(); 
					fStartX = pos.X(); fStartY = pos.Y(); fStartZ = pos.Z();
					std::vector<double> vecdedx = sh.dEdx();
					if (vecdedx.size() == 3)
					{
						fDedxZ = vecdedx[0]; fDedxV = vecdedx[1]; fDedxU = vecdedx[2];
					}
				}
					
				fShTree->Fill();	
			}
		}
		// compute the crossing point

		//cut from mc and clusters
	
		if (mc.IsInside1() && mc.IsInside2() && (fGmom1 > 0.1) && (fGmom2 > 0.1) && (!mc.IsCompton()) && convCluster(e))
		{
			fCountreco = 1;
			if (fNGroups == 2) fEv2Groups++;	
			if ((fNGroups == 2) && (fCountph == 2)) fEv2Good++;
			// cut from reco
			//if (countph == 2)
			if (fNGroups == 2)
			{
				std::vector< std::pair<TVector3, TVector3> > lines;
				const recob::Shower& sh1 = (*shsListHandle)[0];
				const recob::Shower& sh2 = (*shsListHandle)[1];

				std::pair<TVector3, TVector3> frontback1(sh1.ShowerStart(), sh1.ShowerStart() + sh1.Direction());
				std::pair<TVector3, TVector3> frontback2(sh2.ShowerStart(), sh2.ShowerStart() + sh2.Direction());
				lines.push_back(frontback1); lines.push_back(frontback2);

				TVector3 result;
				pma::SolveLeastSquares3D(lines, result); // mse.

				double dist1_0 = pma::Dist2(result, sh1.ShowerStart());
				double dist2_0 = pma::Dist2(result, sh1.ShowerStart() + sh1.Direction());
				double dist1_1 = pma::Dist2(result, sh2.ShowerStart());
				double dist2_1 = pma::Dist2(result, sh2.ShowerStart() + sh2.Direction());
				if ((dist1_0 > dist2_0) || (dist1_1 > dist2_1)) {}
				else
				{
					fMCrecovtx = std::sqrt(pma::Dist2(pospi0, result));

					if (fCountph == 2) fMCrecovtxgood = fMCrecovtx;

					double cosine_reco = sh1.Direction() * sh2.Direction();
					fRecth = 180.0F * (std::acos(cosine_reco)) / TMath::Pi();

					fGdirmcreco1 = fGdir1 * sh1.Direction();
					fGdirmcreco2 = fGdir2 * sh2.Direction();
					if (fCountph == 2)
					{
						fGdirmcreco1good = fGdirmcreco1;
						fGdirmcreco2good = fGdirmcreco2;
					}

					if (fCountph == 2) fRecthgood = fRecth;

					fMCrecoTh = fRecth - fMcth;

					if (fCountph == 2) fMCrecoThgood = fMCrecoTh;					

					fEvReco++;
					fRecoTree->Fill();
				}
			}
			fEvInput++;
			//fRecoTree->Fill();
		}
	}

	fEvTree->Fill();		
}

// true if there are clusters corresponding to mc conversion points
bool ems::MultiEMShowers::convCluster(art::Event const & evt)
{
	ems::MCinfo mc(evt);
	TVector3 convp[2]; 
	convp[0] = mc.GetPosgamma1();
	convp[1] = mc.GetPosgamma2();

	double vtx[3] = {convp[0].X(), convp[0].Y(), convp[0].Z()};

	art::ServiceHandle<geo::Geometry> geom;
	geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);
	size_t cryoid = geom->FindCryostatAtPosition(vtx);

	art::Handle< std::vector<recob::Hit> > hitListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle;

	//map: conversion point, vec of id clusters in each view
	std::map < size_t, std::vector< size_t > > used;

	double maxdist = 1.0; // 1 cm
	if (geom->HasTPC(idtpc))
	{
		const geo::CryostatGeo& cryostat = geom->Cryostat(cryoid);
		if (evt.getByLabel(fHitsModuleLabel, hitListHandle) &&
	    		evt.getByLabel(fCluModuleLabel, cluListHandle))
		{
			size_t conv = 0;
			while (conv < 2)
			{
				for (size_t view = 0; view < cryostat.MaxPlanes(); view++)
				{
					art::FindManyP< recob::Hit > fbc(cluListHandle, evt, fCluModuleLabel);

					double mindist = maxdist; int clid = 0;
					for (size_t c = 0; c < cluListHandle->size(); ++c)
					{
					
						bool exist = false;
						for (auto const & ids : used)
							for (auto i : ids.second)
								if (i == c) exist = true;
						if (exist) continue;

						std::vector< art::Ptr<recob::Hit> > hits = fbc.at(c);
						if (hits.size() < 20) continue;
						if (hits[0]->WireID().Plane != view) continue;
				
						double dist = getMinDist(hits, convp[conv], view, idtpc.TPC, cryoid);
						if (dist < mindist)
						{
							mindist = dist;
							clid = c;
						}					
					}
					if (mindist < maxdist) used[conv].push_back(clid);			
				}
				conv++;
			}
		}
	}
	bool result = false;
	
	if (used.size() > 1)	 
		for (auto const & ids : used)
		{
			if (ids.second.size() > 1) result = true;
			else {result = false; break;}
		}
	
	return result;
}

double ems::MultiEMShowers::getMinDist(std::vector< art::Ptr<recob::Hit> > const & v, 
								TVector3 const & convmc, 
								size_t view, size_t tpc, size_t cryo)
{
	double mindist = 9999;
	// MC vertex projected to view
	TVector2 proj = pma::GetProjectionToPlane(convmc, view, tpc, cryo);
	
	// loop over hits to find the closest to MC 2d vtx
	for (size_t h = 0; h < v.size(); ++h)
	{
		if ((v[h]->WireID().Plane == view) &&
			(v[h]->WireID().TPC == tpc))
		{
			TVector2 hpoint = pma::WireDriftToCm(v[h]->WireID().Wire, v[h]->PeakTime(), view, tpc, cryo);
			
			double dist = pma::Dist2(proj, hpoint);
			if (dist < mindist) 
			{
				mindist = dist;
			}
		}
	}
	
	return mindist;
}

DEFINE_ART_MODULE(ems::MultiEMShowers)

