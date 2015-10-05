////////////////////////////////////////////////////////////////////////
// Class:       MergeEMShower3D
// Module Type: producer
// File:        MergeEMShower3D_module.cc
// Authors: dorota.stefan@cern.ch robert.sulej@cern.ch	
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
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

namespace ems
{
	class MergeEMShower3D;
	class ShowerInfo;
	class bDistLess;
	class ShowersCollection;
}

class ems::ShowerInfo
{
	public:
	ShowerInfo(int key, int gid, bool hasvtx, double adcsum, recob::Track const& trk);
	double Pointsto(ems::ShowerInfo const& s1) const;
	double Angleto(TVector3 const& pmapoint) const;

	bool HasConPoint(void) const {return fHasVtx;}
	int GetKey(void) const {return fKey;}
	int GetGid(void) const {return fGId;}
	double GetAdcSum(void) const {return fAdcSum;}
	TVector3 GetFront(void) const {return fFront;}
	TVector3 GetEnd(void) const {return fEnd;}
	TVector3 GetDir(void) const {return fDir;}
	std::vector<double> GetDedx(void) const { return fVdqdx; }

	double GetP0Dist(void) const { return fP0Dist; }
	void SetP0Dist(const TVector3 p)
  	{
		if (fHasVtx) fP0Dist = std::sqrt( pma::Dist2(p, fFront) );
		else fP0Dist = std::sqrt( pma::Dist2(p, 0.5 * (fFront + fEnd)) );
	}
 	
	private:
	int fKey;
	int fGId;
	bool fHasVtx;
	double fAdcSum;
	double fP0Dist;

	TVector3 fFront;
	TVector3 fEnd; 
	TVector3 fDir;

	std::vector<double> fVdqdx;
};

ems::ShowerInfo::ShowerInfo(int key, int gid,  bool hasvtx, double adcsum, recob::Track const& trk) :
fKey(key),
fGId(gid),
fAdcSum(adcsum),
fP0Dist(0)
{
	fDir = trk.VertexDirection();
	fFront = trk.Vertex();
	fEnd = trk.End();			
	fHasVtx = hasvtx;

	fVdqdx.push_back(trk.DQdxAtPoint(0, geo::kZ)); 
	fVdqdx.push_back(trk.DQdxAtPoint(0, geo::kV)); 
	fVdqdx.push_back(trk.DQdxAtPoint(0, geo::kU)); 
}

double ems::ShowerInfo::Pointsto(ems::ShowerInfo const& s1) const
{
	double cosine0 = (this->fFront - s1.fFront) * (1.0 / (this->fFront - s1.fFront).Mag()) * s1.fDir;
	double th0 = 180.0F * (std::acos(cosine0)) / TMath::Pi();

	std::vector< std::pair<TVector3, TVector3> > lines;
	lines.push_back(std::pair<TVector3, TVector3>(this->fFront, this->fFront + this->fDir));
	lines.push_back(std::pair<TVector3, TVector3>(s1.fFront, s1.fFront + s1.fDir));

	TVector3 isect, pm;
	pma::SolveLeastSquares3D(lines, isect);
	pm = pma::GetProjectionToSegment(isect, this->fFront, this->fFront +  this->fDir); 
	double cosine1 = (pm - s1.fFront) * (1.0/(pm - s1.fFront).Mag()) * s1.fDir;
	double th1 = 180.0F * (std::acos(cosine1)) / TMath::Pi();

	double thmin = th0;
	if (th1 < th0) thmin = th1;	

	return thmin;
}

double ems::ShowerInfo::Angleto(TVector3 const& pmapoint) const
{
	double cosine = (pmapoint - this->fFront) * (1.0 / (pmapoint - this->fFront).Mag()) * this->fDir;
	double th = 180.0F * (std::acos(cosine)) / TMath::Pi();
	
	return th;	
}

class ems::bDistLess :
	public std::binary_function<const ems::ShowerInfo&, const ems::ShowerInfo&, bool>
{
	public:
	bool operator() (const ems::ShowerInfo& s1, const ems::ShowerInfo& s2)
	{
		return (s1.GetP0Dist() < s2.GetP0Dist());
	}
};

class ems::ShowersCollection
{
	public:
		ShowersCollection(ShowerInfo const& part);

		double MinAngle(ShowersCollection const& col2);
		void Merge(ShowersCollection const& col2);
		void Merge(ShowerInfo const& src);

		double Angle(TVector3 p0, TVector3 test);

		double Angle(TVector3 test);

		size_t Size() {return fParts.size();}

		std::vector< ShowerInfo > const & GetParts(void) const { return fParts; }

		bool IsClean() {return Clean;}

		TVector3 Front(); 
		TVector3 Dir(); 
		std::vector<double> DeDx(void);

		ShowerInfo first;
	
	private:
		bool Clean;	
		int GId;
		std::vector< ShowerInfo > fParts;
};

ems::ShowersCollection::ShowersCollection(ShowerInfo const& part) :
first(part)
{
	Clean = true;
	GId = part.GetGid();
	if (GId == 9999)
		Clean = false;
	
	fParts.push_back(part);
}


TVector3 ems::ShowersCollection::Front()
{
	if (fParts.size()) return fParts.front().GetFront();
	else return TVector3(0, 0, 0);
}

TVector3 ems::ShowersCollection::Dir()
{
	if (fParts.size()) return fParts.front().GetDir();
	else return TVector3(0, 0, 0);
}

std::vector<double> ems::ShowersCollection::DeDx()
{
	if (fParts.size()) return fParts.front().GetDedx();
	else return std::vector<double>(0.0);
}

double ems::ShowersCollection::Angle(TVector3 p0, TVector3 test)
{
	p0 = Front() - p0;
	p0 *= 1.0 / p0.Mag();
	test -= Front();
	test *= 1.0 / test.Mag();
	double c = p0 * test;
	if (c > 1.0) c = 1.0;
	return 180.0 * acos(c) / TMath::Pi();
}

double ems::ShowersCollection::Angle(TVector3 test)
{
	test -= Front();
	test *= 1.0 / test.Mag();
	double c = Dir() * test;
	if (c > 1.0) c = 1.0;
	return 180.0 * acos(c) / TMath::Pi();
}


double ems::ShowersCollection::MinAngle(ShowersCollection const& col2)
{
	double a, a_min = 360.0;
	for (size_t i = 0; i < fParts.size(); ++i)
	{
		const ShowerInfo& temp = fParts[i];
		a = col2.first.Pointsto(temp);
		if (a < a_min) a_min = a;
	}
	return a_min;
}

void ems::ShowersCollection::Merge(ShowersCollection const& col2)
{
	for (size_t i = 0; i < col2.fParts.size(); ++i)
		for (size_t j = 0; j < fParts.size(); ++j)
			if (fParts[j].GetGid() != col2.fParts[i].GetGid())	
			{
				Clean = false;
				break;
			}

	for (size_t i = 0; i < col2.fParts.size(); ++i)
		fParts.push_back(col2.fParts[i]);
	
}

void ems::ShowersCollection::Merge(ShowerInfo const& src)
{
	if (fParts.size())
	{
		for (auto const & p : fParts)
			if (p.GetGid() != src.GetGid())
			{
				Clean = false;
				break;
			}
	}
	else
	{
		Clean = true;
		GId = src.GetGid();
	}
	fParts.push_back(src);
}

class ems::MergeEMShower3D : public art::EDProducer {
public:
  explicit MergeEMShower3D(fhicl::ParameterSet const & p);

  MergeEMShower3D(MergeEMShower3D const &) = delete;
  MergeEMShower3D(MergeEMShower3D &&) = delete;
  MergeEMShower3D & operator = (MergeEMShower3D const &) = delete;
  MergeEMShower3D & operator = (MergeEMShower3D &&) = delete;

	void beginJob();

	void produce(art::Event & e) override;

	void reconfigure(fhicl::ParameterSet const& p);

	void mcinfo(art::Event & evt);

	int getClusterBestId(const std::vector< art::Ptr<recob::Hit> >& v);

	int getGammaId(art::Event & evt, const size_t t);

	double getCos3D(const TVector3& p0, const recob::Track& trk);

	double getClusterAdcSum(const std::vector< art::Ptr<recob::Hit> >& v);

	std::vector< ems::ShowersCollection > collectshowers(art::Event & evt, std::vector< ems::ShowerInfo > & showers, bool refpoint);

	TVector3 getBestPoint(
	const std::vector<recob::Track>& tracks,
	const std::vector< ShowerInfo >& showers,
	const TVector3& p0, const TVector3& p1,
	double step);

private:

	TTree* fEvTree;
	TTree* fClTree;

	int fVtxIndex;

	int fNParts[2];
	int fNPMA; int fNConv; int fNTot;
	int fHasConvPt;
	int fWhat, fEvWhat;
	int fNMerged, fNCleanMerged;

	double fDistvtxmcreco;
	double fMcMom, fTrkLen, fAdcSum;
	int fEvNumber;

	TVector3 fPi0vtx;
	TVector3 fRefPoint;

	std::string fHitsModuleLabel;
	std::string fCluModuleLabel;
	std::string fTrk3DModuleLabel;
	std::string fVtxModuleLabel;

	double fNarrowConeAngle;
	double fWideConeAngle;

};

ems::MergeEMShower3D::MergeEMShower3D(fhicl::ParameterSet const & p) 
{
	reconfigure(p);

	produces< std::vector<recob::Shower> >();
	produces< std::vector<recob::Vertex> >();

	produces< art::Assns<recob::Shower, recob::Vertex> >(); 
	produces< art::Assns<recob::Shower, recob::Cluster> >();
	produces< art::Assns<recob::Shower, recob::Hit> >();
}

void ems::MergeEMShower3D::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;

	fEvTree = tfs->make<TTree>("ShowerTestEv", "pi0 tests");
	fEvTree->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fEvTree->Branch("fNPartsA", &fNParts[0], "fNPartsA/I");
	fEvTree->Branch("fNPartsB", &fNParts[1], "fNPartsB/I");
	fEvTree->Branch("fNTot", &fNTot, "fNTot/I");
	fEvTree->Branch("fNPMA", &fNPMA, "fNPMA/I");
	fEvTree->Branch("fNConv", &fNConv, "fNConv/I");
	fEvTree->Branch("fEvWhat", &fEvWhat, "fEvWhat/I");
	fEvTree->Branch("fMcMom", &fMcMom, "fMcMom/D");

	fEvTree->Branch("fDistvtxmcreco", &fDistvtxmcreco, "fDistvtxmcreco/D");
	fEvTree->Branch("fNMerged", &fNMerged, "fNMerged/I");
	fEvTree->Branch("fNCleanMerged", &fNCleanMerged, "fNCleanMerged/I");
}

void ems::MergeEMShower3D::reconfigure(fhicl::ParameterSet const & p)
{
	fHitsModuleLabel = p.get< std::string >("HitsModuleLabel");
	fCluModuleLabel = p.get< std::string >("ClustersModuleLabel");
	fTrk3DModuleLabel = p.get< std::string >("Trk3DModuleLabel");
	fVtxModuleLabel = p.get< std::string >("VtxModuleLabel");

	fNarrowConeAngle = p.get< double >("NarrowConeAngle");
	fWideConeAngle = p.get< double >("WideConeAngle");

	return;
}

void ems::MergeEMShower3D::produce(art::Event & evt)
{
	std::unique_ptr< std::vector< recob::Shower > > cascades(new std::vector< recob::Shower >);
	std::unique_ptr< std::vector< recob::Vertex > > vertices(new std::vector< recob::Vertex >);

	std::unique_ptr< art::Assns< recob::Shower, recob::Vertex > > shs2vtx(new art::Assns< recob::Shower, recob::Vertex >);
	std::unique_ptr< art::Assns< recob::Shower, recob::Cluster > > shs2cl(new art::Assns< recob::Shower, recob::Cluster >);
	std::unique_ptr< art::Assns< recob::Shower, recob::Hit > > shs2hit(new art::Assns< recob::Shower, recob::Hit >);

	std::vector< ShowersCollection > gammawithconv;

	art::Handle< std::vector<recob::Track> > trkListHandle;
	art::Handle< std::vector<recob::Vertex> > vtxListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	art::Handle< std::vector<recob::Hit> > hitListHandle;
	if (evt.getByLabel(fTrk3DModuleLabel, trkListHandle) &&
	    evt.getByLabel(fVtxModuleLabel, vtxListHandle) &&
	    evt.getByLabel(fCluModuleLabel, cluListHandle) &&
	    evt.getByLabel(fHitsModuleLabel, hitListHandle))
	{

		fEvNumber = evt.id().event();
		fNMerged = 0; fNCleanMerged = 0; fNParts[0] = 0; fNParts[1] = 0;
		fNPMA = 0; fNConv = 0; fNTot = 0; fMcMom = 0; fWhat = 0;
		fDistvtxmcreco = -10.0;

		art::FindManyP< recob::Cluster > cluFromTrk(trkListHandle, evt, fTrk3DModuleLabel);
		art::FindManyP< recob::Vertex > vtxFromTrk(trkListHandle, evt, fVtxModuleLabel);
		art::FindManyP< recob::Hit > hitFromClu(cluListHandle, evt, fCluModuleLabel);

		// all reconstructed shower fragments
		std::vector< ShowerInfo > showers;
		for (size_t t = 0; t < trkListHandle->size(); ++t)
		{
			const recob::Track& trk = (*trkListHandle)[t];

			auto src_clu_list = cluFromTrk.at(t);
			double adcsum = 0.0;
			for (size_t c = 0; c < src_clu_list.size(); ++c)  
			{
				std::vector< art::Ptr<recob::Hit> > v = hitFromClu.at(src_clu_list[c].key());
				adcsum += getClusterAdcSum(v);
			}

			auto cnv = vtxFromTrk.at(t);
			if (cnv.size()) 
			{
				int gid = getGammaId(evt, t);
				ShowerInfo si(t, gid, true, adcsum, trk);
				showers.push_back(si);	
			}
			else	
			{	
				int gid = getGammaId(evt, t);
				ShowerInfo si(t, gid, false, adcsum, trk);
				showers.push_back(si);
			}
		}

		fNTot = showers.size();

		//// collect shower fragments 
		gammawithconv = collectshowers(evt, showers, true); // true switches on refpoint
	
		// with pma segments reconstructed
		// procceed with pma segments when two conversions 
		std::vector< ShowerInfo > pmaseg;
		for (size_t i = 0; i < showers.size(); ++i)
			if (!showers[i].HasConPoint())
				pmaseg.push_back(showers[i]);
			
		
		fNPMA = pmaseg.size();

		const double bigcone = fWideConeAngle; // degree.

		if (gammawithconv.size())
			for (size_t i = 0; i < pmaseg.size(); i++)
			{
				double a_min = bigcone; size_t best = 0;
				for (size_t j = 0; j < gammawithconv.size(); j++)
				{
					TVector3 halfpoint = (pmaseg[i].GetFront() + pmaseg[i].GetEnd()) * 0.5;
					double a = gammawithconv[j].first.Angleto(halfpoint);
					if (a < a_min)
					{
						a_min = a; best = j;
					}
				}
				if (a_min < bigcone)
					gammawithconv[best].Merge(pmaseg[i]);		
			}
		
		mcinfo(evt);
		

		for (size_t i = 0; i < gammawithconv.size(); ++i)
				if (gammawithconv[i].IsClean()) fNCleanMerged++;
		
		fNMerged = gammawithconv.size(); 
		if (gammawithconv.size() == 2)
		{
			fNParts[0] = gammawithconv[0].Size();
			fNParts[1] = gammawithconv[1].Size();
		}

		fDistvtxmcreco = std::sqrt(pma::Dist2(fRefPoint, fPi0vtx));	
		fEvTree->Fill();

		fVtxIndex = 0;
		for (size_t i = 0; i < gammawithconv.size(); ++i)
		{
			int id = i;
			TVector3 v0(0., 0., 0.);
			TVector3 dir = gammawithconv[i].Dir();
			TVector3 front = gammawithconv[i].Front();
			 
			std::vector<double> dedx = gammawithconv[i].DeDx();
			std::vector< double > vd;
			recob::Shower cas(
				dir, v0, front, v0,
				vd, vd, dedx, vd, 0, id);

			cascades->push_back(cas);

			std::vector< art::Ptr<recob::Cluster> > cls;
			std::vector< art::Ptr<recob::Hit> > hits;
			for (size_t p = 0; p < gammawithconv[i].Size(); p++)
			{
				int trkKey = gammawithconv[i].GetParts()[p].GetKey();
				auto src_clu_list = cluFromTrk.at(trkKey);

				for (size_t c = 0; c < src_clu_list.size(); c++)
				{
					cls.push_back(src_clu_list[c]);

					auto v = hitFromClu.at(src_clu_list[c].key());
					for (size_t h = 0; h < v.size(); h++) hits.push_back(v[h]);
				}

				auto ver_list = vtxFromTrk.at(trkKey);
			}

			double vtx_pos[3] = {front.X(), front.Y(), front.Z()};
			vertices->push_back(recob::Vertex(vtx_pos, fVtxIndex));
			fVtxIndex++;

			if (vertices->size())
			{
				size_t vtx_idx = (size_t)(vertices->size() - 1);
				util::CreateAssn(*this, evt, *cascades, *vertices, *shs2vtx, vtx_idx, vtx_idx + 1);
			}

			util::CreateAssn(*this, evt, *cascades, cls, *shs2cl);
			util::CreateAssn(*this, evt, *cascades, hits, *shs2hit);

		}

	}

	evt.put(std::move(cascades));
	evt.put(std::move(vertices));

	evt.put(std::move(shs2vtx));
	evt.put(std::move(shs2cl));
	evt.put(std::move(shs2hit));
}

std::vector< ems::ShowersCollection > ems::MergeEMShower3D::collectshowers(art::Event & evt, std::vector< ems::ShowerInfo > & showers, bool refpoint)
{
	std::vector< ems::ShowersCollection > gammawithconv;

	const double smallcone = fNarrowConeAngle; // degree.
	bool merge = true;
	if (refpoint)
	{
		art::ServiceHandle< geo::Geometry > geom;
		art::Handle< std::vector<recob::Track> > trkListHandle;		
		if (evt.getByLabel(fTrk3DModuleLabel, trkListHandle)) 
		{
			double dsize[6];
			geom->CryostatBoundaries(dsize, 0);
			TVector3 geoP0(dsize[0], dsize[2], dsize[4]), geoP1(dsize[1], dsize[3], dsize[5]);
			TVector3 pov = getBestPoint(*trkListHandle, showers, geoP0, geoP1, 5.0);
			fRefPoint = pov;

			for (size_t is = 0; is < showers.size(); is++)
				showers[is].SetP0Dist(pov);
			std::sort(showers.begin(), showers.end(), ems::bDistLess());

			for (size_t is = 0; is < showers.size(); is++)
			{
				size_t m_idx = 0;
				double a, a_min = 360.0;
				for (size_t m = 0; m < gammawithconv.size(); ++m)
				{
					a = gammawithconv[m].Angle(showers[is].GetFront());
					if ((a < fNarrowConeAngle) && (a < a_min))
					{
						m_idx = m; a_min = a; 
					}

					std::vector< std::pair<TVector3, TVector3> > lines;
					lines.push_back(std::pair<TVector3, TVector3>(
							showers[is].GetFront(), showers[is].GetFront() + showers[is].GetDir()));
					lines.push_back(std::pair<TVector3, TVector3>(
							gammawithconv[m].Front(), gammawithconv[m].Front() +  gammawithconv[m].Dir()));

					TVector3 isect, pm;
					pma::SolveLeastSquares3D(lines, isect);
					pm = pma::GetProjectionToSegment(isect,
							gammawithconv[m].Front(), gammawithconv[m].Front() +  gammawithconv[m].Dir());

					if (pma::Dist2(pov, gammawithconv[m].Front()) < pma::Dist2(pov, pm))
					{
						a = gammawithconv[m].Angle(isect);
						if ((a < fNarrowConeAngle) && (a < a_min))
						{
							m_idx = m; a_min = a;
						}
					}
				}
				if (a_min < fNarrowConeAngle) 
				{
					gammawithconv[m_idx].Merge(showers[is]);
				}
				else if (showers[is].HasConPoint())
				{
					ems::ShowersCollection sc(showers[is]);
					gammawithconv.push_back(sc);
				}

			}
		}
	}
	else
	{
		for (size_t i = 0; i < showers.size(); ++i)
			if (showers[i].HasConPoint())
			{
				ems::ShowersCollection sc(showers[i]);
				gammawithconv.push_back(sc);
			}

		fNConv = gammawithconv.size();
		while (merge)
		{
			merge = false;
			size_t i = 0;
			while (i < gammawithconv.size())
			{
				size_t best = 0; double a_min = smallcone; 
				for (size_t j = 0; j < gammawithconv.size(); j++)
					if (i != j)
					{
						double a = gammawithconv[i].MinAngle(gammawithconv[j]);
						if (a < a_min)
						{
							a_min = a; best = j; merge = true;
						}
					}
				if (merge)
				{
					gammawithconv[i].Merge(gammawithconv[best]);
					gammawithconv.erase(gammawithconv.begin() + best);
					break;
				}
				i++;
			}
		}
	}
	return  gammawithconv;
}


int ems::MergeEMShower3D::getClusterBestId(const std::vector< art::Ptr<recob::Hit> >& v)
{
	art::ServiceHandle< cheat::BackTracker > bt;

	std::map< int, size_t > ids;
	for (auto ptr : v)
	{
		auto hid = bt->HitToTrackID(ptr);
		if (hid.size()) ids[hid.front().trackID]++;
	}

	int best_id = 9999;
	size_t max = 0;
	for (auto p : ids)
	{
		if ((p.second > (int)(0.8 * v.size())) && (p.second > max))
		{
			max = p.second; best_id = p.first;
		}
	}

	return best_id;
}

void ems::MergeEMShower3D::mcinfo(art::Event & evt)
{
	int pi0_idx = -1;

	art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector< art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel("generator", mctruthListHandle))
	{
		art::fill_ptr_vector(mclist, mctruthListHandle);
		if (mclist.size())
		{
			size_t np = mclist[0]->NParticles();
			for (size_t i = 0; i < np; ++i)
				if (mclist[0]->GetParticle(i).PdgCode() == 111)
					{pi0_idx = i; break;}
		}
	}

	if (pi0_idx < 0) return; 

	const simb::MCParticle & pi0 = mclist[0]->GetParticle(pi0_idx);
	TVector3 pi0_vtx(pi0.Vx(), pi0.Vy(), pi0.Vz());
	fPi0vtx = pi0_vtx;
	fMcMom = pi0.P();
}

int ems::MergeEMShower3D::getGammaId(art::Event & evt, const size_t t)
{
	fWhat = 0;
	int cid = 0, gid = 0;
	art::Handle< std::vector<recob::Track> > trkListHandle;
	art::Handle< std::vector<recob::Vertex> > vtxListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	art::Handle< std::vector<recob::Hit> > hitListHandle;
	if (evt.getByLabel(fTrk3DModuleLabel, trkListHandle) &&
	    evt.getByLabel(fVtxModuleLabel, vtxListHandle) &&
	    evt.getByLabel(fCluModuleLabel, cluListHandle) &&
	    evt.getByLabel(fHitsModuleLabel, hitListHandle))
	{
		art::FindManyP< recob::Cluster > cluFromTrk(trkListHandle, evt, fTrk3DModuleLabel);
		art::FindManyP< recob::Vertex > vtxFromTrk(trkListHandle, evt, fVtxModuleLabel);
		art::FindManyP< recob::Hit > hitFromClu(cluListHandle, evt, fCluModuleLabel);

		auto src_clu_list = cluFromTrk.at(t);
		for (size_t c = 0; c < src_clu_list.size(); ++c)  
		{
			std::vector< art::Ptr<recob::Hit> > v = hitFromClu.at(src_clu_list[c].key());
			cid = getClusterBestId(v);

			if ((fWhat == 0) && (cid == 9999)) fWhat = 1; // confused 2D
			if ((fWhat == 0) && (gid == 0)) gid = cid;    
			if ((fWhat == 0) && (cid != gid)) fWhat = 2;  // confused 3D 
		}	
	}

	if (fWhat == 0) return gid;
	else return 9999;
}

TVector3 ems::MergeEMShower3D::getBestPoint(
	const std::vector< recob::Track >& tracks,
	const std::vector< ShowerInfo >& showers,
	const TVector3& p0, const TVector3& p1, double step)
{
	TVector3 best, p;

	double f, fmin = 1.0e10;

	double dx = step; 
	double dy = step; 
	double dz = step; 

	double x0 = p0.X();
	while (x0 < p1.X())
	{
		double y0 = p0.Y();
		while (y0 < p1.Y())
		{
			double z0 = p0.Z();
			while (z0 < p1.Z())
			{
				f = 0.0;
				//size_t nOK = 0;
				TVector3 mid(0., 0., 0.);

				p.SetXYZ(x0, y0, z0);
				for (size_t t = 0; t < tracks.size(); t++)
				//	if (showers[t].OK)
				{
					auto const & trk = tracks[t];

					double cos = -acos( getCos3D(p, trk) ) / (0.5*TMath::Pi());
					//double cos = getCos3D(p, trk);

					if (showers[t].HasConPoint()) cos *= 3.0;
					cos *= sqrt( showers[t].GetAdcSum() );

					mid += 0.5 * (trk.Vertex() + trk.End());

					f += cos;
				//	nOK++;
				}
				//if (!nOK) return best;

			//	f /= nOK;
			//	mid *= 1.0 / nOK;

				f = -f + 0.0001 * sqrt(pma::Dist2(p, mid));
				if (f < fmin)
				{
					fmin = f; best = p;
				}

				z0 += dz;
			}
			y0 += dy;
		}
		x0 += dx;
	}
	return best;
}

double ems::MergeEMShower3D::getCos3D(const TVector3& p0, const recob::Track& trk)
{
	TVector3 p1, dir;

	if (pma::Dist2(p0, trk.Vertex()) < pma::Dist2(p0, trk.End()))
		p1 = trk.Vertex();
	else
		p1 = trk.End();

	p1 -= p0;
	double m = p1.Mag();
	if (m > 0.0)
	{
		p1 *= 1.0 / m;
		double c = fabs(p1 * trk.VertexDirection());
		if (c > 1.0) c = 1.0;
		return c;
	}
	else return 0.0;
}

double ems::MergeEMShower3D::getClusterAdcSum(std::vector< art::Ptr<recob::Hit> > const & v)
{
	double sum = 0.0;
	for (auto ptr : v)
	{
		sum += ptr->SummedADC();
	}
	return sum;
}

DEFINE_ART_MODULE(ems::MergeEMShower3D)
