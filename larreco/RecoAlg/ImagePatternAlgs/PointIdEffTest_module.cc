////////////////////////////////////////////////////////////////////////
// Class:       PointIdEffTest
// Module Type: analyzer
// File:        PointIdEffTest_module.cc
//
// Author: dorota.stefan@cern.ch
//
// Generated at Fri Apr 29 06:42:27 2016 by Dorota Stefan using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/Utilities/DatabaseUtil.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/ArtDataHelper/MVAReader.h"

#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>

#include <cmath>

#define MVA_LENGTH 4

namespace nnet {

class PointIdEffTest : public art::EDAnalyzer {
public:
    enum EId { kShower = 0, kTrack = 1, kMichel = 2 };

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<calo::CalorimetryAlg::Config> CalorimetryAlg {
			Name("CalorimetryAlg"), Comment("Used to calculate electron lifetime correction.")
		};

		fhicl::Atom<art::InputTag> SimModuleLabel { Name("SimModuleLabel"), Comment("Simulation producer") };

		fhicl::Atom<art::InputTag> PfpModuleLabel { Name("PfpModuleLabel"), Comment("PFP producer tag, to compare with NNet results") };

		fhicl::Atom<art::InputTag> NNetModuleLabel { Name("NNetModuleLabel"), Comment("NNet outputs tag") };

		fhicl::Atom<bool> SaveHitsFile { Name("SaveHitsFile"), Comment("Dump hits info to text file") };

		fhicl::Atom<unsigned int> View { Name("View"), Comment("Which view is evaluated") };
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit PointIdEffTest(Parameters const& config);

    PointIdEffTest(PointIdEffTest const &) = delete;
    PointIdEffTest(PointIdEffTest &&) = delete;
    PointIdEffTest & operator = (PointIdEffTest const &) = delete;
    PointIdEffTest & operator = (PointIdEffTest &&) = delete;

	virtual void beginRun(const art::Run& run) override;

	virtual void beginJob() override;
	virtual void endJob() override;

    virtual void analyze(art::Event const & e) override;

private:
    void cleanup(void);

    void countTruthDep(
        const std::vector< sim::SimChannel > & channels,
        float & emLike, float & trackLike) const;

    void countPfpDep(
        const std::vector< recob::PFParticle > & pfparticles,
        const art::FindManyP<recob::Cluster> & pfpclus,
        const art::FindManyP<recob::Hit> & cluhits,
        float & emLike, float & trackLike) const;

    bool isMuonDecaying(
        const simb::MCParticle & particle,
        const std::unordered_map< int, const simb::MCParticle* > & particleMap) const;

	int testCNN(
	    const std::vector< sim::SimChannel > & channels,
        const std::vector< art::Ptr<recob::Hit> > & hits,
        const std::array<float, MVA_LENGTH> & cnn_out,
        const std::vector< anab::FeatureVector<MVA_LENGTH> > & hit_outs,
        size_t cidx);

	int fRun, fEvent;
    float fMcDepEM, fMcDepTrack, fMcFractionEM;
    float fPfpDepEM, fPfpDepTrack;
    float fHitEM_0p5, fHitTrack_0p5, fHitMichel_0p5, fHitMcFractionEM;
    float fHitEM_mc, fpEM;
    float fHitMichel_mc, fpMichel_hit, fpMichel_cl;
    float fOutTrk, fOutEM, fOutNone;
    float fHitEM_0p85, fHitTrack_0p85;
    float fTotHit, fCleanHit;

    //const int kEffSize = 100;
    float fHitsEM_OK_0p5[100], fHitsTrack_OK_0p5[100];
    float fHitsMichel_OK_0p5[100], fHitsMichel_False_0p5[100];
    float fHitsEM_OK_0p85[100], fHitsTrack_OK_0p85[100];
    float fHitRecoEM[100], fHitRecoFractionEM[100];

	int fMcPid;
	int fClSize;
	int fPure;
	float fPidValue; // P(track-like)
	

	int fTrkOk[100], fTrkBad[100];
	int fShOk[100], fShBad[100];
	int fNone, fTotal;

	double fElectronsToGeV;
	
	int fTrkLikeIdx, fEmLikeIdx, fNoneIdx, fMichelLikeIdx;

	TTree *fEventTree, *fClusterTree, *fHitTree;

	std::ofstream fHitsOutFile;

	unsigned int fView;

	geo::GeometryCore const* fGeometry;

	std::unordered_map< int, const simb::MCParticle* > fParticleMap;

    calo::CalorimetryAlg fCalorimetryAlg;
	art::InputTag fSimulationProducerLabel;
	art::InputTag fPfpModuleLabel;
	art::InputTag fNNetModuleLabel;
	bool fSaveHitsFile;
};

} // namespace nnet

nnet::PointIdEffTest::PointIdEffTest(nnet::PointIdEffTest::Parameters const& config) : art::EDAnalyzer(config),
	fMcPid(-1), fClSize(0),

    fTrkLikeIdx(-1), fEmLikeIdx(-1), fNoneIdx(-1), fMichelLikeIdx(-1),

	fView(config().View()),
	fCalorimetryAlg(config().CalorimetryAlg()),
	fSimulationProducerLabel(config().SimModuleLabel()),
	fPfpModuleLabel(config().PfpModuleLabel()),
	fNNetModuleLabel(config().NNetModuleLabel()),
	fSaveHitsFile(config().SaveHitsFile())
{
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());
}

void nnet::PointIdEffTest::beginRun(const art::Run&)
{
	art::ServiceHandle< sim::LArG4Parameters > larParameters;
	fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

void nnet::PointIdEffTest::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;

    fEventTree = tfs->make<TTree>("event","event info");
    fEventTree->Branch("fRun", &fRun, "fRun/I");
    fEventTree->Branch("fEvent", &fEvent, "fEvent/I");
    fEventTree->Branch("fMcDepEM", &fMcDepEM, "fMcDepEM/F");
    fEventTree->Branch("fMcDepTrack", &fMcDepTrack, "fMcDepTrack/F");
    fEventTree->Branch("fMcFractionEM", &fMcFractionEM, "fMcFractionEM/F");
    fEventTree->Branch("fPfpDepEM", &fPfpDepEM, "fPfpDepEM/F");
    fEventTree->Branch("fPfpDepTrack", &fPfpDepTrack, "fPfpDepTrack/F");
    fEventTree->Branch("fHitEM_0p5", &fHitEM_0p5, "fHitEM_0p5/F");
    fEventTree->Branch("fHitMichel_0p5", &fHitMichel_0p5, "fHitMichel_0p5/F");
    fEventTree->Branch("fHitTrack_0p5", &fHitTrack_0p5, "fHitTrack_0p5/F");
    fEventTree->Branch("fHitMcFractionEM", &fHitMcFractionEM, "fHitMcFractionEM/F");
    fEventTree->Branch("fHitEM_0p85", &fHitEM_0p85, "fHitEM_0p85/F");
    fEventTree->Branch("fHitTrack_0p85", &fHitTrack_0p85, "fHitTrack_0p85/F");
    fEventTree->Branch("fCleanHit", &fCleanHit, "fCleanHit/F");

    fEventTree->Branch("fHitsEM_OK_0p5", fHitsEM_OK_0p5, "fHitsEM_OK_0p5[100]/F");
    fEventTree->Branch("fHitsTrack_OK_0p5", fHitsTrack_OK_0p5, "fHitsTrack_OK_0p5[100]/F");
    fEventTree->Branch("fHitsMichel_OK_0p5", fHitsMichel_OK_0p5, "fHitsMichel_OK_0p5[100]/F");
    fEventTree->Branch("fHitsMichel_False_0p5", fHitsMichel_False_0p5, "fHitsMichel_False_0p5[100]/F");

    fEventTree->Branch("fHitsEM_OK_0p85", fHitsEM_OK_0p85, "fHitsEM_OK_0p85[100]/F");
    fEventTree->Branch("fHitsTrack_OK_0p85", fHitsTrack_OK_0p85, "fHitsTrack_OK_0p85[100]/F");

    fEventTree->Branch("fHitRecoEM", fHitRecoEM, "fHitRecoEM[100]/F");
    fEventTree->Branch("fHitRecoFractionEM", fHitRecoFractionEM, "fHitRecoFractionEM[100]/F");


	fClusterTree = tfs->make<TTree>("cluster","clusters info");
	fClusterTree->Branch("fMcPid", &fMcPid, "fMcPid/I");
	fClusterTree->Branch("fClSize", &fClSize, "fClSize/I");
	fClusterTree->Branch("fPidValue", &fPidValue, "fPidValue/F");
	fClusterTree->Branch("fpMichel_cl", &fpMichel_cl, "fpMichel_cl/F");

    fHitTree = tfs->make<TTree>("hits","hits info");
    fHitTree->Branch("fRun", &fRun, "fRun/I");
    fHitTree->Branch("fEvent", &fEvent, "fEvent/I");
    fHitTree->Branch("fHitEM_mc", &fHitEM_mc, "fHitEM_mc/F");
    fHitTree->Branch("fpEM", &fpEM, "fpEM/F");
    fHitTree->Branch("fPidValue", &fPidValue, "fPidValue/F");
    fHitTree->Branch("fHitMichel_mc", &fHitMichel_mc, "fHitMichel_mc/F");
    fHitTree->Branch("fpMichel_hit", &fpMichel_hit, "fpMichel_hit/F");
    fHitTree->Branch("fOutTrk", &fOutTrk, "fOutTrk/F");
    fHitTree->Branch("fOutEM", &fOutEM, "fOutEM/F");
    fHitTree->Branch("fOutNone", &fOutNone, "fOutNone/F");

	if (fSaveHitsFile) fHitsOutFile.open("hits_pid.prn");

    fNone = 0; fTotal = 0;
    for (size_t i = 0; i < 100; ++i)
    {
        fShOk[i] = 0; fShBad[i] = 0;
        fTrkOk[i] = 0; fTrkBad[i] = 0;
    }
}

void nnet::PointIdEffTest::endJob()
{
	if (fSaveHitsFile) fHitsOutFile.close();

    art::ServiceHandle<art::TFileService> tfs;

    TTree *thrTree = tfs->make<TTree>("threshold","error rate vs threshold");

    float thr, shErr, trkErr;
    thrTree->Branch("thr", &thr, "thr/F");
    thrTree->Branch("shErr", &shErr, "shErr/F");
    thrTree->Branch("trkErr", &trkErr, "trkErr/F");

    for (size_t i = 0; i < 100; ++i)
    {
        thr = 0.01 * i;
        shErr = fShBad[i] / float(fShBad[i] + fShOk[i]);
        trkErr = fTrkBad[i] / float(fTrkBad[i] + fTrkOk[i]);
        thrTree->Fill();

        //std::cout << "Threshold " << thr << "  fShErr " << shErr << " fTrkErr " << trkErr << " sum:" << shErr + trkErr <<  std::endl;
	}
	// std::cout << "Total " << fTotal << std::endl;
}

void nnet::PointIdEffTest::cleanup(void)
{
    fParticleMap.clear();

    fMcDepEM = 0; fMcDepTrack = 0;
    fMcFractionEM = 0;
    fPfpDepEM = 0; fPfpDepTrack = 0;
    fHitEM_0p5 = 0; fHitTrack_0p5 = 0;
    fHitEM_0p85 = 0; fHitTrack_0p85 = 0;
    fHitMichel_0p5 = 0;
    fHitMcFractionEM = 0;
    fTotHit = 0; fCleanHit = 0;

    for (size_t i = 0; i < 100; ++i)
    {
        fHitsEM_OK_0p5[i] = 0; fHitsTrack_OK_0p5[i] = 0;
        fHitsEM_OK_0p85[i] = 0; fHitsTrack_OK_0p85[i] = 0;
        fHitsMichel_OK_0p5[i] = 0; fHitsMichel_False_0p5[i] = 0;
        fHitRecoEM[i] = 0; fHitRecoFractionEM[i] = 0;
    }

    fTrkLikeIdx = -1; fEmLikeIdx = -1; fNoneIdx = -1; fMichelLikeIdx = -1;
}

void nnet::PointIdEffTest::analyze(art::Event const & e)
{
    cleanup(); // remove everything from members

	fRun = e.run();
	fEvent = e.id().event();
	// std::cout << "event " << fEvent << std::endl;

	// access to MC information
	
	// MC particles list
	auto particleHandle = e.getValidHandle< std::vector<simb::MCParticle> >(fSimulationProducerLabel);
	for (auto const& particle : *particleHandle) { fParticleMap[particle.TrackId()] = &particle; }

	// SimChannels
	auto simChannelHandle = e.getValidHandle< std::vector<sim::SimChannel> >(fSimulationProducerLabel);
	countTruthDep(*simChannelHandle, fMcDepEM, fMcDepTrack);

    // PFParticle selection results
	art::Handle< std::vector<recob::PFParticle> > pfpHandle;
	if (e.getByLabel(fPfpModuleLabel, pfpHandle))
	{
	    auto cluHandle = e.getValidHandle< std::vector<recob::Cluster> >(fPfpModuleLabel);
	    const art::FindManyP<recob::Cluster> clusFromPfps(pfpHandle, e, fPfpModuleLabel);
	    const art::FindManyP<recob::Hit> hitsFromClus(cluHandle, e, fPfpModuleLabel);
	    countPfpDep(*pfpHandle, clusFromPfps, hitsFromClus, fPfpDepEM, fPfpDepTrack);
    }

    // output from cnn's

    anab::MVAReader<recob::Hit, MVA_LENGTH> hitResults(e, fNNetModuleLabel);                     // hit-by-hit outpus just to be dumped to file for debugging
    fTrkLikeIdx = hitResults.getIndex("track");
    fEmLikeIdx = hitResults.getIndex("em");
    fNoneIdx = hitResults.getIndex("none");
    fMichelLikeIdx = hitResults.getIndex("michel");
    if ((fTrkLikeIdx < 0) || (fEmLikeIdx < 0))
    {
        throw cet::exception("PointIdEffTest") << "No em/track labeled columns in MVA data products." << std::endl;
    }

    auto cluResults = anab::MVAReader<recob::Cluster, MVA_LENGTH>::create(e, fNNetModuleLabel);  // outputs for clusters recovered in not-throwing way 
    if (cluResults)
    {
    	const art::FindManyP<recob::Hit> hitsFromClusters(cluResults->dataHandle(), e, cluResults->dataTag());

	    for (size_t c = 0; c < cluResults->size(); ++c)
	    {
	        const recob::Cluster & clu = cluResults->item(c);

            if (clu.Plane().Plane != fView) continue;

    	    const std::vector< art::Ptr<recob::Hit> > & hits = hitsFromClusters.at(c);
    	    std::array<float, MVA_LENGTH> cnn_out = cluResults->getOutput(c);

            testCNN(*simChannelHandle, hits, cnn_out, hitResults.outputs(), c); // test hits in the cluster
    	}

        if (fTotHit > 0) fCleanHit = fCleanHit / fTotHit;
        else fCleanHit = 0;

        double totMcDep = fMcDepEM + fMcDepTrack;
        if (totMcDep) fMcFractionEM = fMcDepEM / totMcDep;
        else fMcFractionEM = 0;

        double totEmTrk0p5 = fHitEM_0p5 + fHitTrack_0p5;
        if (totEmTrk0p5 > 0) fHitMcFractionEM = fHitEM_0p5 / totEmTrk0p5;
        else fHitMcFractionEM = 0;

        for (size_t i = 0; i < 100; ++i)
        {
            if (fHitEM_0p5 > 0) fHitsEM_OK_0p5[i] /= fHitEM_0p5;
            else fHitsEM_OK_0p5[i] = 0;

            if (fHitTrack_0p5 > 0) fHitsTrack_OK_0p5[i] /= fHitTrack_0p5;
            else fHitsTrack_OK_0p5[i] = 0;

            if (fHitEM_0p85 > 0) fHitsEM_OK_0p85[i] /= fHitEM_0p85;
            else fHitsEM_OK_0p85[i] = 0;

            if (fHitTrack_0p85 > 0) fHitsTrack_OK_0p85[i] /= fHitTrack_0p85;
            else fHitsTrack_OK_0p85[i] = 0;


            if (totEmTrk0p5 > 0) fHitRecoFractionEM[i] = fHitRecoEM[i] / totEmTrk0p5;
            else fHitRecoFractionEM[i] = 0;
        }

    	fEventTree->Fill();
    }
    else { mf::LogWarning("TrainingDataAlg") << "MVA FOR CLUSTERS NOT FOUND"; }

	cleanup(); // remove everything from members
}
/******************************************/

void nnet::PointIdEffTest::countTruthDep(
    const std::vector< sim::SimChannel > & channels,
    float & emLike, float & trackLike) const
{
    emLike = 0; trackLike = 0;
	for (auto const& channel : channels)
	{
		// for every time slice in this channel:
		auto const& timeSlices = channel.TDCIDEMap();
		for (auto const& timeSlice : timeSlices)
		{
			// loop over the energy deposits.
			auto const& energyDeposits = timeSlice.second;
			for (auto const& energyDeposit : energyDeposits)
			{
				int trackID = energyDeposit.trackID;

				double energy = energyDeposit.numElectrons * fElectronsToGeV * 1000;

				if (trackID < 0)
				{
					emLike += energy;
				}
				else if (trackID > 0)
				{
					auto search = fParticleMap.find(trackID);
					bool found = true;
					if (search == fParticleMap.end())
					{
						mf::LogWarning("TrainingDataAlg") << "PARTICLE NOT FOUND";
						found = false;
					}

					int pdg = 0;
					if (found)
					{
						const simb::MCParticle& particle = *((*search).second);
						if (!pdg) pdg = particle.PdgCode(); // not EM activity so read what PDG it is
					}

					if ((pdg == 11) || (pdg == -11) || (pdg == 22)) emLike += energy;
					else trackLike += energy;
				}
			}
		}
	}
}
/******************************************/

void nnet::PointIdEffTest::countPfpDep(
        const std::vector< recob::PFParticle > & pfparticles,
        const art::FindManyP<recob::Cluster> & pfpclus,
        const art::FindManyP<recob::Hit> & cluhits,
        float & emLike, float & trackLike) const
{
    emLike = 0; trackLike = 0;
    for (size_t i = 0; i < pfparticles.size(); ++i)
    {
        const auto & pfp = pfparticles[i];
        const auto & clus = pfpclus.at(i);

        float hitdep = 0;
        for (const auto & c : clus)
        {
            const auto & hits = cluhits.at(c.key());
            for (const auto & h : hits)
            {
                if (h->View() == fView) { hitdep += h->SummedADC() * fCalorimetryAlg.LifetimeCorrection(h->PeakTime()); }
            }
        }

        if ((pfp.PdgCode() == 11) || pfp.PdgCode() == -11) { emLike += hitdep; }
        else { trackLike += hitdep; }
    }
}
/******************************************/

bool nnet::PointIdEffTest::isMuonDecaying(const simb::MCParticle & particle,
    const std::unordered_map< int, const simb::MCParticle* > & particleMap) const
{
    bool hasElectron = false, hasNuMu = false, hasNuE = false;

    int pdg = abs(particle.PdgCode());
	if ((pdg == 13) && (particle.EndProcess() == "FastScintillation")) // potential muon decay at rest
	{
		unsigned int nSec = particle.NumberDaughters();
		for (size_t d = 0; d < nSec; ++d)
		{
			auto d_search = particleMap.find(particle.Daughter(d));
			if (d_search != particleMap.end())
			{
				auto const & daughter = *((*d_search).second);
				int d_pdg = abs(daughter.PdgCode());
				if (d_pdg == 11) hasElectron = true;
				else if (d_pdg == 14) hasNuMu = true;
				else if (d_pdg == 12) hasNuE = true;
			}
		}
	}

	return (hasElectron && hasNuMu && hasNuE);
}
/******************************************/

int nnet::PointIdEffTest::testCNN(
    const std::vector< sim::SimChannel > & channels,
    const std::vector< art::Ptr<recob::Hit> > & hits,
    const std::array<float, MVA_LENGTH> & cnn_out,
    const std::vector< anab::FeatureVector<MVA_LENGTH> > & hit_outs,
    size_t cidx)
{
	fClSize = hits.size();

    std::unordered_map<int, int> mcHitPid;

	fPidValue = 0;
	double p_trk_or_sh = cnn_out[fTrkLikeIdx] + cnn_out[fEmLikeIdx];
	if (p_trk_or_sh > 0) { fPidValue = cnn_out[fTrkLikeIdx] / p_trk_or_sh; }

    double p_michel = 0;
    if (fMichelLikeIdx >= 0) { fpMichel_cl = cnn_out[fMichelLikeIdx]; }

	double totEnSh = 0, totEnTrk = 0, totEnMichel = 0;
	for (auto const & hit: hits)
	{
		// the channel associated with this hit.
		auto hitChannelNumber = hit->Channel();

		double hitEn = 0, hitEnSh = 0, hitEnTrk = 0, hitEnMichel = 0;
        
		auto const & vout = hit_outs[hit.key()];
		fOutTrk = vout[fTrkLikeIdx];
		fOutEM = vout[fEmLikeIdx];
		if (fNoneIdx >= 0) { fOutNone = vout[fNoneIdx]; }
		if (fMichelLikeIdx >= 0) { p_michel = vout[fMichelLikeIdx];  }

		for (auto const & channel : channels)
		{
			if (channel.Channel() != hitChannelNumber) continue;

			// for every time slice in this channel:
			auto const& timeSlices = channel.TDCIDEMap();
			for (auto const& timeSlice : timeSlices)
			{
				int time = timeSlice.first;
				if (std::abs(hit->TimeDistanceAsRMS(time)) < 1.0)
				{
					// loop over the energy deposits.
					auto const & energyDeposits = timeSlice.second;
		
					for (auto const & energyDeposit : energyDeposits)
					{
						int trackID = energyDeposit.trackID;

						double energy = energyDeposit.numElectrons * fElectronsToGeV * 1000;
						hitEn += energy;

						if (trackID < 0) { hitEnSh += energy; } // EM activity
						else if (trackID > 0)
						{
							auto search = fParticleMap.find(trackID);
							if (search != fParticleMap.end())
							{
								const simb::MCParticle & particle = *((*search).second);
								int pdg = particle.PdgCode(); // not EM activity so read what PDG it is

							    if ((pdg == 11) || (pdg == -11) || (pdg == 22)) hitEnSh += energy;
							    else hitEnTrk += energy;

							    if (pdg == 11) // electron, check if it is Michel
                                {
                  	                auto msearch = fParticleMap.find(particle.Mother());
	    					        if (msearch != fParticleMap.end())
	    					        {
	    					    	    auto const & mother = *((*msearch).second);
	    		                        if (isMuonDecaying(mother, fParticleMap)) { hitEnMichel += energy; }
	    					        }
                                }
                            }
							else { mf::LogWarning("TrainingDataAlg") << "PARTICLE NOT FOUND"; }
						}
					}
				}
			}
		}
		totEnSh += hitEnSh;
		totEnTrk += hitEnTrk;
		totEnMichel += hitEnMichel;

        double hitAdc = hit->SummedADC() * fCalorimetryAlg.LifetimeCorrection(hit->PeakTime());
		fTotHit += hitAdc;

        int hitPidMc_0p5 = -1;
		if (hitEnSh > hitEnTrk)
		{
		    fHitEM_0p5 += hitAdc;
		    hitPidMc_0p5 = nnet::PointIdEffTest::kShower;
		}
		else
		{
		    fHitTrack_0p5 += hitAdc;
		    hitPidMc_0p5 = nnet::PointIdEffTest::kTrack;
		}
		mcHitPid[hit.key()] = hitPidMc_0p5;
		auto const & hout = hit_outs[hit.key()];
	    fpEM = 0;
        float hit_trk_or_sh = hout[fTrkLikeIdx] + hout[fEmLikeIdx];
        if (hit_trk_or_sh > 0) fpEM = hout[fEmLikeIdx] / hit_trk_or_sh;
		fHitEM_mc = hitEnSh / (hitEnSh + hitEnTrk);

        int hitPidMc_0p85 = -1;
        double hitDep = hitEnSh + hitEnTrk;
		if (hitEnSh > 0.85 * hitDep)
		{
		    fHitEM_0p85 += hitAdc; fCleanHit += hitAdc;
		    hitPidMc_0p85 = nnet::PointIdEffTest::kShower;
		}
		else if (hitEnTrk > 0.85 * hitDep)
		{
		    fHitTrack_0p85 += hitAdc; fCleanHit += hitAdc;
		    hitPidMc_0p85 = nnet::PointIdEffTest::kTrack;
		}


        bool mcMichel = false;
        fpMichel_hit = p_michel;
        fHitMichel_mc = hitEnMichel / hitEn;
        if (fHitMichel_mc > 0.5)
        {
            fHitMichel_0p5 += hitAdc; mcMichel = true;
        }

        fHitTree->Fill();

		for (size_t i = 0; i < 100; ++i)
		{
		    double thr = 0.01 * i;

            if (p_michel > thr)
            {
                if (mcMichel) { fHitsMichel_OK_0p5[i] += hitAdc; }
                else { fHitsMichel_False_0p5[i] += hitAdc; }
            }

            int recoPid = -1;
    	    if (fPidValue < thr)
    	    {
    	        recoPid = nnet::PointIdEffTest::kShower;
    	        fHitRecoEM[i] += hitAdc;
    	    }
	        else recoPid = nnet::PointIdEffTest::kTrack;

	        if ((recoPid == nnet::PointIdEffTest::kShower) && (hitPidMc_0p5 == nnet::PointIdEffTest::kShower))
	        {
	        	fHitsEM_OK_0p5[i] += hitAdc;
	        }
	        else if ((recoPid == nnet::PointIdEffTest::kTrack) && (hitPidMc_0p5 == nnet::PointIdEffTest::kTrack))
	        {
	        	fHitsTrack_OK_0p5[i] += hitAdc;
	        }

	        if ((recoPid == nnet::PointIdEffTest::kShower) && (hitPidMc_0p85 == nnet::PointIdEffTest::kShower))
	        {
	        	fHitsEM_OK_0p85[i] += hitAdc;
	        }
	        else if ((recoPid == nnet::PointIdEffTest::kTrack) && (hitPidMc_0p85 == nnet::PointIdEffTest::kTrack))
	        {
	        	fHitsTrack_OK_0p85[i] += hitAdc;
	        }
		}
	}

    // ************ count clusters *************

	fMcPid = -1;
	if (totEnSh > 1.5 * totEnTrk) // major energy deposit from EM activity
	{
		fMcPid = nnet::PointIdEffTest::kShower;
	}
	else if (totEnTrk > 1.5 * totEnSh)
	{
		fMcPid = nnet::PointIdEffTest::kTrack;
	}

    for (size_t i = 0; i < 100; ++i)
    {
        double thr = 0.01 * i;

        int recoPid = -1;
    	if (fPidValue < thr) recoPid = nnet::PointIdEffTest::kShower;
	    else recoPid = nnet::PointIdEffTest::kTrack;

	    if ((recoPid == nnet::PointIdEffTest::kShower) && (fMcPid == nnet::PointIdEffTest::kShower))
	    {
	    	fShOk[i] += fClSize;
	    }
	    else if ((recoPid == nnet::PointIdEffTest::kTrack) && (fMcPid == nnet::PointIdEffTest::kTrack))
	    {
	    	fTrkOk[i] += fClSize;
	    }	
	    else if ((recoPid == nnet::PointIdEffTest::kShower) && (fMcPid == nnet::PointIdEffTest::kTrack))
	    {
	    	fTrkBad[i] += fClSize;
	    }
	    else if ((recoPid == nnet::PointIdEffTest::kTrack) && (fMcPid == nnet::PointIdEffTest::kShower))
	    {
	    	fShBad[i] += fClSize;
	    }
	    else
	    {
	    	fNone++;
	    }

	}
	fTotal++;

	if (fSaveHitsFile)
	{
		for (auto const & h : hits)
		{
		    auto const & vout = hit_outs[h.key()];
	    	double hitPidValue = 0;
        	double h_trk_or_sh = vout[fTrkLikeIdx] + vout[fEmLikeIdx];
        	if (h_trk_or_sh > 0) hitPidValue = vout[fTrkLikeIdx] / h_trk_or_sh;

			fHitsOutFile << fRun << " " << fEvent << " "
				<< h->WireID().TPC  << " " << h->WireID().Wire << " " << h->PeakTime() << " "
				<< h->SummedADC() * fCalorimetryAlg.LifetimeCorrection(h->PeakTime()) << " "
				<< mcHitPid[h.key()] << " " << fPidValue << " " << hitPidValue;

			if (fMichelLikeIdx >= 0)
			{
			    fHitsOutFile << " " << vout[fMichelLikeIdx]; // is michel?
			}

			fHitsOutFile << " " << cidx << std::endl;
		}
	}

    fClusterTree->Fill();
	return fMcPid;
}
/******************************************/

DEFINE_ART_MODULE(nnet::PointIdEffTest)
