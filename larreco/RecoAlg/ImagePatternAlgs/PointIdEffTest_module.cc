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
#include "larcore/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
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

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/MVAReader.h"

#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>

#include <cmath>

namespace nnet {

class PointIdEffTest : public art::EDAnalyzer {
public:
    enum EId { kShower = 0, kTrack = 1 };

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<calo::CalorimetryAlg::Config> CalorimetryAlg {
			Name("CalorimetryAlg"),
			Comment("Used to calculate electron lifetime correction.")
		};

		fhicl::Atom<art::InputTag> SimModuleLabel {
			Name("SimModuleLabel"),
			Comment("...")
		};

		fhicl::Atom<art::InputTag> NNetModuleLabel {
			Name("NNetModuleLabel"),
			Comment("NNet outputs tag")
		};

		fhicl::Atom<bool> SaveHitsFile {
			Name("SaveHitsFile"),
			Comment("...")
		};

		fhicl::Atom<unsigned int> View {
			Name("View"),
			Comment("...")
		};
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

	int testCNN(
	    const std::vector< sim::SimChannel > & channels,
        const std::vector< art::Ptr<recob::Hit> > & hits,
        const std::array<float, 3> & cnn_out,
        const std::vector< anab::FeatureVector<3> > & hit_outs);

	int fRun, fEvent;
    float fMcDepEM, fMcDepTrack, fMcFractionEM;
    float fHitEM_0p5, fHitTrack_0p5, fHitMcFractionEM;
    float fHitEM_0p85, fHitTrack_0p85;
    float fTotHit, fCleanHit;

    //const int kEffSize = 100;
    float fHitsEM_OK_0p5[100], fHitsTrack_OK_0p5[100];
    float fHitsEM_OK_0p85[100], fHitsTrack_OK_0p85[100];
    float fHitRecoEM[100], fHitRecoFractionEM[100];

	int fMcPid;
	int fClSize;
	int fPure;
	double fPidValue;

	int fTrkOk[100], fTrkBad[100];
	int fShOk[100], fShBad[100];
	int fNone, fTotal;

	double fElectronsToGeV;

	TTree *fEventTree, *fClusterTree;

	std::ofstream fHitsOutFile;

	unsigned int fView;

	geo::GeometryCore const* fGeometry;

	std::unordered_map< int, const simb::MCParticle* > fParticleMap;

    calo::CalorimetryAlg fCalorimetryAlg;
	art::InputTag fSimulationProducerLabel;
	art::InputTag fNNetModuleLabel;
	bool fSaveHitsFile;
};

} // namespace nnet

nnet::PointIdEffTest::PointIdEffTest(nnet::PointIdEffTest::Parameters const& config) : art::EDAnalyzer(config),
	fMcPid(-1), fClSize(0),

	fView(config().View()),
	fCalorimetryAlg(config().CalorimetryAlg()),
	fSimulationProducerLabel(config().SimModuleLabel()),
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
    fEventTree->Branch("fHitEM_0p5", &fHitEM_0p5, "fHitEM_0p5/F");
    fEventTree->Branch("fHitTrack_0p5", &fHitTrack_0p5, "fHitTrack_0p5/F");
    fEventTree->Branch("fHitMcFractionEM", &fHitMcFractionEM, "fHitMcFractionEM/F");
    fEventTree->Branch("fHitEM_0p85", &fHitEM_0p85, "fHitEM_0p85/F");
    fEventTree->Branch("fHitTrack_0p85", &fHitTrack_0p85, "fHitTrack_0p85/F");
    fEventTree->Branch("fCleanHit", &fCleanHit, "fCleanHit/F");

    fEventTree->Branch("fHitsEM_OK_0p5", fHitsEM_OK_0p5, "fHitsEM_OK_0p5[100]/F");
    fEventTree->Branch("fHitsTrack_OK_0p5", fHitsTrack_OK_0p5, "fHitsTrack_OK_0p5[100]/F");

    fEventTree->Branch("fHitsEM_OK_0p85", fHitsEM_OK_0p85, "fHitsEM_OK_0p85[100]/F");
    fEventTree->Branch("fHitsTrack_OK_0p85", fHitsTrack_OK_0p85, "fHitsTrack_OK_0p85[100]/F");

    fEventTree->Branch("fHitRecoEM", fHitRecoEM, "fHitRecoEM[100]/F");
    fEventTree->Branch("fHitRecoFractionEM", fHitRecoFractionEM, "fHitRecoFractionEM[100]/F");


	fClusterTree = tfs->make<TTree>("cluster","clusters info");
	fClusterTree->Branch("fMcPid", &fMcPid, "fMcPid/I");
	fClusterTree->Branch("fClSize", &fClSize, "fClSize/I");
	fClusterTree->Branch("fPidValue", &fPidValue, "fPidValue/D");

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
	std::cout << "Total " << fTotal << std::endl;
}

void nnet::PointIdEffTest::cleanup(void)
{
    fParticleMap.clear();

    fMcDepEM = 0; fMcDepTrack = 0;
    fMcFractionEM = 0;
    fHitEM_0p5 = 0; fHitTrack_0p5 = 0;
    fHitEM_0p85 = 0; fHitTrack_0p85 = 0;
    fHitMcFractionEM = 0;
    fTotHit = 0; fCleanHit = 0;

    for (size_t i = 0; i < 100; ++i)
    {
        fHitsEM_OK_0p5[i] = 0; fHitsTrack_OK_0p5[i] = 0;
        fHitsEM_OK_0p85[i] = 0; fHitsTrack_OK_0p85[i] = 0;
        fHitRecoEM[i] = 0; fHitRecoFractionEM[i] = 0;
    }
}

void nnet::PointIdEffTest::analyze(art::Event const & e)
{
    cleanup(); // remove everything from members

	fRun = e.run();
	fEvent = e.id().event();
	std::cout << "event " << fEvent << std::endl;

	// access to MC information
	
	// MC particles list
	auto particleHandle = e.getValidHandle< std::vector<simb::MCParticle> >(fSimulationProducerLabel);
	for (auto const& particle : *particleHandle) { fParticleMap[particle.TrackId()] = &particle; }

	// SimChannels
	auto simChannelHandle = e.getValidHandle< std::vector<sim::SimChannel> >(fSimulationProducerLabel);

	countTruthDep(*simChannelHandle, fMcDepEM, fMcDepTrack);

    // output from cnn's

    anab::MVAReader<recob::Hit, 3> hitResults(e, fNNetModuleLabel);     // hit-by-hit outpus just dumped to file for debugging
    anab::MVAReader<recob::Cluster, 3> cluResults(e, fNNetModuleLabel); // outputs for clusters used as a EM/track decision
 
	const art::FindManyP<recob::Hit> hitsFromClusters(cluResults.dataHandle(), e, cluResults.dataTag());

	for (size_t c = 0; c < cluResults.size(); ++c)
	{
	    const recob::Cluster & clu = cluResults.item(c);

        if (clu.Plane().Plane != fView) continue;

	    const std::vector< art::Ptr<recob::Hit> > & hits = hitsFromClusters.at(c);
	    std::array<float, 3> cnn_out = cluResults.getOutput(c);

        testCNN(*simChannelHandle, hits, cnn_out, hitResults.outputs()); // test hits in the cluster
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

int nnet::PointIdEffTest::testCNN(
    const std::vector< sim::SimChannel > & channels,
    const std::vector< art::Ptr<recob::Hit> > & hits,
    const std::array<float, 3> & cnn_out,
    const std::vector< anab::FeatureVector<3> > & hit_outs)
{
	fClSize = hits.size();

	fPidValue = 0;
	double p_trk_or_sh = cnn_out[0] + cnn_out[1];
	if (p_trk_or_sh > 0) fPidValue = cnn_out[0] / p_trk_or_sh;

	double totEnSh = 0, totEnTrk = 0;
	for (auto const & hit: hits)
	{
		// the channel associated with this hit.
		auto hitChannelNumber = hit->Channel();

        double hitEnSh = 0, hitEnTrk = 0;
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

						if (trackID < 0)
						{
							hitEnSh += energy;
						}
						else if (trackID > 0)
						{
							auto search = fParticleMap.find(trackID);
							bool found = (search != fParticleMap.end());

							int pdg = 0;
							if (found)
							{
								const simb::MCParticle & particle = *((*search).second);
								if (!pdg) pdg = particle.PdgCode(); // not EM activity so read what PDG it is
							}
							else
							{
								mf::LogWarning("TrainingDataAlg") << "PARTICLE NOT FOUND";
							}

							if ((pdg == 11) || (pdg == -11) || (pdg == 22)) hitEnSh += energy;
							else hitEnTrk += energy;
						}
					}
				}
			}
		}
		totEnSh += hitEnSh;
		totEnTrk += hitEnTrk;

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

		for (size_t i = 0; i < 100; ++i)
		{
		    double thr = 0.01 * i;
		    
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
        	double h_trk_or_sh = vout[0] + vout[1];
        	if (h_trk_or_sh > 0) hitPidValue = vout[0] / h_trk_or_sh;

			fHitsOutFile << fRun << " " << fEvent << " "
				<< h->WireID().TPC  << " " << h->WireID().Wire << " " << h->PeakTime() << " "
				<< h->SummedADC() * fCalorimetryAlg.LifetimeCorrection(h->PeakTime()) << " "
				<< fMcPid << " " << fPidValue << " " << hitPidValue << std::endl;
		}
	}

    fClusterTree->Fill();
	return fMcPid;
}

/******************************************/

DEFINE_ART_MODULE(nnet::PointIdEffTest)
