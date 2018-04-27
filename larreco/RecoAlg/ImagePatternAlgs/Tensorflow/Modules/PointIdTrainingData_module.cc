////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdTrainingData
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan, May 2016
//
// Training data for PointIdAlg
//
//      We use this to dump deconv. ADC for preparation of various classifiers.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PointIdTrainingData_Module
#define PointIdTrainingData_Module

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "larreco/RecoAlg/ImagePatternAlgs/Tensorflow/PointIdAlg/PointIdAlg.h"

// Framework includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
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

#include "TH2I.h" // PDG+vertex info map
#include "TH2F.h" // ADC and deposit maps

namespace nnet	 {

  class PointIdTrainingData : public art::EDAnalyzer
  {
  public:
 
 	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<nnet::TrainingDataAlg::Config> TrainingDataAlg {
			Name("TrainingDataAlg")
		};

		fhicl::Atom<std::string> OutTextFilePath {
			Name("OutTextFilePath"),
			Comment("Text files with all needed data dumped.")
		};

		fhicl::Atom<bool> DumpToRoot {
			Name("DumpToRoot"),
			Comment("Dump to ROOT histogram file (replaces the text files)")
		};

		fhicl::Sequence<int> SelectedTPC {
			Name("SelectedTPC"),
			Comment("use selected views only, or all views if empty list")
		};

		fhicl::Sequence<int> SelectedView {
			Name("SelectedView"),
			Comment("use selected tpc's only, or all tpc's if empty list")
		};

		fhicl::Atom<bool> Crop {
			Name("Crop"),
			Comment("Crop the projection to the event region plus margin")
		};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit PointIdTrainingData(Parameters const& config);

    void analyze(const art::Event& event) override;

  private:

    nnet::TrainingDataAlg fTrainingDataAlg;

	std::string fOutTextFilePath;
	bool fDumpToRoot;

	std::vector<int> fSelectedTPC;
	std::vector<int> fSelectedPlane;

	int fEvent;     /// number of the event being processed
	int fRun;       /// number of the run being processed
	int fSubRun;    /// number of the sub-run being processed

    bool fCrop;     /// crop data to event (set to false when dumping noise!)

	geo::GeometryCore const* fGeometry;
  };

  //-----------------------------------------------------------------------
  PointIdTrainingData::PointIdTrainingData(PointIdTrainingData::Parameters const& config) : art::EDAnalyzer(config),
	fTrainingDataAlg(config().TrainingDataAlg()),
	fOutTextFilePath(config().OutTextFilePath()),
	fDumpToRoot(config().DumpToRoot()),
	fSelectedTPC(config().SelectedTPC()),
	fSelectedPlane(config().SelectedView()),
	fCrop(config().Crop())
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());

	const size_t TPC_CNT = (size_t)fGeometry->NTPC(0);
	if (fSelectedTPC.empty())
	{
		for (size_t tpc = 0; tpc < TPC_CNT; ++tpc)
			fSelectedTPC.push_back(tpc);
	}

	if (fSelectedPlane.empty())
	{
	    for (size_t p = 0; p < fGeometry->MaxPlanes(); ++p )
	        fSelectedPlane.push_back(p);
	}
  }

  //-----------------------------------------------------------------------
  void PointIdTrainingData::analyze(const art::Event& event) 
  {
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    bool saveSim = fTrainingDataAlg.saveSimInfo() && !event.isRealData();

	std::ostringstream os;
	os << "event_" << fEvent << "_run_" << fRun << "_subrun_" << fSubRun;

	std::cout << "analyze " << os.str() << std::endl;

	for (size_t i = 0; i < fSelectedTPC.size(); ++i)
		for (size_t v = 0; v < fSelectedPlane.size(); ++v)
	{
		fTrainingDataAlg.setEventData(event, fSelectedPlane[v], fSelectedTPC[i], 0);

        unsigned int w0, w1, d0, d1;
        if (fCrop && saveSim)
        {
            if (fTrainingDataAlg.findCrop(0.004F, w0, w1, d0, d1))
            {
                std::cout << "   crop: " << w0 << " " << w1 << " " << d0 << " " << d1 << std::endl;
            }
            else
            {
                std::cout << "   skip empty tpc:" << fSelectedTPC[i] << " / view:" << fSelectedPlane[v] << std::endl;
                continue;
            }
        }
        else
        {
            w0 = 0; w1 = fTrainingDataAlg.NWires();
            d0 = 0; d1 = fTrainingDataAlg.NScaledDrifts();
        }

        if (fDumpToRoot)
        {
	        std::ostringstream ss1;
	        ss1 << "raw_" << os.str() << "_tpc_" << fSelectedTPC[i] << "_view_" << fSelectedPlane[v]; // TH2's name

            art::ServiceHandle<art::TFileService> tfs;
            TH2F* rawHist = tfs->make<TH2F>((ss1.str() + "_raw").c_str(), "ADC", w1 - w0, w0, w1, d1 - d0, d0, d1);
            TH2F* depHist = 0;
            TH2I* pdgHist = 0;
            if (saveSim)
            {
                depHist = tfs->make<TH2F>((ss1.str() + "_deposit").c_str(), "Deposit", w1 - w0, w0, w1, d1 - d0, d0, d1);
                pdgHist = tfs->make<TH2I>((ss1.str() + "_pdg").c_str(), "PDG", w1 - w0, w0, w1, d1 - d0, d0, d1);
            }

            for (size_t w = w0; w < w1; ++w)
            {
                auto const & raw = fTrainingDataAlg.wireData(w);
                for (size_t d = d0; d < d1; ++d) { rawHist->Fill(w, d, raw[d]); }

                if (saveSim)
                {
    		        auto const & edep = fTrainingDataAlg.wireEdep(w);
	    	        for (size_t d = d0; d < d1; ++d) { depHist->Fill(w, d, edep[d]); }

	    	        auto const & pdg = fTrainingDataAlg.wirePdg(w);
	    	        for (size_t d = d0; d < d1; ++d) { pdgHist->Fill(w, d, pdg[d]); }
	    	    }
            }
        }
        else
        {
	        std::ostringstream ss1;
	        ss1 << fOutTextFilePath << "/raw_" << os.str()
        		<< "_tpc_" << fSelectedTPC[i] << "_view_" << fSelectedPlane[v];

            std::ofstream fout_raw, fout_deposit, fout_pdg;

        	fout_raw.open(ss1.str() + ".raw");
        	if (saveSim)
        	{
    	        fout_deposit.open(ss1.str() + ".deposit");
	            fout_pdg.open(ss1.str() + ".pdg");
	        }

	        for (size_t w = w0; w < w1; ++w)
	        {
		        auto const & raw = fTrainingDataAlg.wireData(w);
		        for (size_t d = d0; d < d1; ++d) { fout_raw << raw[d] << " "; }
		        fout_raw << std::endl;

                if (saveSim)
                {
    		        auto const & edep = fTrainingDataAlg.wireEdep(w);
	    	        for (size_t d = d0; d < d1; ++d) { fout_deposit << edep[d] << " "; }
	    	        fout_deposit << std::endl;

		            auto const & pdg = fTrainingDataAlg.wirePdg(w);
		            for (size_t d = d0; d < d1; ++d) { fout_pdg << pdg[d] << " "; }
		            fout_pdg << std::endl;
		        }
	        }

	        fout_raw.close();
	        if (saveSim)
	        {
    	        fout_deposit.close();
	            fout_pdg.close();
	        }
        }
	}

  } // PointIdTrainingData::analyze()

  DEFINE_ART_MODULE(PointIdTrainingData)

}

#endif // PointIdTrainingData_Module

