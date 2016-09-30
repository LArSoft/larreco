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
#include "larcore/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"

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
			Comment("...")
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
			Comment("...")
		};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit PointIdTrainingData(Parameters const& config);

    virtual void analyze (const art::Event& event) override;

  private:

    nnet::TrainingDataAlg fTrainingDataAlg;

	std::string fOutTextFilePath;

	std::vector<int> fSelectedTPC;
	std::vector<int> fSelectedView;

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
	fSelectedTPC(config().SelectedTPC()),
	fSelectedView(config().SelectedView()),
	fCrop(config().Crop())
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());

	const size_t TPC_CNT = (size_t)fGeometry->NTPC(0);
	if (fSelectedTPC.empty())
	{
		for (size_t tpc = 0; tpc < TPC_CNT; ++tpc)
			fSelectedTPC.push_back(tpc);
	}

	if (fSelectedView.empty())
	{
		if (fGeometry->TPC(fSelectedTPC.front(), 0).HasPlane(geo::kU))
			fSelectedView.push_back((int)geo::kU);
		if (fGeometry->TPC(fSelectedTPC.front(), 0).HasPlane(geo::kV))
			fSelectedView.push_back((int)geo::kV);
		if (fGeometry->TPC(fSelectedTPC.front(), 0).HasPlane(geo::kZ))
			fSelectedView.push_back((int)geo::kZ);
	}
  }

  //-----------------------------------------------------------------------
  void PointIdTrainingData::analyze(const art::Event& event) 
  {
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

	std::ostringstream os;
	os << "event_" << fEvent << "_run_" << fRun << "_subrun_" << fSubRun;

	std::cout << "analyze " << os.str() << std::endl;

	std::ofstream fout_raw, fout_deposit, fout_pdg;

	for (size_t i = 0; i < fSelectedTPC.size(); ++i)
		for (size_t v = 0; v < fSelectedView.size(); ++v)
	{
		fTrainingDataAlg.setEventData(event, fSelectedView[v], fSelectedTPC[i], 0);

        unsigned int w0, w1, d0, d1;
        if (fCrop)
        {
            if (fTrainingDataAlg.findCrop(0.004F, w0, w1, d0, d1))
            {
                std::cout << "   crop: " << w0 << " " << w1 << " " << d0 << " " << d1 << std::endl;
            }
            else
            {
                std::cout << "   skip empty tpc:" << fSelectedTPC[i] << " / view:" << fSelectedView[v] << std::endl;
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

		std::ostringstream ss1;
		ss1 << fOutTextFilePath << "/raw_" << os.str()
			<< "_tpc_" << fSelectedTPC[i]
			<< "_view_" << fSelectedView[v];

		fout_raw.open(ss1.str() + ".raw");
		fout_deposit.open(ss1.str() + ".deposit");
		fout_pdg.open(ss1.str() + ".pdg");

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

		fout_raw.close();
		fout_deposit.close();
		fout_pdg.close();
	}

  } // Raw2DRegionID::analyze()

  DEFINE_ART_MODULE(PointIdTrainingData)

}

#endif // PointIdTrainingData_Module

