////////////////////////////////////////////////////////////////////////
// Class:       MCHitFinder
// Module Type: producer
// File:        MCHitFinder_module.cc
//
// Generated at Tue Jun 24 07:30:08 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>

#include "lardataobj/Simulation/SimChannel.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/MCBase/MCDataHolder.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "lardataobj/MCBase/MCWireCollection.h"

namespace hit {

  class MCHitFinder;

  class MCHitFinder : public art::EDProducer {
  public:
    explicit MCHitFinder(fhicl::ParameterSet const& p);

  private:
    void produce(art::Event& e) override;

    std::string fLArG4ModuleName;
    bool fVerbose;
    bool fMakeMCWire;
  };

  MCHitFinder::MCHitFinder(fhicl::ParameterSet const& p) : EDProducer{p}
  {
    fLArG4ModuleName = p.get<std::string>("LArG4ModuleName");
    fMakeMCWire = p.get<bool>("MakeMCWire");
    fVerbose = p.get<bool>("Verbose");
    produces<std::vector<sim::MCHitCollection>>();
    if (fMakeMCWire) produces<std::vector<sim::MCWireCollection>>();
  }

  void MCHitFinder::produce(art::Event& e)
  {

    art::ServiceHandle<geo::Geometry const> geo;

    const unsigned int nch = geo->Nchannels();

    std::unique_ptr<std::vector<sim::MCHitCollection>> hits_v(
      new std::vector<sim::MCHitCollection>());
    std::unique_ptr<std::vector<sim::MCWireCollection>> wires_v(
      new std::vector<sim::MCWireCollection>());

    hits_v->reserve(nch);
    wires_v->reserve(nch);
    for (size_t ch = 0; ch < nch; ++ch) {

      hits_v->push_back(sim::MCHitCollection(ch));
      wires_v->push_back(sim::MCWireCollection(ch));
    }

    art::Handle<std::vector<sim::SimChannel>> simchArray;
    e.getByLabel(fLArG4ModuleName, simchArray);
    if (!simchArray.isValid())
      throw cet::exception(__PRETTY_FUNCTION__)
        << "Did not find sim::SimChannel with a label: " << fLArG4ModuleName.c_str() << std::endl;

    // Loop over SimChannel
    for (size_t simch_index = 0; simch_index < simchArray->size(); ++simch_index) {

      const art::Ptr<sim::SimChannel> simch_ptr(simchArray, simch_index);

      size_t ch = simch_ptr->Channel();

      if (ch >= hits_v->size())

        throw cet::exception(__PRETTY_FUNCTION__)
          << "Channel number " << ch << " exceeds total # of channels: " << nch << std::endl;

      auto& mchits = hits_v->at(ch);
      auto& mcwires = wires_v->at(ch);

      std::map<sim::MCEnDep, sim::MCWire> edep_wire_map;
      auto tdc_ide_map = simch_ptr->TDCIDEMap();

      ////////////////////////////////////////////////////////////
      // Loop over stored data for a particular sim::SimChannel //
      ////////////////////////////////////////////////////////////

      //
      // Read & Convert
      //
      if (fVerbose) std::cout << std::endl << "Processing Ch: " << ch << std::endl;

      for (auto const& tdc_ide_pair : tdc_ide_map) {

        auto const& tdc = tdc_ide_pair.first;
        auto const& ide_v = tdc_ide_pair.second;

        for (auto const& ide : ide_v) {
          sim::MCEnDep edep;
          edep.SetVertex(ide.x, ide.y, ide.z);
          edep.SetEnergy(ide.energy);
          edep.SetTrackId(ide.trackID);

          sim::MCWire wire;
          wire.SetStartTDC(tdc);

          auto edep_iter = edep_wire_map.insert(std::make_pair(edep, wire));

          //auto edep_iter = edep_wire_map.find(edep);

          auto last_tdc =
            (edep_iter).first->second.StartTDC() + (edep_iter).first->second.size() - 1;

          if (fVerbose) {

            if (edep_iter.second) std::cout << std::endl;

            std::cout << "  Track: " << ide.trackID << " Vtx: " << ide.x << " " << ide.y << " "
                      << ide.z << " " << ide.energy << " ...  @ TDC = " << tdc << " ... "
                      << ide.numElectrons << std::endl;
          }

          if (!(edep_iter.second)) {

            if (last_tdc + 1 != tdc) {
              /*
              std::cerr
                << "Found discontinuous TDC! "
                << " Last (ADC @ TDC): " << (*(edep_iter.first->second.rbegin())) << " @ " << last_tdc
                << " while current (ADC @ TDC): " << (ide.numElectrons * detp->ElectronsToADC()) << " @ " << tdc
                << " ... skipping to store!"
                << std::endl;
              */
              continue;
            }
          }
          (edep_iter).first->second.push_back(ide.numElectrons);
        } // Done looping over IDEs in a particular TDC
      }   // Done looping over TDC

      //
      // Store
      //
      for (auto const& edep_wire_pair : edep_wire_map) {

        auto const& edep = edep_wire_pair.first;
        auto const& wire = edep_wire_pair.second;

        // Create MCHit
        ::sim::MCHit hit;
        float vtx[3] = {float(edep.Vertex()[0]), float(edep.Vertex()[1]), float(edep.Vertex()[2])};
        hit.SetParticleInfo(vtx, edep.Energy(), edep.TrackId());

        double max_time = 0;
        double qsum = 0;
        double qmax = 0;
        for (size_t wire_index = 0; wire_index < wire.size(); ++wire_index) {

          auto q = wire.at(wire_index);

          qsum += q;

          if (q > qmax) {
            qmax = q;
            max_time = wire.StartTDC() + wire_index;
          }
        }
        hit.SetCharge(qsum, qmax);
        hit.SetTime(max_time, 0);

        mchits.push_back(hit);

        if (!fMakeMCWire) continue;

        mcwires.push_back(wire);
      } // End looping over all MCEnDep-MCWire pairs
    }   // End looping over all SimChannels

    simchArray.removeProduct();

    std::sort((*hits_v).begin(), (*hits_v).end());
    e.put(std::move(hits_v));

    if (fMakeMCWire) {
      std::sort((*wires_v).begin(), (*wires_v).end());
      e.put(std::move(wires_v));
    }
  }
}
DEFINE_ART_MODULE(hit::MCHitFinder)
