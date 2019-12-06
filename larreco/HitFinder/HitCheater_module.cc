////////////////////////////////////////////////////////////////////////
// Class:       HitCheater
// Module Type: producer
// File:        HitCheater_module.cc
//
// Generated at Tue Nov  8 09:41:20 2011 by Brian Rebel using artmod
// from art v1_00_02.
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <cmath> // std::sqrt()
#include <string>
#include <map>

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/Utilities/StatCollector.h"

namespace hit {
  class HitCheater;
}

class hit::HitCheater : public art::EDProducer {
public:
  explicit HitCheater(fhicl::ParameterSet const & p);

private:
  void produce(art::Event & e) override;

  void FindHitsOnChannel(
    const sim::SimChannel*   sc,
    std::vector<recob::Hit>& hits,
    int                      spill
    );


  std::string         fG4ModuleLabel;              ///< label name for module making sim::SimChannels
  std::string         fWireModuleLabel;                 ///< label name for module making recob::Wires
  double              fMinCharge;                       ///< Minimum charge required to make a hit
  double              fElectronsToADC;             ///< Conversion factor of electrons to ADC counts
  std::string         fCalDataProductInstanceName; ///< label name for module making recob::Wires
  int                 fReadOutWindowSize;          ///< Number of samples in a readout window; NOT total samples
  int                 fNumberTimeSamples;          ///< Number of total time samples (N*readoutwindowsize)
  double              fSamplingRate;               ///< from detinfo::DetectorPropertiesService
  int                 fTriggerOffset;              ///< from detinfo::DetectorPropertiesService
  int                 fNewHitTDCGap;               ///< gap allowed in tdcs without charge before making a new hit
};

//-------------------------------------------------------------------
hit::HitCheater::HitCheater(fhicl::ParameterSet const & p)
  : EDProducer{p}
{
  fG4ModuleLabel   = p.get< std::string >("G4ModuleLabel",   "largeant");
  fWireModuleLabel = p.get< std::string >("WireModuleLabel", "caldata" );
  fMinCharge       = p.get< double      >("MinimumCharge",   5.        );
  fNewHitTDCGap    = p.get< int         >("NewHitTDCGap",    1         );

  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fElectronsToADC = detprop->ElectronsToADC();
  fSamplingRate   = detprop->SamplingRate();
  fTriggerOffset  = detprop->TriggerOffset();

   fCalDataProductInstanceName="";
   size_t pos = fWireModuleLabel.find(":");
   if( pos!=std::string::npos ) {
     fCalDataProductInstanceName = fWireModuleLabel.substr( pos+1 );
     fWireModuleLabel = fWireModuleLabel.substr( 0, pos );
   }

   fReadOutWindowSize  = detprop->ReadOutWindowSize();
   fNumberTimeSamples  = detprop->NumberTimeSamples();

  // let HitCollectionCreator declare that we are going to produce
  // hits and associations with wires and raw digits
  // (with no particular product label)
   recob::HitCollectionCreator::declare_products(producesCollector());
}

//-------------------------------------------------------------------
void hit::HitCheater::produce(art::Event & e)
{
  // this object contains the hit collection
  // and its associations to wires and raw digits:
  recob::HitCollectionCreator hits(e);

  // Read in the wire List object(s).
  art::InputTag WireInputTag(fWireModuleLabel, fCalDataProductInstanceName);
  art::ValidHandle< std::vector<recob::Wire>> wHandle =
    e.getValidHandle<std::vector<recob::Wire>>(WireInputTag);

  int whatSpill = 1;
  if( !fCalDataProductInstanceName.empty() ) {
    if( fCalDataProductInstanceName.find("ost") != std::string::npos) whatSpill=2;
    else whatSpill=0;
  }

  // also get the raw digits associated with the wires;
  // we assume they have been created by the same module as the wires
  art::FindOneP<raw::RawDigit> WireToRawDigits(wHandle, e, WireInputTag);

  // make a map of wires to channel numbers
  std::map<raw::ChannelID_t, art::Ptr<recob::Wire>> wireMap;

  for(size_t wc = 0; wc < wHandle->size(); ++wc){
    art::Ptr<recob::Wire> wire(wHandle, wc);
    wireMap[wire->Channel()] = wire;
  }

  // get the sim::SimChannels out of the event
  std::vector<const sim::SimChannel*> sccol;
  e.getView(fG4ModuleLabel, sccol);

  // find the hits on each channel
  for(sim::SimChannel const* sc: sccol) {
    std::vector<recob::Hit> hits_on_channel;

    FindHitsOnChannel(sc, hits_on_channel, whatSpill);

    art::Ptr<recob::Wire> const& wire = wireMap[sc->Channel()];
    art::Ptr<raw::RawDigit> rawdigits; // null by default
    if (wire.isNonnull()) rawdigits = WireToRawDigits.at(wire.key());

    // add all the hits found on this channel to the data product,
    // all associated to the same hit and wire
    for (recob::Hit& hit: hits_on_channel)
      hits.emplace_back(std::move(hit), wire, rawdigits);

  }// end loop over SimChannels

  // put the cheated hits into the event
  MF_LOG_DEBUG("HitCheater") << "putting " << hits.size() << " hits into the event";
  hits.put_into(e);

  return;
}

//-------------------------------------------------------------------
void hit::HitCheater::FindHitsOnChannel(const sim::SimChannel*   sc,
                                        std::vector<recob::Hit>& hits,
                                        int                      spill)
{
  art::ServiceHandle<geo::Geometry const> geo;

  raw::ChannelID_t channel = sc->Channel();
  geo::SigType_t signal_type = geo->SignalType(channel);
  geo::View_t view = geo->View(channel);


  // determine the possible geo::WireIDs for this particular channel
  // then make a map of tdc to electrons for each one of those geo::WireIDs
  // then find hits on each geo::WireID
  std::vector<geo::WireID> wireids = geo->ChannelToWire(channel);

  std::map<geo::WireID, std::map< unsigned int, double> > wireIDSignals;

  auto const& idemap = sc->TDCIDEMap();

  for(auto const& mapitr : idemap){
    unsigned short tdc = mapitr.first;

    if( fReadOutWindowSize != fNumberTimeSamples ) {
      if( tdc < spill*fReadOutWindowSize ||
          tdc > (spill+1)*fReadOutWindowSize )  continue;
     } else {
      if ( tdc < 0 || tdc > fReadOutWindowSize) continue;
    }

    // figure out which TPC we are in for each voxel

    for(auto const& ideitr : mapitr.second){

      const float edep = ideitr.numElectrons;

      std::array<double, 3> pos;
      pos[0] = ideitr.x;
      pos[1] = ideitr.y;
      pos[2] = ideitr.z;

      geo::TPCID tpcID = geo->FindTPCAtPosition(pos.data());
      if (!tpcID.isValid) {
        mf::LogWarning("HitCheater") << "TPC for position ( "
          << pos[0] << " ; " << pos[1] << " ; " << pos[2] << " )"
          << " in no TPC; move on to the next sim::IDE";
        continue;
      }
      const unsigned int tpc   = tpcID.TPC;
      const unsigned int cstat = tpcID.Cryostat;

      for( auto const& wid : wireids){
        if(wid.TPC == tpc && wid.Cryostat == cstat){
          // in the right TPC, now figure out which wire we want
          // this works because there is only one plane option for
          // each WireID in each TPC
          if(wid.Wire == geo->NearestWire(pos.data(), wid.Plane, wid.TPC, wid.Cryostat))
            wireIDSignals[wid][tdc] += edep;
        }// end if in the correct TPC and Cryostat
      }// end loop over wireids for this channel
    }// end loop over ides for this channel
  }// end loop over tdcs for this channel

  // now loop over each wire ID and determine where the hits are
  for( auto const& widitr : wireIDSignals){

    // get the first tdc in the
    unsigned short prev         = widitr.second.begin()->first;
    unsigned short startTime    = prev;
    double         totCharge    = 0.;
    double         maxCharge    = -1.;
    double         peakTime     = 0.;
    int            multiplicity =  1 ;
    lar::util::StatCollector<double> time; // reduce rounding errors

    // loop over all the tdcs for this geo::WireID
    for( auto tdcitr : widitr.second){
      unsigned short tdc = tdcitr.first;
      if (tdc < prev) {
        throw art::Exception(art::errors::LogicError)
          << "SimChannel TDCs going backward!";
      }

      // more than a one tdc gap between times with
      // signal, start a new hit
      if(tdc - prev > fNewHitTDCGap){

        if(totCharge > fMinCharge){
          hits.emplace_back(
            channel,       // channel
            (raw::TDCtick_t) startTime,       // start_tick
            (raw::TDCtick_t) prev,            // end_tick
            peakTime,                         // peak_time
            1.,                               // sigma_peak_time
            time.RMS(),                       // RMS
            maxCharge,                        // peak_amplitude
            std::sqrt(maxCharge),             // sigma_peak_amplitude
            totCharge,                        // summedADC
            totCharge,                        // hit_integral
            std::sqrt(totCharge),             // hit_sigma_integral
            multiplicity,                     // multiplicity
            0,                                // local_index
            1.,                               // goodness_of_fit
            0.,                               // dof
            view,                             // view
            signal_type,                      // signal_type
            widitr.first                      // wireID
            );

          MF_LOG_DEBUG("HitCheater") << "new hit is " << hits.back();

        }// end if charge is large enough

        // reset the variables for each hit
        startTime = tdc;
        peakTime  = tdc;
        totCharge = 0.;
        maxCharge = -1.;
        time.clear();

      }// end if need to start a new hit

      const double adc = fElectronsToADC*tdcitr.second;

      totCharge += adc;
      // use ADC as weight; the shift reduces the precision needed;
      // average would need to be shifted back: time.Averatge() + startTime
      // but RMS is not affected
      time.add(tdc - startTime, adc);
      if(adc > maxCharge){
        maxCharge = adc;
        peakTime = tdc;
      }

      prev = tdc;

    }// end loop over tdc values for the current geo::WireID


    // We might have missed the last hit, so do it now
    if(totCharge > fMinCharge){
      hits.emplace_back(
        channel,       // channel
        (raw::TDCtick_t) startTime,       // start_tick
        (raw::TDCtick_t) prev + 1,        // end_tick; prev included in the hit
        peakTime,                         // peak_time
        1.,                               // sigma_peak_time
        time.RMS(),                       // RMS
        maxCharge,                        // peak_amplitude
        std::sqrt(maxCharge),             // sigma_peak_amplitude
        totCharge,                        // summedADC
        totCharge,                        // hit_integral
        std::sqrt(totCharge),             // hit_sigma_integral
        multiplicity,                     // multiplicity
        0,                                // local_index
        1.,                               // goodness_of_fit
        0.,                               // dof
        view,                             // view
        signal_type,                      // signal_type
        widitr.first                      // wireID
        );

      MF_LOG_DEBUG("HitCheater") << "last hit is " << hits.back();

    }// end if charge is large enough

  }// end loop over map of geo::WireID to map<tdc,electrons>


  return;
}

DEFINE_ART_MODULE(hit::HitCheater)
