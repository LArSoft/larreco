/*!
 * Title:   TTHitFinder class
 * Author:  wketchum@lanl.gov
 * Inputs:  recob::Wire (calibrated)
 * Outputs: recob::Hit
 *
 * Description:
 * This module, TimeTickHitFinder (or TTHitFinder for short) is designed to
 * produce a minimal hit object, that is simple a time tick above threshold.
 * There is intention to allow for overlap of hits, with a downstream app
 * that will need to clean it up.
 */
#include <string>
#include <math.h>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindOneP.h"

// LArSoft Includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/HitCreator.h"

namespace hit{

  class TTHitFinder : public art::EDProducer {

  public:

    explicit TTHitFinder(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;

    std::string    fCalDataModuleLabel; /// Input caldata module name
    float          fMinSigPeakInd;      /// Induction wire signal height threshold at peak
    float          fMinSigPeakCol;      /// Collection wire signal height threshold at peak
    float          fMinSigTailInd;      /// Induction wire signal height threshold outside peak
    float          fMinSigTailCol;      /// Collection wire signal height threshold outside peak
    int            fIndWidth;           /// Induction wire hit width (in time ticks)
    int            fColWidth;           /// Collection wire hit width (in time ticks)

    float getTotalCharge(const float*,int,float);

  }; // class TTHitFinder

  //-------------------------------------------------
  TTHitFinder::TTHitFinder(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
  {
    fCalDataModuleLabel = pset.get< std::string >("CalDataModuleLabel");
    fMinSigPeakInd      = pset.get< float       >("MinSigPeakInd");
    fMinSigPeakCol      = pset.get< float       >("MinSigPeakCol");
    fMinSigTailInd      = pset.get< float       >("MinSigTailInd",-99); //defaulting to well-below zero
    fMinSigTailCol      = pset.get< float       >("MinSigTailCol",-99); //defaulting to well-below zero
    fIndWidth           = pset.get< int         >("IndWidth", 3); //defaulting to 3
    fColWidth           = pset.get< int         >("ColWidth", 3); //defaulting to 3

    //enforce a minimum width
    if(fIndWidth<1){
      mf::LogWarning("TTHitFinder") << "IndWidth must be 1 at minimum. Resetting width to one time tick";
      fIndWidth = 1;
    }
    if(fColWidth<1){
      mf::LogWarning("TTHitFinder") << "ColWidth must be 1 at minimum. Resetting width to one time tick";
      fColWidth = 1;
    }

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    recob::HitCollectionCreator::declare_products(producesCollector(), "uhits");
    recob::HitCollectionCreator::declare_products(producesCollector(), "vhits");
    recob::HitCollectionCreator::declare_products(producesCollector(), "yhits");

  }

  //-------------------------------------------------
  void TTHitFinder::produce(art::Event& evt)
  {

    // these objects contain the hit collections
    // and their associations to wires and raw digits:
    recob::HitCollectionCreator hitCollection_U(evt, "uhits");
    recob::HitCollectionCreator hitCollection_V(evt, "vhits");
    recob::HitCollectionCreator hitCollection_Y(evt, "yhits");

    // Read in the wire List object(s).
    art::Handle< std::vector<recob::Wire> > wireVecHandle;
    evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
    std::vector<recob::Wire> const& wireVec(*wireVecHandle);

    // also get the raw digits associated with wires
    art::FindOneP<raw::RawDigit> WireToRawDigits
      (wireVecHandle, evt, fCalDataModuleLabel);

    art::ServiceHandle<geo::Geometry const> geom;

    //initialize some variables that will be in the loop.
    float threshold_peak = 0;
    float threshold_tail = -99;
    int   width = 3;

    //Loop over wires
    for(unsigned int wireIter = 0; wireIter < wireVec.size(); wireIter++) {

      //get our wire
      art::Ptr<recob::Wire> wire(wireVecHandle, wireIter);
      art::Ptr<raw::RawDigit> const& rawdigits = WireToRawDigits.at(wireIter);

      std::vector<float> signal(wire->Signal());
      std::vector<float>::iterator timeIter;   // iterator for time bins
      geo::WireID wire_id = (geom->ChannelToWire(wire->Channel())).at(0); //just grabbing the first one


      //set the thresholds and widths based on wire type
      geo::SigType_t sigType = geom->SignalType(wire->Channel());
      if(sigType == geo::kInduction){
	threshold_peak = fMinSigPeakInd;
	threshold_tail = fMinSigTailInd;
	width          = fIndWidth;
      }
      else if(sigType == geo::kCollection){
	threshold_peak = fMinSigPeakCol;
	threshold_tail = fMinSigTailCol;
	width          = fColWidth;
      }

      //make a half_width variable to be the search window around each time tick.
      float half_width = ((float)width-1)/2.;

      //now do the loop over the time ticks on the wire
      int time_bin = -1;
      float peak_val = 0;
      for(timeIter = signal.begin(); timeIter < signal.end(); timeIter++){
	time_bin++;

	//set the peak value, taking average between ticks if desired total width is even
	if(width%2==1) peak_val = *timeIter;
	else if(width%2==0) peak_val = 0.5 * (*timeIter + *(timeIter+1));

	//continue immediately if we are not above the threshold
	if(peak_val < threshold_peak) continue;

	//continue if we are too close to the edge
	if( time_bin-half_width < 0 ) continue;
	if( time_bin+half_width > signal.size() ) continue;

	//if necessary, do loop over hit width, and check tail thresholds
	int begin_tail_tick = std::floor(time_bin-half_width);
	float totalCharge = getTotalCharge(&signal[begin_tail_tick],width,threshold_tail);
	if(totalCharge==-999) {
	  MF_LOG_DEBUG("TTHitFinder") << "Rejecting would be hit at (plane,wire,time_bin,first_bin,last_bin)=("
				   << wire_id.Plane << "," << wire_id.Wire << "," << time_bin << "," << begin_tail_tick << "," << begin_tail_tick+width-1 << "): "
				   << signal.at(time_bin-1) << " "
				   << signal.at(time_bin) << " "
				   << signal.at(time_bin+1);
	  continue;
	}

	//OK, if we've passed all tests up to this point, we have a hit!

	float hit_time = time_bin;
	if(width%2==0) hit_time = time_bin+0.5;

	// hit time region is 2 widths (4 RMS) wide
	const raw::TDCtick_t start_tick = hit_time - width,
	  end_tick = hit_time + width;

	// make the hit
	recob::HitCreator hit(
	  *wire,                     // wire
	  wire_id,                   // wireID
	  start_tick,                // start_tick
	  end_tick,                  // end_tick
	  width / 2.,                // rms
	  hit_time,                  // peak_time
	  0.,                        // sigma_peak_time
	  peak_val,                  // peak_amplitude
	  0.,                        // sigma_peak_amplitude
	  totalCharge,               // hit_integral
	  0.,                        // hit_sigma_integral
	  totalCharge,               // summedADC
	  1,                         // multiplicity (dummy value)
	  0,                         // local_index (dummy value)
	  1.,                        // goodness_of_fit (dummy value)
	  0                          // dof
	  );
	if(wire_id.Plane==0)
	  hitCollection_U.emplace_back(hit.move(), wire, rawdigits);
	else if(wire_id.Plane==1)
	  hitCollection_V.emplace_back(hit.move(), wire, rawdigits);
	else if(wire_id.Plane==2)
	  hitCollection_Y.emplace_back(hit.move(), wire, rawdigits);

      }//End loop over time ticks on wire

      MF_LOG_DEBUG("TTHitFinder") << "Finished wire " << wire_id.Wire << " (plane " << wire_id.Plane << ")"
			       << "\tTotal hits (U,V,Y)= ("
			       << hitCollection_U.size() << ","
			       << hitCollection_V.size() << ","
			       << hitCollection_Y.size() << ")";

    }//End loop over all wires

    // put the hit collection and associations into the event
    mf::LogInfo("TTHitFinder") << "Total TTHitFinder hits (U,V,Y)=("
			       << hitCollection_U.size() << ","
			       << hitCollection_V.size() << ","
			       << hitCollection_Y.size() << ")";
    hitCollection_U.put_into(evt); // reminder: instance name was "uhits"
    hitCollection_V.put_into(evt); // reminder: instance name was "vhits"
    hitCollection_Y.put_into(evt); // reminder: instance name was "yhits"

  } // End of produce()


  //-------------------------------------------------
  float TTHitFinder::getTotalCharge(const float* signal_vector,int width=3,float threshold=-99){

    float totalCharge = 0;
    for(int tick=0; tick<width; tick++){
      if(signal_vector[tick] < threshold){
	totalCharge = -999; //special value for being below threshold
	break;
      }
      totalCharge += signal_vector[tick];
    }
    return totalCharge;

  }//end getTotalCharge method

  DEFINE_ART_MODULE(TTHitFinder)

} // end of hit namespace
