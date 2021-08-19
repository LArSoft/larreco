////////////////////////////////////////////////////////////////////
//HIT FINDER THAT RUNS ON RAW SIGNALS INSTEAD OF DECONVOLUTED ONES.
//WRITTEN INITIALLY FOR DUNE 35T ONLINE FILTER.
//abooth@fnal.gov, mstancar@fnal.gov.
////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>

//Framework
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"

//LArSoft
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

//LArSoft From FFT

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"

namespace hit {

  class RawHitFinder : public art::EDProducer {

    public:

      explicit RawHitFinder(fhicl::ParameterSet const& pset);

    private:
      void produce(art::Event& evt) override;


    unsigned int  fDataSize;                  //SIZE OF RAW DATA ON ONE WIRE.
    art::InputTag fDigitModuleLabel;          //MODULE THAT MADE DIGITS.
    std::string   fSpillName;                 //NOMINAL SPILL IS AN EMPTY STRING.

    //FFT COPIED VARIABLES.
    std::string         fCalDataModuleLabel;
    std::string         fHitLabelName;
    double              fMinSigInd;           //INDUCTION SIGNAL HEIGHT THRESHOLD.
    double              fMinSigCol;           //COLLECTION SIGNAL HEIGHT THRESHOLD.
    double              fIndWidth;            //INITIAL WIDTH FOR COLLECTION FIT.
    double              fColWidth;            //INITIAL WIDTH FOR COLLECTION FIT.
    double              fIndMinWidth;         //MINIMUM INDUCTION HIT WIDTH.
    double              fColMinWidth;         //MINIMUM COLLECTION HIT WIDTH.
    int                 fMaxMultiHit;         //MAXIMUM HITS FOR MULTIFIT.
    int                 fAreaMethod;          //TYPE OF AREA CALCULATION.
    std::vector<double> fAreaNorms;           //FACTORS FOR CONVERTING AREA TO SAME UNIT AS PEAK HEIGHT.
    bool                fUncompressWithPed;   //OPTION TO UNCOMPRESS WITH PEDESTAL.
    bool fSkipInd;  //OPTION TO SKIP INDUCTION PLANES
    double              fIncludeMoreTail;     //FRACTION OF THE HIT WIDTH INCLUDED BELOW THRESHOLD IN
                                              //    THE CALCULATION OF CHARGE.  COLLECTION PLANE ONLY
    int              fColMinWindow;       // Minimum length of integration window for charge in ticks
    int              fIndCutoff;          //MAX WIDTH FOR EARLY SIDE OF INDUCTION HIT IN TICKS

  }; // class RawHitFinder


  //-------------------------------------------------
  RawHitFinder::RawHitFinder(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
  {
    fDigitModuleLabel   = pset.get< art::InputTag >("DigitModuleLabel", "daq");
    fCalDataModuleLabel = pset.get< std::string  >("CalDataModuleLabel");
    fMinSigInd          = pset.get< double       >("MinSigInd");
    fMinSigCol          = pset.get< double       >("MinSigCol");
    fIncludeMoreTail    = pset.get< double       >("IncludeMoreTail", 0.);
    fIndWidth           = pset.get< int          >("IndWidth",20);
    fColWidth           = pset.get< double       >("ColWidth");
    fIndMinWidth        = pset.get< double       >("IndMinWidth");
    fColMinWidth        = pset.get< double       >("ColMinWidth",0.);
    fMaxMultiHit        = pset.get< int          >("MaxMultiHit");
    fAreaMethod         = pset.get< int          >("AreaMethod");
    fAreaNorms          = pset.get< std::vector< double > >("AreaNorms");
    fUncompressWithPed  = pset.get< bool         >("UncompressWithPed", true);
    fSkipInd = pset.get< bool         >("SkipInd", false);
    fColMinWindow        = pset.get< int       >("ColMinWindow",0);
    fIndCutoff        = pset.get< int       >("IndCutoff",20);
    mf::LogInfo("RawHitFinder_module") << "fDigitModuleLabel: " << fDigitModuleLabel << std::endl;

    //LET HITCOLLECTIONCREATOR DECLARE THAT WE ARE GOING TO PRODUCE
    //HITS AND ASSOCIATIONS TO RAW DIGITS BUT NOT ASSOCIATIONS TO WIRES
    //(WITH NO PARTICULAR PRODUCT LABEL).
    recob::HitCollectionCreator::declare_products(producesCollector(),
        /*instance_name*/"",
        /*doWireAssns*/false,
        /*doRawDigitAssns*/true);
  }

  //-------------------------------------------------
  void RawHitFinder::produce(art::Event& evt)
  {
    //GET THE GEOMETRY.
    art::ServiceHandle<geo::Geometry const> geom;

    //READ IN THE DIGIT LIST OBJECTS(S).
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;

    bool retVal = evt.getByLabel(fDigitModuleLabel, digitVecHandle);
    if(retVal==true)
      mf::LogInfo("RawHitFinder_module")    << "I got fDigitModuleLabel: "         << fDigitModuleLabel << std::endl;
    else
      mf::LogWarning("RawHitFinder_module") << "Could not get fDigitModuleLabel: " << fDigitModuleLabel << std::endl;

    std::vector<float> holder;      //HOLDS SIGNAL DATA.
    std::vector<short> rawadc;      //UNCOMPRESSED ADC VALUES.

    // ###############################################
    // ### Making a ptr vector to put on the event ###
    // ###############################################
    // THIS CONTAINS THE HIT COLLECTION AND ITS ASSOCIATIONS TO WIRES AND RAW DIGITS.
    recob::HitCollectionCreator hcol(evt, false /* doWireAssns */, true /* doRawDigitAssns */);

    std::vector<float> startTimes;  //STORES TIME OF WINDOW START.
    std::vector<float> maxTimes;    //STORES TIME OF LOCAL MAXIMUM.
    std::vector<float> endTimes;    //STORES TIME OF WINDOW END.
    std::vector<float> peakHeight;  //STORES ADC COUNT AT THE MAXIMUM.
    std::vector<float> hitrms;    	 //STORES CHARGE WEIGHTED RMS OF TIME ACROSS THE HIT.
    std::vector<double> charge;     //STORES THE TOTAL CHARGE ASSOCIATED WITH THE HIT.

    uint32_t channel   = 0;         //CHANNEL NUMBER.
    double   threshold = 0;         //MINIMUM SIGNAL SIZE FOR ID'ING A HIT.
    double   totSig    = 0;
    double   myrms     = 0;
    double   mynorm    = 0;
    geo::SigType_t sigType = geo::kInduction;
    std::stringstream numConv;

    hcol.reserve(digitVecHandle->size());
    auto ts = evt.time().value();
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){
      holder.clear();

      //GET THE REFERENCE TO THE CURRENT raw::RawDigit.
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel   = digitVec->Channel();
      fDataSize = digitVec->Samples();

      rawadc.resize(fDataSize);
      holder.resize(fDataSize);

      //UNCOMPRESS THE DATA.
      if (fUncompressWithPed){
        int pedestal = (int)digitVec->GetPedestal();
        raw::Uncompress(digitVec->ADCs(), rawadc, pedestal, digitVec->Compression());
      }
      else{
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
      }

      //GET THE LIST OF BAD CHANNELS.
      lariov::ChannelStatusProvider const& channelStatus
        = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

      lariov::ChannelSet_t const BadChannels
        = channelStatus.BadChannels(ts);

      for(unsigned int bin = 0; bin < fDataSize; ++bin){
        holder[bin]=(rawadc[bin]-digitVec->GetPedestal());
      }

      sigType = geom->SignalType(channel);

      peakHeight.clear();
      endTimes.clear();
      startTimes.clear();
      maxTimes.clear();
      charge.clear();
      hitrms.clear();

      bool channelSwitch = false;

      for(auto it = BadChannels.begin(); it != BadChannels.end(); it++)
      {
        if(channel==*it)
        {
          channelSwitch = true;
          break;
        }
      }

      if(channelSwitch==false)
      {
        // ###############################################
        // ###             Induction Planes            ###
        // ###############################################

        //THE INDUCTION PLANE METHOD HAS NOT YET BEEN MODIFIED AND TESTED FOR REAL DATA.
	// Or for detectors without a grid plane
	//
        if(sigType == geo::kInduction && !fSkipInd){
          threshold = fMinSigInd;
          //	std::cout<< "Threshold is " << threshold << std::endl;
          // fitWidth = fIndWidth;
          // minWidth = fIndMinWidth;
          //	continue;
          float negthr=-1.0*threshold;
          unsigned int bin =1;
          float minadc=0;

          // find the dips
          while (bin<(fDataSize-1)) {  // loop over ticks
            float thisadc = holder[bin]; float nextadc = holder[bin+1];
            if (thisadc<negthr && nextadc < negthr) { // new region, require two ticks above threshold
	      //              	    std::cout << "new region" << bin << " " << thisadc << std::endl;
              // step back to find zero crossing
		unsigned int place = bin;
		while (thisadc<=0 && bin>0) {
		  //		std::cout << bin << " " << thisadc << std::endl;
		  bin--;
		  thisadc=holder[bin];
		}
		float hittime = bin+thisadc/(thisadc-holder[bin+1]);
		maxTimes.push_back(hittime);

              // step back more to find the hit start time
	      uint32_t stop;
	      if (fIndCutoff<(int)bin) {stop=bin-fIndCutoff;} else {stop=0;}
	      while (thisadc<threshold && bin>stop) {
		//		std::cout << bin << " " << thisadc << std::endl;
                bin--;
                thisadc=holder[bin];
              }
              if (bin>=2) bin-=2;
	      while (thisadc>threshold && bin>stop) {
		//		std::cout << bin << " " << thisadc << std::endl;
                bin--;
                thisadc=holder[bin];
              }
              startTimes.push_back(bin+1);
              // now step forward from hit time to find end time, area of dip
              bin=place;
              thisadc=holder[bin];
              minadc=thisadc;
	      bin++;
              totSig = fabs(thisadc);
              while (thisadc<negthr && bin<fDataSize) {
                totSig += fabs(thisadc);
                thisadc=holder[bin];
                if (thisadc<minadc) minadc=thisadc;
                bin++;
              }
              endTimes.push_back(bin-1);
              peakHeight.push_back(-1.0*minadc);
              charge.push_back(totSig);
              hitrms.push_back(5.0);
              //	    std::cout << "TOTAL SIGNAL INDUCTION " << totSig << "  5.0" << std::endl;
              // std::cout << "filled end times " << bin-1 << "peak height vector size " << peakHeight.size() << std::endl;

              // don't look for a new hit until it returns to baseline
              while (thisadc<0 && bin<fDataSize) {
                //	      std::cout << bin << " " << thisadc << std::endl;
                bin++;
                if (bin == fDataSize) break;
                thisadc=holder[bin];
              }
            } // end region
            bin++;
          }// loop over ticks
        }

        // ###############################################
        // ###             Collection Plane            ###
        // ###############################################

        else if(sigType == geo::kCollection)
        {
          threshold = fMinSigCol;

          float madc = threshold;
          int ibin   = 0;
          int start  = 0;
          int end    = 0;
          unsigned int bin = 0;

          while (bin<fDataSize)
          {
            float thisadc = holder[bin];
            madc = threshold;
            ibin = 0;

            if (thisadc>madc)
            {
              start = bin;

              if(thisadc>threshold && bin<fDataSize)
              {
                while (thisadc>threshold && bin<fDataSize)
                {
                  if (thisadc>madc)
                  {
                    ibin=bin;
                    madc=thisadc;
                  }
                  bin++;
                  if (bin == fDataSize) break;
                  thisadc=holder[bin];
                }
              }
              else
              {
                bin++;
              }

              end = bin-1;

              if(start!=end)
              {
                maxTimes.push_back(ibin);
                peakHeight.push_back(madc);
                startTimes.push_back(start);
                endTimes.push_back(end);

                totSig = 0;
                myrms  = 0;
                mynorm = 0;

                int moreTail = std::ceil(fIncludeMoreTail*(end-start));
		if (moreTail<fColMinWindow) moreTail=fColMinWindow;

                for(int i = start-moreTail; i <= end+moreTail; i++)
                {
                  if(i<(int)(holder.size()) && i>=0)
                  {
                    float temp = ibin-i;
                    myrms += temp*temp*holder[i];

                    totSig += holder[i];
                  }
                }

                charge.push_back(totSig);
                mynorm = totSig;
                myrms/=mynorm;
                hitrms.push_back(sqrt(myrms));

                //PRE CHANGES MADE 04/14/16. A BOOTH, DUNE 35T.
                /*
                   int moreTail = std::ceil(fIncludeMoreTail*(end-start));

                   for(int i = start-moreTail; i <= end+moreTail; i++)
                   {
                   totSig += holder[i];
                   float temp2 = holder[i]*holder[i];
                   mynorm += temp2;
                   float temp = ibin-i;
                   myrms += temp*temp*temp2;
                   }

                   charge.push_back(totSig);
                   myrms/=mynorm;
                   if((end-start+2*moreTail+1)!=0)
                   {
                   myrms/=(float)(end-start+2*moreTail+1);
                   hitrms.push_back(sqrt(myrms));
                   }
                   else
                   {
                   hitrms.push_back(sqrt(myrms));
                   }*/
              }
            }
            start = 0;
            end = 0;
            bin++;
          }
        }
      }

      int    numHits(0);                       //NUMBER OF CONSECUTIVE HITS BEING FITTED.
      int    hitIndex(0);                      //INDEX OF CURRENT HIT IN SEQUENCE.
      double amplitude(0), position(0);        //FIT PARAMETERS.
      double start(0), end(0);
      double amplitudeErr(0), positionErr(0);  //FIT ERRORS.
      double goodnessOfFit(0), chargeErr(0);   //CHI2/NDF and error on charge.
      double hrms(0);

      numHits = maxTimes.size();
      for (int i = 0; i < numHits; ++i)
      {
        amplitude     = peakHeight[i];
        position      = maxTimes[i];
        start         = startTimes[i];
        end           = endTimes[i];
        hrms          = hitrms[i];
        amplitudeErr  = -1;
        positionErr   = 1.0;
        goodnessOfFit = -1;
        chargeErr     = -1;
        totSig        = charge[i];


        std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        geo::WireID wid = wids[0];

        if (start>=end)
        {
          mf::LogWarning("RawHitFinder_module") << "Hit start " << start << " is >= hit end " << end;
          continue;
        }

        recob::HitCreator hit(
            *digitVec,                                                                     //RAW DIGIT REFERENCE.
            wid,                                                                           //WIRE ID.
            start,                                                                         //START TICK.
            end,                                                                           //END TICK.
            hrms,                                                                          //RMS.
            position,                                                                      //PEAK_TIME.
            positionErr,                                                                   //SIGMA_PEAK_TIME.
            amplitude,                                                                     //PEAK_AMPLITUDE.
            amplitudeErr,                                                                  //SIGMA_PEAK_AMPLITUDE.
            totSig,                                                                        //HIT_INTEGRAL.
            chargeErr,                                                                     //HIT_SIGMA_INTEGRAL.
            std::accumulate(holder.begin() + (int) start, holder.begin() + (int) end, 0.), //SUMMED CHARGE.
            1,                                                                             //MULTIPLICITY.
            -1,                                                                            //LOCAL_INDEX.
            goodnessOfFit,                                                                 //WIRE ID.
            int(end-start+1)                                                               //DEGREES OF FREEDOM.
            );
        hcol.emplace_back(hit.move(), digitVec);

        ++hitIndex;
      }
    }

    hcol.put_into(evt);
  }

  DEFINE_ART_MODULE(RawHitFinder)

} // end namespace hit
