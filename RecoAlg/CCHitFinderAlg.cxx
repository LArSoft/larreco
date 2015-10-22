//////////////////////////////////////////////////////////////////////
///
/// CCHitFinder class
///
/// Bruce Baller, baller@fnal.gov
///
/// Find hits for ClusterCrawler and put them in a temporary struct.
/// These hits may be modified by ClusterCrawler before saving them
/// in the event
///
////////////////////////////////////////////////////////////////////////


// class header
#include "RecoAlg/CCHitFinderAlg.h"

// C/C++ standard libraries
#include <cmath> // std::sqrt(), std::abs()
#include <iostream>
#include <iomanip>
#include <sstream>
#include <array>
#include <vector>
#include <utility> // std::pair<>, std::make_pair()
#include <algorithm> // std::sort(), std::copy()

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft Includes
#include "SimpleTypesAndConstants/RawTypes.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/SimpleFits.h" // lar::util::GaussianFit<>
#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"

// ROOT Includes
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"


namespace hit {
  
  constexpr unsigned int CCHitFinderAlg::MaxGaussians; // definition
  
//------------------------------------------------------------------------------
  CCHitFinderAlg::CCHitFinderAlg(fhicl::ParameterSet const& pset):
    FitCache(new
      GausFitCache // run-time, on demand TFormula cache
    //  CompiledGausFitCache<MaxGaussians> // precompiled Gaussian set
    //  CompiledTruncatedGausFitCache<MaxGaussians, 4> // precompiled truncated Gaussian set
    //  CompiledTruncatedGausFitCache<MaxGaussians, 5> // precompiled truncated Gaussian set
      ("GausFitCache_CCHitFinderAlg")
      )
  {
    this->reconfigure(pset);
  }

  void CCHitFinderAlg::reconfigure(fhicl::ParameterSet const& pset)
  {
    fMinPeak            = pset.get<std::vector<float>>("MinPeak");
    fMinRMS             = pset.get<std::vector<float>>("MinRMS");
    fMaxBumps           = pset.get<unsigned short>("MaxBumps");
    fMaxXtraHits        = pset.get<unsigned short>("MaxXtraHits");
    fChiSplit           = pset.get<float>("ChiSplit");
    fChiNorms           = pset.get<std::vector< float > >("ChiNorms");
    fUseFastFit         = pset.get<bool>("UseFastFit", false);
    fUseChannelFilter   = pset.get<bool>("UseChannelFilter", true);
    fStudyHits          = pset.get<bool>("StudyHits", false);
    // The following variables are only used in StudyHits mode
    fUWireRange         = pset.get< std::vector< short >>("UWireRange");
    fUTickRange         = pset.get< std::vector< short >>("UTickRange");
    fVWireRange         = pset.get< std::vector< short >>("VWireRange");
    fVTickRange         = pset.get< std::vector< short >>("VTickRange");
    fWWireRange         = pset.get< std::vector< short >>("WWireRange");
    fWTickRange         = pset.get< std::vector< short >>("WTickRange");
    
    if(fMinPeak.size() != fMinRMS.size()) {
      mf::LogError("CCTF")<<"MinPeak size != MinRMS size";
      return;
    }
    
    if (fMaxBumps > MaxGaussians) {
      // LOG_WARNING will point the user to this line of code.
      // Any value of MaxGaussians can be used.
      // That value is defined in the header file.
      LOG_WARNING("CCHitFinderAlg")
        << "CCHitFinder algorithm is currently hard-coded to support at most "
        << MaxGaussians << " bumps per region of interest, but " << fMaxBumps
        << " have been requested.\n"
        << "We are forcing the parameter to " << MaxGaussians
        << ". If this is not acceptable, increase CCHitFinderAlg::MaxGaussians"
        << " value and recompile.";
      fMaxBumps = MaxGaussians;
    } // if too many gaussians
    
    FinalFitStats.Reset(MaxGaussians);
    TriedFitStats.Reset(MaxGaussians);
    
    // sanity check for StudyHits mode
    if(fStudyHits) {
      if(fUWireRange.size() != 2 || fUTickRange.size() != 2 ||
         fVWireRange.size() != 2 || fVTickRange.size() != 2 ||
         fWWireRange.size() != 2 || fWTickRange.size() != 2) {
        mf::LogError("CCHF")<<"Invalid vector size for StudyHits. Must be 2";
        return;
      }
    } // fStudyHits

  }

//------------------------------------------------------------------------------
  CCHitFinderAlg::HitChannelInfo_t::HitChannelInfo_t
    (recob::Wire const* w, geo::WireID wid, geo::Geometry const& geom):
    wire(w),
    wireID(wid),
    sigType(geom.SignalType(w->Channel()))
    {}
  
//------------------------------------------------------------------------------
  void CCHitFinderAlg::RunCCHitFinder(std::vector<recob::Wire> const& Wires) {
  
    allhits.clear();

    unsigned short maxticks = 1000;
    float *ticks = new float[maxticks];
    // define the ticks array used for fitting 
    for(unsigned short ii = 0; ii < maxticks; ++ii) {
      ticks[ii] = ii;
    }
    float *signl = new float[maxticks];
    float adcsum = 0;
    // initialize the vectors for the hit study
    if(fStudyHits) StudyHits(0);
    bool first;

//    prt = false;
   
    lariov::IChannelStatusProvider const& channelStatus
      = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();

    for(size_t wireIter = 0; wireIter < Wires.size(); wireIter++){

      recob::Wire const& theWire = Wires[wireIter];
      theChannel = theWire.Channel();
      // ignore bad channels
      if(channelStatus.IsBad(theChannel)) continue;
/*
      geo::SigType_t SigType = geom->SignalType(theChannel);
      minSig = 0.;
      minRMS = 0.;
      if(SigType == geo::kInduction){
        minSig = fMinSigInd;
        minRMS = fMinRMSInd;
      }//<-- End if Induction Plane
      else if(SigType == geo::kCollection){
        minSig = fMinSigCol;
        minRMS  = fMinRMSCol;
      }//<-- End if Collection Plane
*/

      std::vector<geo::WireID> wids = geom->ChannelToWire(theChannel);
      thePlane = wids[0].Plane;
      if(thePlane > fMinPeak.size() - 1) {
        mf::LogError("CCHF")<<"MinPeak vector too small for plane "<<thePlane;
        return;
      }
      theWireNum = wids[0].Wire;
      HitChannelInfo_t WireInfo(&theWire, wids[0], *geom);
      
      // minimum number of time samples
      unsigned short minSamples = 2 * fMinRMS[thePlane];

      // factor used to normalize the chi/dof fits for each plane
      chinorm = fChiNorms[thePlane];

      // edit this line to debug hit fitting on a particular plane/wire
//      prt = (thePlane == 1 && theWireNum == 839);
      std::vector<float> signal(theWire.Signal());
      
      unsigned short nabove = 0;
      unsigned short tstart = 0;
      unsigned short maxtime = signal.size() - 2;
      // find the min time when the signal is below threshold
      unsigned short mintime = 3;
      for(unsigned short time = 3; time < maxtime; ++time) {
        if(signal[time] < fMinPeak[thePlane]) {
          mintime = time;
          break;
        }
      }
      for(unsigned short time = mintime; time < maxtime; ++time) {
        if(signal[time] > fMinPeak[thePlane]) {
          if(nabove == 0) tstart = time;
          ++nabove;
        } else {
          // check for a wide enough signal above threshold
          if(nabove > minSamples) {
            // skip this wire if the RAT is too long
            if(nabove > maxticks) mf::LogError("CCHitFinder")
              <<"Long RAT "<<nabove<<" "<<maxticks
              <<" No signal on wire "<<theWireNum<<" after time "<<time;
            if(nabove > maxticks) break;
            unsigned short npt = 0;
            // look for bumps to inform the fit
            bumps.clear();
            adcsum = 0;
            for(unsigned short ii = tstart; ii < time; ++ii) {
              signl[npt] = signal[ii];
              adcsum += signl[npt];
              if(signal[ii    ] > signal[ii - 1] &&
                 signal[ii - 1] > signal[ii - 2] &&
                 signal[ii    ] > signal[ii + 1] &&
                 signal[ii + 1] > signal[ii + 2]) bumps.push_back(npt);
//  if(prt) mf::LogVerbatim("CCHitFinder")<<"signl "<<ii<<" "<<signl[npt];
              ++npt;
            }
            // decide if this RAT should be studied
            if(fStudyHits) StudyHits(1, npt, ticks, signl, tstart);
            // just make a crude hit if too many bumps
            if(bumps.size() > fMaxBumps) {
              MakeCrudeHit(npt, ticks, signl);
              StoreHits(tstart, npt, WireInfo, adcsum);
              nabove = 0;
              continue;
            }
            // start looking for hits with the found bumps
            unsigned short nHitsFit = bumps.size();
            unsigned short nfit = 0;
            chidof = 0.;
            dof = -1;
            bool HitStored = false;
            unsigned short nMaxFit = bumps.size() + fMaxXtraHits;
            // only used in StudyHits mode
            first = true;
            while(nHitsFit <= nMaxFit) {
      
              FitNG(nHitsFit, npt, ticks, signl);
              if(fStudyHits && first && SelRAT) {
                first = false;
                StudyHits(2, npt, ticks, signl, tstart);
              }
              // good chisq so store it
              if(chidof < fChiSplit) {
                StoreHits(tstart, npt, WireInfo, adcsum);
                HitStored = true;
                break;
              }
              // the previous fit was better, so revert to it and
              // store it
              ++nHitsFit;
              ++nfit;
            } // nHitsFit < fMaxXtraHits
            if( !HitStored && npt < maxticks) {
              // failed all fitting. Make a crude hit
              MakeCrudeHit(npt, ticks, signl);
              StoreHits(tstart, npt, WireInfo, adcsum);
            }
            else if (nHitsFit > 0) FinalFitStats.AddMultiGaus(nHitsFit);
          } // nabove > minSamples
          nabove = 0;
        } // signal < fMinPeak
      } // time
    } // wireIter

    // print out
    if(fStudyHits) StudyHits(4);

    delete[] ticks;
    delete[] signl;

  } //RunCCHitFinder


/////////////////////////////////////////
  bool CCHitFinderAlg::FastGaussianFit(
    unsigned short npt, float const*ticks, float const*signl,
    std::array<double, 3>& params,
    std::array<double, 3>& paramerrors,
    float& chidof
  ) {
    // parameters: amplitude, mean, sigma
    
    lar::util::GaussianFit<double> fitter; // probably "double" is overdoing
    
    // apply a time shift so that the center of the time interval is 0
    const float time_shift = (ticks[npt - 1] + ticks[0]) / 2.F;
    
    // fill the input data, no uncertainty
    for (size_t i = 0; i < npt; ++i) {
      if (signl[i] <= 0) {
        LOG_DEBUG("CCHitFinderAlg")
          << "Non-positive charge encountered. Backing up to ROOT fit.";
        return false;
      }
      // we could freely add a Poisson uncertainty (as third parameter)
      fitter.add(ticks[i] - time_shift, signl[i]);
    } // for
    
    // we might have found that we don't want the fast fit after all...
    if (!fitter.FillResults(params, paramerrors)) {
      // something went wrong...
      LOG_DEBUG("CCHitFinderAlg") << "Fast Gaussian fit failed.";
      return false;
    }
    
    // note that this is not the full chi^2, but it is the chi^2 of the
    // parabolic fit underlying the Gaussian one
    const double chi2 = fitter.ChiSquare();
    chidof = chi2 / fitter.NDF();
    
    // remove the time shift
    params[1] += time_shift; // mean
    
    // GP: inflate the uncertainties on the fit parameters according to chi2/NDF
    // (not sure if this is in any way justified)
    if (chidof > 1.)
      for (double& par: paramerrors) par *= std::sqrt(chidof);
    
    return true;
  } // FastGaussianFit()
  
  
/////////////////////////////////////////
  void CCHitFinderAlg::FitNG(unsigned short nGaus, unsigned short npt, 
    float *ticks, float *signl)
  {
    // Fit the signal to n Gaussians

    dof = npt - 3 * nGaus;
    
    chidof = 9999.;

    if(dof < 3) return;
    if(bumps.size() == 0) return;
    
    // load the fit into a temp vector
    std::vector<double> partmp;
    std::vector<double> partmperr;
    
    //
    // if it is possible, we try first with the quick single Gaussian fit
    //
    TriedFitStats.AddMultiGaus(nGaus);
    
    bool bNeedROOTfit = (nGaus > 1) || !fUseFastFit;
    if (!bNeedROOTfit) {
      // so, we need only one puny Gaussian;
      std::array<double, 3> params, paramerrors;
      
      TriedFitStats.AddFast();
      
      if (FastGaussianFit(npt, ticks, signl, params, paramerrors, chidof)) {
        // success? copy the results in the proper structures
        partmp.resize(3);
        std::copy(params.begin(), params.end(), partmp.begin());
        partmperr.resize(3);
        std::copy(paramerrors.begin(), paramerrors.end(), partmperr.begin());
      }
      else bNeedROOTfit = true; // if we fail, let's schedule ROOT to back us up
      
      if (!bNeedROOTfit) FinalFitStats.AddFast();
      
    } // if we don't need ROOT to fit
    
    if (bNeedROOTfit) {
      // we may land here either because the simple Gaussian fit did not work
      // (either failed, or we chose not to trust it)
      // or because the fit is multi-Gaussian
      
      // define the fit string to pass to TF1
      
      std::stringstream numConv;
      std::string eqn = "gaus";
      if(nGaus > 1) eqn = "gaus(0)";
      for(unsigned short ii = 3; ii < nGaus*3; ii+=3){
        eqn.append(" + gaus(");
        numConv.str("");
        numConv << ii;
        eqn.append(numConv.str());
        eqn.append(")");
      }
      
      std::unique_ptr<TF1> Gn(new TF1("gn",eqn.c_str()));
      /*
      TF1* Gn = FitCache->Get(nGaus);
      */
      TGraph *fitn = new TGraph(npt, ticks, signl);
  /*
    if(prt) mf::LogVerbatim("CCHitFinder")
      <<"FitNG nGaus "<<nGaus<<" nBumps "<<bumps.size();
  */
      // put in the bump parameters. Assume that nGaus >= bumps.size()
      for(unsigned short ii = 0; ii < bumps.size(); ++ii) {
        unsigned short index = ii * 3;
        unsigned short bumptime = bumps[ii];
        double amp = signl[bumptime];
        Gn->SetParameter(index    , amp);
        Gn->SetParLimits(index, 0., 9999.);
        Gn->SetParameter(index + 1, (double)bumptime);
        Gn->SetParLimits(index + 1, 0, (double)npt);
        Gn->SetParameter(index + 2, (double)fMinRMS[thePlane]);
        Gn->SetParLimits(index + 2, 1., 3*(double)fMinRMS[thePlane]);
  /*
    if(prt) mf::LogVerbatim("CCHitFinder")<<"Bump params "<<ii<<" "<<(short)amp
      <<" "<<(int)bumptime<<" "<<(int)fMinRMS[thePlane];
  */
      } // ii bumps
  
      // search for other bumps that may be hidden by the already found ones
      for(unsigned short ii = bumps.size(); ii < nGaus; ++ii) {
        // bump height must exceed fMinPeak
        float big = fMinPeak[thePlane];
        unsigned short imbig = 0;
        for(unsigned short jj = 0; jj < npt; ++jj) {
          float diff = signl[jj] - Gn->Eval((Double_t)jj, 0, 0, 0);
          if(diff > big) {
            big = diff;
            imbig = jj;
          }
        } // jj
        if(imbig > 0) {
  /*
    if(prt) mf::LogVerbatim("CCHitFinder")<<"Found bump "<<ii<<" "<<(short)big
      <<" "<<imbig;
  */
          // set the parameters for the bump
          unsigned short index = ii * 3;
          Gn->SetParameter(index    , (double)big);
          Gn->SetParLimits(index, 0., 9999.);
          Gn->SetParameter(index + 1, (double)imbig);
          Gn->SetParLimits(index + 1, 0, (double)npt);
          Gn->SetParameter(index + 2, (double)fMinRMS[thePlane]);
          Gn->SetParLimits(index + 2, 1., 5*(double)fMinRMS[thePlane]);
        } // imbig > 0
      } // ii 
      
      // W = set weights to 1, N = no drawing or storing, Q = quiet
      // B = bounded parameters
      fitn->Fit(&*Gn,"WNQB");
      
      for(unsigned short ipar = 0; ipar < 3 * nGaus; ++ipar) {
        partmp.push_back(Gn->GetParameter(ipar));
        partmperr.push_back(Gn->GetParError(ipar));
      }
      chidof = Gn->GetChisquare() / ( dof * chinorm);
      
      delete fitn;
    //  delete Gn;
      
    } // if ROOT fit
    
    // Sort by increasing time if necessary
    if(nGaus > 1) {
      std::vector< std::pair<unsigned short, unsigned short> > times;
      // fill the sort vector
      for(unsigned short ii = 0; ii < nGaus; ++ii) {
        unsigned short index = ii * 3;
        times.push_back(std::make_pair(partmp[index + 1],ii));
      } // ii
      std::sort(times.begin(), times.end());
      // see if re-arranging is necessary
      bool sortem = false;
      for(unsigned short ii = 0; ii < nGaus; ++ii) {
        if(times[ii].second != ii) {
          sortem = true;
          break;
        }
      } // ii
      if(sortem) {
        // temp temp vectors for putting things in the right time order
        std::vector<double> partmpt;
        std::vector<double> partmperrt;
        for(unsigned short ii = 0; ii < nGaus; ++ii) {
          unsigned short index = times[ii].second * 3;
          partmpt.push_back(partmp[index]);
          partmpt.push_back(partmp[index+1]);
          partmpt.push_back(partmp[index+2]);
          partmperrt.push_back(partmperr[index]);
          partmperrt.push_back(partmperr[index+1]);
          partmperrt.push_back(partmperr[index+2]);
        } // ii
        partmp = partmpt;
        partmperr = partmperrt;
      } // sortem
    } // nGaus > 1
/*
  if(prt) {
    mf::LogVerbatim("CCHitFinder")<<"Fit "<<nGaus<<" chi "<<chidof
      <<" npars "<<partmp.size();
    mf::LogVerbatim("CCHitFinder")<<"pars    errs ";
    for(unsigned short ii = 0; ii < partmp.size(); ++ii) {
      mf::LogVerbatim("CCHitFinder")<<ii<<" "<<partmp[ii]<<" "
        <<partmperr[ii];
    }
  }
*/
    // ensure that the fit is reasonable
    bool fitok = true;
    for(unsigned short ii = 0; ii < nGaus; ++ii) {
      unsigned short index = ii * 3;
      // ensure that the fitted time is within the signal bounds
      short fittime = partmp[index + 1];
      if(fittime < 0 || fittime > npt - 1) {
        fitok = false;
        break;
      }
      // ensure that the signal peak is large enough
      if(partmp[index] < fMinPeak[thePlane]) {
        fitok = false;
        break;
      }
      // ensure that the RMS is large enough but not too large
      float rms = partmp[index + 2];
      if(rms < 0.5 * fMinRMS[thePlane] || rms > 5 * fMinRMS[thePlane]) {
        fitok = false;
        break;
      }
      // ensure that the hits are not too similar in time (< 2 ticks)
      for(unsigned short jj = 0; jj < nGaus; ++jj) {
        if(jj == ii) continue;
        unsigned short jndex = jj * 3;
        float timediff = std::abs(partmp[jndex + 1] - partmp[index + 1]);
        if(timediff < 2.) {
          fitok = false;
          break;
        }
      }
      if(!fitok) break;
    }

    if(fitok) {
      par = partmp;
      parerr = partmperr;
    } else {
      chidof = 9999.;
      dof = -1;
//      if(prt) mf::LogVerbatim("CCHitFinder")<<"Bad fit parameters";
    }
    
    return;
  } // FitNG

/////////////////////////////////////////
  void CCHitFinderAlg::MakeCrudeHit(unsigned short npt, 
    float *ticks, float *signl)
  {
    // make a single crude hit if fitting failed
    float sumS = 0.;
    float sumST = 0.;
    for(unsigned short ii = 0; ii < npt; ++ii) {
      sumS  += signl[ii];
      sumST += signl[ii] * ticks[ii];
    }
    float mean = sumST / sumS;
    float rms = 0.;
    for(unsigned short ii = 0; ii < npt; ++ii) {
      float arg = ticks[ii] - mean;
      rms += signl[ii] * arg * arg;
    }
    rms = std::sqrt(rms / sumS);
    float amp = sumS / (Sqrt2Pi * rms);
    par.clear();
/*
  if(prt) mf::LogVerbatim("CCHitFinder")<<"Crude hit Amp "<<(int)amp<<" mean "
    <<(int)mean<<" rms "<<rms;
*/
    par.push_back(amp);
    par.push_back(mean);
    par.push_back(rms);
    // need to do the errors better
    parerr.clear();
    float amperr = npt;
    float meanerr = std::sqrt(1/sumS);
    float rmserr = 0.2 * rms;
    parerr.push_back(amperr);
    parerr.push_back(meanerr);
    parerr.push_back(rmserr);
/*
  if(prt) mf::LogVerbatim("CCHitFinder")<<" errors Amp "<<amperr<<" mean "
    <<meanerr<<" rms "<<rmserr;
*/
    chidof = 9999.;
    dof = -1;
  } // MakeCrudeHit


/////////////////////////////////////////
  void CCHitFinderAlg::StoreHits(unsigned short TStart, unsigned short npt,
    HitChannelInfo_t info, float adcsum
  ) {
    // store the hits in the struct
    size_t nhits = par.size() / 3;
    
    if(allhits.max_size() - allhits.size() < nhits) {
      mf::LogError("CCHitFinder")
        << "Too many hits: existing " << allhits.size() << " plus new " << nhits
        << " beyond the maximum " << allhits.max_size();
      return;
    }
    
    if(nhits == 0) return;

    // fill RMS for single hits
    if(fStudyHits) StudyHits(3);

    const float loTime = TStart;
    const float hiTime = TStart + npt;

    // Find sum of the areas of all Gaussians
    float gsum = 0.;
    for(size_t hit = 0; hit < nhits; ++hit) {
      const unsigned short index = 3 * hit;
      gsum += Sqrt2Pi * par[index] * par[index + 2];
    }
    for(size_t hit = 0; hit < nhits; ++hit) {
      const size_t index = 3 * hit;
      const float charge = Sqrt2Pi * par[index] * par[index + 2];
      const float charge_err = SqrtPi
        * (parerr[index] * par[index + 2] + par[index] * parerr[index + 2]);
      
      allhits.emplace_back(
        info.wire->Channel(),     // channel
        loTime,                   // start_tick
        hiTime,                   // end_tick
        par[index + 1] + TStart,  // peak_time
        parerr[index + 1],        // sigma_peak_time
        par[index + 2],           // rms
        par[index],               // peak_amplitude
        parerr[index],            // sigma_peak_amplitude
        adcsum * charge / gsum,   // summedADC
        charge,                   // hit_integral
        charge_err,               // hit_sigma_integral
        nhits,                    // multiplicity
        hit,                      // local_index
        chidof,                   // goodness_of_fit
        dof,                      // dof
        info.wire->View(),        // view
        info.sigType,             // signal_type
        info.wireID               // wireID
        );
/*
  if(prt) {
    mf::LogVerbatim("CCHitFinder")<<"W:T "<<allhits.back().WireID().Wire
      <<":"<<(short)allhits.back().PeakTime()
      <<" Chg "<<(short)allhits.back().Integral()
      <<" RMS "<<allhits.back().RMS()
      <<" lo ID "<<allhits.back().LocalIndex()
      <<" numHits "<<allhits.back().Multiplicity()
      <<" loTime "<<allhits.back().StartTick()<<" hiTime "<<allhits.back().EndTick()
      <<" chidof "<<allhits.back().GoodnessOfFit()
      << " DOF " << allhits.back().DegreesOfFreedom();
  }
*/
    } // hit
  } // StoreHits


//////////////////////////////////////////////////
  void CCHitFinderAlg::StudyHits(unsigned short flag, unsigned short npt,
      float *ticks, float *signl, unsigned short tstart) {
    // study hits in user-selected ranges of wires and ticks in each plane. The user should identify
    // a shallow-angle isolated track, e.g. using the event display, to determine the wire/tick ranges.
    // One hit should be reconstructed on each wire when the hit finding fcl parameters are set correctly.
    // The intent of this study is to determine the correct fcl parameters. The flag variable determines
    // the operation performed.
    // flag = 0: Initialize study vectors
    // flag = 1: Set SelRat true if the Region Above Threshold resides within a wire/hit range
    // flag = 2: Find the maximum signal and calculate the RMS. Also find the low and high ticks of signals
    //           in the wire range to allow a later calculation of the track angle. This isn't strictly 
    //           necessary for the study and presumes that the user has selected compatible regions in each plane.
    // flag = 3: Accumulate the RMS from the first Gaussian fit
    // flag = 4: Calculate recommended fcl parameters and print the results to the screen

    // init
    if(flag == 0) {
      for(unsigned short ipl = 0; ipl < 3; ++ipl) {
        // Average chisq of the first fit on a single bump in each plane
        bumpChi.push_back(0.);
        // Average RMS of the dump
        bumpRMS.push_back(0.);
        // The number of single bumps in each plane
        bumpCnt.push_back(0.);
        // number of RATs
        RATCnt.push_back(0);
        // The number of single hits found in each plane
        hitCnt.push_back(0.);
        // Average reconstructed hit RMS
        hitRMS.push_back(0.);
        // lo/hi wire/time
        loWire.push_back(9999.);
        loTime.push_back(0.);
        hiWire.push_back(-1.);
        hiTime.push_back(0.);
      } // ii
      return;
    } // flag == 0
    
    if(flag == 1) {
      SelRAT = false;
      if(thePlane == 0) {
        if(theWireNum > fUWireRange[0] && theWireNum < fUWireRange[1] && 
           tstart     > fUTickRange[0] && tstart     < fUTickRange[1]) {
          SelRAT = true;
          RATCnt[thePlane] += 1;
        }
        return;
      } // thePlane == 0
      if(thePlane == 1) {
        if(theWireNum > fVWireRange[0] && theWireNum < fVWireRange[1] && 
           tstart     > fVTickRange[0] && tstart     < fVTickRange[1]) {
          SelRAT = true;
          RATCnt[thePlane] += 1;
        }
        return;
      } // thePlane == 1
      if(thePlane == 2) {
        if(theWireNum > fWWireRange[0] && theWireNum < fWWireRange[1] && 
           tstart     > fWTickRange[0] && tstart     < fWTickRange[1]) {
          SelRAT = true;
          RATCnt[thePlane] += 1;
        }
        return;
      } // thePlane == 2
    } // flag == 1
    
    if(flag == 2) {
      if(!SelRAT) return;
      // in this section we find the low/hi wire/time for a signal. This can be used to calculate
      // the slope dT/dW to study hit width, fraction of crude hits, etc vs dT/dW
      float big = 0.;
      float imbig = 0.;
      for(unsigned short ii = 0; ii < npt; ++ii) {
        if(signl[ii] > big) {
          big = signl[ii];
          imbig = ii;
        }
      } // ii
      // require a significant PH 
      if(big > fMinPeak[0]) {
        // get the Lo info
        if(theWireNum < loWire[thePlane]) {
          loWire[thePlane] = theWireNum;
          loTime[thePlane] = tstart + imbig;
        }
        // get the Hi info
        if(theWireNum > hiWire[thePlane]) {
          hiWire[thePlane] = theWireNum;
          hiTime[thePlane] = tstart + imbig;
        }
      } // big > fMinPeak[0]
      if(bumps.size() == 1 && chidof < 9999.) {
        bumpCnt[thePlane] += bumps.size();
        bumpChi[thePlane] += chidof;
        // calculate the average bin
        float sumt = 0.;
        float sum = 0.;
        for(unsigned short ii = 0; ii < npt; ++ii) {
          sum  += signl[ii];
          sumt += signl[ii] * ii;
        } // ii
        float aveb = sumt / sum;
        // now calculate the RMS
        sumt = 0.;
        for(unsigned short ii = 0; ii < npt; ++ii) {
          float dbin = (float)ii - aveb;
          sumt += signl[ii] * dbin * dbin;
        } // ii
        bumpRMS[thePlane] += std::sqrt(sumt / sum);
      } // bumps.size() == 1 && chidof < 9999.
      return;
    } // flag == 2    

    // fill info for single hits
    if(flag == 3) {
      if(!SelRAT) return;
      if(par.size() == 3) {
        hitCnt[thePlane] += 1;
        hitRMS[thePlane] += par[2];
      }
      return;
    }


    if(flag == 4) {
      // The goal is to adjust the fcl inputs so that the number of single 
      // hits found is ~equal to the number of single bumps found for shallow
      // angle tracks. The ChiNorm inputs should be adjusted so the average
      //  chisq/DOF is ~1 in each plane.
      std::cout<<"Check lo and hi W/T for each plane"<<std::endl;
      for(unsigned short ipl = 0; ipl < 3; ++ipl) {
        std::cout<<ipl<<" lo "<<loWire[ipl]<<" "<<loTime[ipl]
          <<" hi "<<hiWire[ipl]<<" "<<hiTime[ipl]<<std::endl;
      }
      std::cout<<" ipl nRAT bCnt   bChi   bRMS hCnt   hRMS  dT/dW New_ChiNorm"<<std::endl;
      for(unsigned short ipl = 0; ipl < 3; ++ipl) {
        if(bumpCnt[ipl] > 0) {
          bumpChi[ipl] = bumpChi[ipl] / (float)bumpCnt[ipl];
          bumpRMS[ipl] = bumpRMS[ipl] / (float)bumpCnt[ipl];
          hitRMS[ipl]  = hitRMS[ipl]  / (float)hitCnt[ipl];
          // calculate the slope
          float dTdW = std::abs((hiTime[ipl] - loTime[ipl]) / (hiWire[ipl] - loWire[ipl]));
          std::cout<<ipl<<std::right<<std::setw(5)<<RATCnt[ipl]
            <<std::setw(5)<<bumpCnt[ipl]
            <<std::setw(7)<<std::fixed<<std::setprecision(2)<<bumpChi[ipl]
            <<std::setw(7)<<bumpRMS[ipl]
            <<std::setw(7)<<hitCnt[ipl]
            <<std::setw(7)<<std::setprecision(1)<<hitRMS[ipl]
            <<std::setw(7)<<dTdW
            <<std::setw(7)<<std::setprecision(2)
            <<bumpChi[ipl]*fChiNorms[ipl]
            <<std::endl;
        } // 
      } // ipl
      std::cout<<"nRAT is the number of Regions Above Threshold (RAT) used in the study.\n";
      std::cout<<"bCnt is the number of single bumps that were successfully fitted \n";
      std::cout<<"bChi is the average chisq/DOF of the first fit\n";
      std::cout<<"bRMS is the average calculated RMS of the bumps\n";
      std::cout<<"hCnt is the number of RATs that have a single hit\n";
      std::cout<<"hRMS is the average RMS from the Gaussian fit -> use this value for fMinRMS[plane] in the fcl file\n";
      std::cout<<"dTdW is the slope of the track\n";
      std::cout<<"New_ChiNorm is the recommended values of ChiNorm that should be used in the fcl file\n";
      bumpChi.clear();
      bumpRMS.clear();
      bumpCnt.clear();
      RATCnt.clear();
      hitRMS.clear();
      hitCnt.clear();
      loWire.clear();
      loTime.clear();
      hiWire.clear();
      hiTime.clear();
    }
  } // StudyHits


//////////////////////////////////////////////////
  void CCHitFinderAlg::FitStats_t::Reset(unsigned int nGaus) {
    if (nGaus == 0) return;
    MultiGausFits.resize(nGaus);
    std::fill(MultiGausFits.begin(), MultiGausFits.end(), 0);
    FastFits = 0;
  } // CCHitFinderAlg::FitStats_t::Reset()
  
  
  void CCHitFinderAlg::FitStats_t::AddMultiGaus(unsigned int nGaus) {
    ++MultiGausFits[std::min(nGaus, (unsigned int) MultiGausFits.size()) - 1];
  } // CCHitFinderAlg::FitStats_t::AddMultiGaus()
  
  
} // namespace hit

