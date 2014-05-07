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

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <stdint.h>
#include <iostream>
#include <iomanip>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h"   
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDProducer.h" 


// LArSoft Includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Hit.h"

// ROOT Includes 
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
// #include "TVirtualFitter.h"

#include "RecoAlg/CCHitFinderAlg.h"

namespace cluster{

//------------------------------------------------------------------------------
  CCHitFinderAlg::CCHitFinderAlg(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  void CCHitFinderAlg::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCalDataModuleLabel = pset.get< std::string  >("CalDataModuleLabel");
    fMinSigInd          = pset.get< float       >("MinSigInd");
    fMinSigCol          = pset.get< float       >("MinSigCol");
    fMinRMSInd          = pset.get< float       >("MinRMSInd");
    fMinRMSCol          = pset.get< float       >("MinRMSCol");
    fMaxBumps           = pset.get< unsigned short >("MaxBumps");
    fMaxXtraHits        = pset.get< unsigned short >("MaxXtraHits");
    fChiSplit           = pset.get< float       >("ChiSplit");
    fChiNorms           = pset.get< std::vector< float > >("ChiNorms");
    fTimeOffsets        = pset.get< std::vector< float > >("TimeOffsets");
    fChgNorms           = pset.get< std::vector< float > >("ChgNorms");

    // stuff these parameters into the hitcut struct so they can be accessed
    // by other CC algs
    hitcuts.MinSigInd = fMinSigInd;
    hitcuts.MinSigCol = fMinSigCol;
    hitcuts.MinRMSInd = fMinRMSInd;
    hitcuts.MinRMSCol = fMinRMSCol;
    hitcuts.ChiSplit  = fChiSplit;
    hitcuts.ChiNorms  = fChiNorms;
    hitcuts.TimeOffsets = fTimeOffsets;
    hitcuts.ChgNorms  = fChgNorms;

  }

//------------------------------------------------------------------------------
  CCHitFinderAlg::~CCHitFinderAlg()
  {
  }
  
  void CCHitFinderAlg::RunCCHitFinder(art::Event & evt) {
  
    allhits.clear();

    // make this accessible to ClusterCrawler_module
    art::Handle< std::vector<recob::Wire> > wireVecHandle;
    evt.getByLabel(fCalDataModuleLabel,wireVecHandle);

    unsigned short maxticks = 1000;
    float *ticks = new float[maxticks];
    // define the ticks array used for fitting 
    for(unsigned short ii = 0; ii < maxticks; ++ii) {
      ticks[ii] = ii;
    }
    float *signl = new float[maxticks];
    // initialize the vectors for the hit study
//  StudyHits(0);

//    prt = false;

    for(size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter++){

      art::Ptr<recob::Wire> theWire(wireVecHandle, wireIter);
      theChannel = theWire->Channel();
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


      // minimum number of time samples
      unsigned short minSamples = 2 * minRMS;

      std::vector<geo::WireID> wids = geom->ChannelToWire(theChannel);
      thePlane = wids[0].Plane;
      theWireNum = wids[0].Wire;

      // factor used to normalize the chi/dof fits for each plane
      chinorm = fChiNorms[thePlane];
      timeoff = fTimeOffsets[thePlane];
      ChgNorm = fChgNorms[thePlane];

      // edit this line to debug hit fitting on a particular plane/wire
//      prt = (thePlane == 1 && theWireNum == 839);
      std::vector<float> signal(theWire->Signal());

      unsigned short nabove = 0;
      unsigned short tstart = 0;
      unsigned short maxtime = signal.size() - 2;
      // find the min time when the signal is below threshold
      unsigned short mintime = 3;
      for(unsigned short time = 3; time < maxtime; ++time) {
        if(signal[time] < minSig) {
          mintime = time;
          break;
        }
      }
      for(unsigned short time = mintime; time < maxtime; ++time) {
        if(signal[time] > minSig) {
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
            for(unsigned short ii = tstart; ii < time; ++ii) {
              signl[npt] = signal[ii];
              if(signal[ii    ] > signal[ii - 1] &&
                 signal[ii - 1] > signal[ii - 2] &&
                 signal[ii    ] > signal[ii + 1] &&
                 signal[ii + 1] > signal[ii + 2]) bumps.push_back(npt);
//  if(prt) mf::LogVerbatim("CCHitFinder")<<"signl "<<ii<<" "<<signl[npt];
              ++npt;
            }
// decide if this RAT should be studied
//  StudyHits(1, npt, ticks, signl, tstart);
            // just make a crude hit if too many bumps
            if(bumps.size() > fMaxBumps) {
              MakeCrudeHit(npt, ticks, signl);
              StoreHits(tstart, npt, theWire);
              nabove = 0;
              continue;
            }
            // start looking for hits with the found bumps
            unsigned short nHitsFit = bumps.size();
            unsigned short nfit = 0;
            chidof = 0.;
            bool HitStored = false;
            unsigned short nMaxFit = bumps.size() + fMaxXtraHits;
//  bool first = true;
            while(nHitsFit <= nMaxFit) {
              FitNG(nHitsFit, npt, ticks, signl);
/*
  if(first && SelRAT) {
    first = false;
    StudyHits(2, npt, ticks, signl, tstart);
  }
*/
              // good chisq so store it
              if(chidof < fChiSplit) {
                StoreHits(tstart, npt, theWire);
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
              StoreHits(tstart, npt, theWire);
            }
          } // nabove > minSamples
          nabove = 0;
        } // signal < minSig
      } // time
    } // wireIter

// print out
//  StudyHits(4);

    delete ticks;
    delete signl;

  } //RunCCHitFinder


/////////////////////////////////////////
  void CCHitFinderAlg::FitNG(unsigned short nGaus, unsigned short npt, 
    float *ticks, float *signl)
  {
    // Fit the signal to n Gaussians

    short ndof = npt - 3 * nGaus;
    
    chidof = 9999.;

    if(ndof < 3) return;
    if(bumps.size() == 0) return;

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
    
    TGraph *fitn = new TGraph(npt, ticks, signl);
    TF1 *Gn = new TF1("gn",eqn.c_str());
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
      Gn->SetParameter(index + 2, (double)minRMS);
      Gn->SetParLimits(index + 2, 1., 3*(double)minRMS);
/*
  if(prt) mf::LogVerbatim("CCHitFinder")<<"Bump params "<<ii<<" "<<(short)amp
    <<" "<<(int)bumptime<<" "<<(int)minRMS;
*/
    } // ii bumps

    // search for other bumps that may be hidden by the already found ones
    for(unsigned short ii = bumps.size(); ii < nGaus; ++ii) {
      // bump height must exceed minSig
      float big = minSig;
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
        Gn->SetParameter(index + 2, (double)minRMS);
        Gn->SetParLimits(index + 2, 1., 5*(double)minRMS);
      } // imbig > 0
    } // ii 
    
    // W = set weights to 1, N = no drawing or storing, Q = quiet
    // B = bounded parameters
    fitn->Fit(Gn,"WNQB");
    
    // load the fit into a temp vector
    std::vector<double> partmp;
    std::vector<double> partmperr;

    for(unsigned short ipar = 0; ipar < 3 * nGaus; ++ipar) {
      partmp.push_back(Gn->GetParameter(ipar));
      partmperr.push_back(Gn->GetParError(ipar));
    }
    chidof = Gn->GetChisquare() / ( ndof * chinorm);

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
      if(partmp[index] < minSig) {
        fitok = false;
        break;
      }
      // ensure that the RMS is large enough but not too large
      float rms = partmp[index + 2];
      if(rms < 0.5 * minRMS || rms > 5 * minRMS) {
        fitok = false;
        break;
      }
      // ensure that the hits are not too similar in time (< 2 ticks)
      for(unsigned short jj = 0; jj < nGaus; ++jj) {
        if(jj == ii) continue;
        unsigned short jndex = jj * 3;
        float timediff = fabs(partmp[jndex + 1] - partmp[index + 1]);
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
//      if(prt) mf::LogVerbatim("CCHitFinder")<<"Bad fit parameters";
    }
    
    delete fitn;
    delete Gn;
    
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
    rms = sqrt(rms / sumS);
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
    float meanerr = sqrt(1/sumS);
    float rmserr = 0.2 * rms;
    parerr.push_back(amperr);
    parerr.push_back(meanerr);
    parerr.push_back(rmserr);
/*
  if(prt) mf::LogVerbatim("CCHitFinder")<<" errors Amp "<<amperr<<" mean "
    <<meanerr<<" rms "<<rmserr;
*/
    chidof = 9999.;
  } // MakeCrudeHit


/////////////////////////////////////////
  void CCHitFinderAlg::StoreHits(unsigned short TStart, unsigned short npt,
    art::Ptr<recob::Wire>& theWire)
  {
    // store the hits in the struct
    unsigned short nhits = par.size() / 3;
    
    if(nhits == 0) return;

  // fill RMS for single hits
//  StudyHits(3);

    float loTime = TStart;
    float hiTime = TStart + npt;
    // check for large separation between hits. These vectors define the
    // boundaries of sub-multiplets
    std::vector<float> loTimes(nhits, loTime);
    std::vector<float> hiTimes(nhits, hiTime);
    std::vector<unsigned short> loHitIDs(nhits,allhits.size());
    std::vector<unsigned short> nMultHits(nhits,nhits);
    unsigned short nhm = 0;
    // scan for large separation
    for(unsigned short hit = 1; hit < nhits; ++hit) {
      // increment the number of hits in this (sub-)multiplet
      ++nhm;
      unsigned short index = 3 * hit;
      // RMS of the previous hit
      float rms1 = par[index - 1];
      // significance of the time separation from the previous hit
      float sep = fabs(par[index + 1] - par[index - 2]) / rms1;
      if(sep > 5) {
        // large separation. re-define the boundaries
        float newLoTime = 0.5 * (par[index + 1] + par[index - 2]);
        // re-define the previous hit boundaries
        for(unsigned short phit = 0; phit < hit; ++phit) {
          // hi times for the previous hits are the new low times for the
          // next hits
          hiTimes[phit] = newLoTime;
          nMultHits[phit] = nhm;
        } // phit
        // reset the hit counter
        nhm = 0;
        // re-define the lower boundary of the next set of hits
        for(unsigned short phit = hit; phit < nhits; ++phit) {
          loTimes[phit] = newLoTime;
          loHitIDs[phit] = hit;
          // set the hit multiplicity = 0 as a flag
          nMultHits[phit] = 0;
        } // phit
      } // sep > 5
    } // hit
    // correct the high boundary for the last set of hits if necessary
    if(nhm < nhits) {
      for(unsigned short hit = nhits - 1; hit > 0; --hit) {
        if(nMultHits[hit] == 0) {
          nMultHits[hit] = nhm;
        } else {
          break;
        }
      }
    }
/*
  if(prt) {
    mf::LogVerbatim("CCHitFinder")<<"hit loTime hiTime loHitIDs nMultHits";
    for(unsigned short hit = 0; hit < nhits; ++hit) {
      mf::LogVerbatim("CCHitFinder")<<hit<<" "<<(int)loTimes[hit]
        <<" "<<(int)hiTimes[hit]
        <<" "<<loHitIDs[hit]<<" "<<nMultHits[hit];
    }
  }
*/
    CCHit onehit;
    // lohitid is the index of the first hit that will be added. Hits with
    // Multiplicity > 1 will reside in a block from
    // lohitid to lohitid + numHits - 1
    unsigned short lohitid = allhits.size();
    for(unsigned short hit = 0; hit < nhits; ++hit) {
      unsigned short index = 3 * hit;
      onehit.Charge = Sqrt2Pi * par[index] * par[index + 2] / ChgNorm;
      onehit.ChargeErr = SqrtPi * (parerr[index] * par[index + 2] +
                                   par[index] * parerr[index + 2]);
      onehit.Amplitude = par[index];
      onehit.AmplitudeErr = parerr[index];
      onehit.Time = par[index + 1] + TStart + timeoff;
      onehit.TimeErr = parerr[index + 1];
      onehit.RMS = par[index + 2];
      onehit.RMSErr = parerr[index + 2];
      onehit.ChiDOF = chidof;
      onehit.Wire = theWire;
      onehit.WireNum = theWireNum;
      onehit.numHits = nhits;
      onehit.LoHitID = lohitid;
      onehit.LoTime = loTime;
      onehit.HiTime = hiTime;
      // set flag indicating hit is not used in a cluster
      onehit.InClus = 0;
/*
  if(prt) {
    mf::LogVerbatim("CCHitFinder")<<"W:T "<<theWireNum<<":"<<(short)onehit.Time
      <<" Chg "<<(short)onehit.Charge
      <<" RMS "<<onehit.RMS
      <<" lo ID "<<onehit.LoHitID
      <<" numHits "<<nhm
      <<" loTime "<<loTime<<" hiTime "<<hiTime
      <<" chidof "<<chidof;
  }
*/
      allhits.push_back(onehit);
    } // hit
  } // StoreHits

/*
//////////////////////////////////////////////////
  void CCHitFinderAlg::StudyHits(unsigned short flag, unsigned short npt,
      float *ticks, float *signl, unsigned short tstart) {
    // study hits

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
      // decide if this RAT should be studied. look for a large PH 
      SelRAT = false;
      for(unsigned short ii = 0; ii < npt; ++ii) {
        if(signl[ii] > 20.) {
          SelRAT = true;
          RATCnt[thePlane] += 1;
          break;
        }
      }  // ii
      return;
    } // flag == 1
    
    if(flag == 2) {
      if(!SelRAT) return;
      // in this section we find the low/hi wire/time for a large PH
      // signal (e.g. proton track). This will be used to calculate
      // the slope dT/dW to study hit width, fraction of crude hits, etc
      // vs dT/dW
      // find the peak value and time of the peak value
      float big = 0.;
      float imbig = 0.;
      for(unsigned short ii = 0; ii < npt; ++ii) {
        if(signl[ii] > big) {
          big = signl[ii];
          imbig = ii;
        }
      } // ii
      // require a significant PH 
      if(big > 20) {
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
      } // big > 20
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
        bumpRMS[thePlane] += sqrt(sumt / sum);
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
      std::cout<<" ipl nRAT bCnt   bChi   bRMS 1hCnt   1hRMS  Theta New_ChiNorm"<<std::endl;
      for(unsigned short ipl = 0; ipl < 3; ++ipl) {
        if(bumpCnt[ipl] > 0) {
          bumpChi[ipl] = bumpChi[ipl] / (float)bumpCnt[ipl];
          bumpRMS[ipl] = bumpRMS[ipl] / (float)bumpCnt[ipl];
          hitRMS[ipl]  = hitRMS[ipl]  / (float)hitCnt[ipl];
          // calculate the slope
          float dTdW = fabs((hiTime[ipl] - loTime[ipl]) / (hiWire[ipl] - loWire[ipl]));
          // scale factor is for MicroBooNE 
          int theta = atan(0.273 * dTdW) * 180. / 3.142;
          std::cout<<"BB "<<ipl<<std::right<<std::setw(5)<<RATCnt[ipl]
            <<std::setw(5)<<bumpCnt[ipl]
            <<std::setw(7)<<std::fixed<<std::setprecision(2)<<bumpChi[ipl]
            <<std::setw(7)<<bumpRMS[ipl]
            <<std::setw(7)<<hitCnt[ipl]
            <<std::setw(7)<<std::setprecision(1)<<hitRMS[ipl]
            <<std::setw(7)<<theta
            <<std::setw(7)<<std::setprecision(1)
            <<bumpChi[ipl]*fChiNorms[ipl]
            <<std::endl;
        } // 
      } // ipl
      std::cout<<"Set MinRMSInd and MinRMSCol in the fcl file to the "
               <<"values of hRMS printed above"<<std::endl;
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
*/

} // namespace cluster

