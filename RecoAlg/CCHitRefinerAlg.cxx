//////////////////////////////////////////////////////////////////////
///
/// CCHitRefiner class
///
/// Bruce Baller, baller@fnal.gov
///
/// Refines hits found by CCHitFinder, using information from
/// ClusterCrawlerAlg
///
////////////////////////////////////////////////////////////////////////


#include <stdint.h>
#include <iostream>
#include <iomanip>

// ROOT Includes
#include "TMinuit.h"

#include "RecoAlg/CCHitRefinerAlg.h"
#include "RecoAlg/CCHitFinderAlg.h"
#include "RecoAlg/ClusterCrawlerAlg.h"

/*
MinuitStruct gMinStruct;

/////////////////////////////////////////
  void fcnW(Int_t &, Double_t *, Double_t &fval,
      double *par, Int_t iflag)
  {
    // Minuit function for fitting the amplitudes of all cluster hits on
    // one Wire + a small adjustment to the hit time
    unsigned short tsize = gMinStruct.WireSignals[0].size();
    if(tsize != gMinStruct.HitSignals[0].size()) {
      mf::LogError("ClusterCrawler")<<"Inconsistent times ";
    }
    unsigned short w = gMinStruct.WireIndex;
    
    // clear the hit signals on this wire
    for(unsigned short t = 0; t< tsize; ++t) {
      gMinStruct.HitSignals[w][t] = 0.;
    }
    
    // now fill the hit signals vector on this wire
    unsigned short ip = 0;
    for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
      // no hit expected on this wire for this cluster?
      if(gMinStruct.vcl[icl].Time[w] < 0) continue;
      float clRMS = gMinStruct.vcl[icl].RMS;
      // add a time offset
      float ptime = gMinStruct.vcl[icl].Time[w] + par[ip+1];
      float ptchk = ptime - 3.5 * clRMS;
      unsigned short ptimeLo = 0;
      if(ptchk > 0.) ptimeLo = ptchk;
      ptchk = ptime + 3.5 * clRMS;
      unsigned short ptimeHi = tsize - 1;
      if(ptchk < ptimeHi) ptimeHi = ptchk;
      for(unsigned short t = ptimeLo; t < ptimeHi; ++t) {
        float arg = (t - ptime) / clRMS;
        gMinStruct.HitSignals[w][t] += par[ip] * exp( -0.5 * arg * arg);
      } // t
      // increment by 2 (Amp, Time offset)
      ip += 2;
    } // icl
    
    // now calculate the difference
    fval = 0.;
    unsigned short nsum = 0.;
    for(unsigned short t = 0; t < tsize; ++t) {
      float arg = (gMinStruct.WireSignals[w][t] 
                 - gMinStruct.HitSignals[w][t]);
      fval += arg * arg;
      nsum += 1;
    } // t
    
    fval = fval / nsum;
    // stash the chisq in the struct
    gMinStruct.fcnVal = fval;
    

    if(iflag == 3) {
      mf::LogVerbatim("ClusterCrawler")<<"fcnW: par "
        <<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]
        <<" fval "<<gMinStruct.fcnVal;
      // print out for debugging
      for(unsigned short t = 0; t < tsize; ++t) {
        mf::LogVerbatim("ClusterCrawler")
          <<t<<std::fixed<<std::setprecision(1)
          <<std::setw(7)<<gMinStruct.WireSignals[w][t]
          <<std::setw(7)<<gMinStruct.HitSignals[w][t]<<std::endl;
      }
    }

  }

/////////////////////////////////////////
  void fcnA(Int_t &, Double_t *, Double_t &fval,
      double *par, Int_t )
  {
    // Minuit function for fitting the vertex position and average hit signals
    // on All wires in the RAT range

    // always calculate fval
    // Fill hit signals using these parameters
    // access the gMinStruct.wireSignals
    unsigned short wsize = gMinStruct.WireSignals.size();
    if(wsize != gMinStruct.HitSignals.size()) {
      mf::LogError("ClusterCrawler")<<"Inconsistent wires ";
    }
    unsigned short tsize = gMinStruct.WireSignals[0].size();
    if(tsize != gMinStruct.HitSignals[0].size()) {
      mf::LogError("ClusterCrawler")<<"Inconsistent times ";
    }
    
    // clear the hitsignals vector
    for(unsigned short w = 0; w < wsize; ++w) {
      for(unsigned short t = 0; t< tsize; ++t) {
        gMinStruct.HitSignals[w][t] = 0.;
      }
    }

    // refill using the input pars for each cluster
    for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
      float clW0  = gMinStruct.vcl[icl].Wire0;
      float clT0  = gMinStruct.vcl[icl].Time0;
      float clRMS = gMinStruct.vcl[icl].RMS;
      float slope = (par[1] - clT0) / (par[0] - clW0);
      for(unsigned short w = 0; w < wsize; ++w) {
        // dead wire?
        if(gMinStruct.WireWght[w] < 0) continue;
        // check for cluster US/DS of the vertex
        short dw = w - (short)par[0];
        // fraction of the cell travelled by the cluster
        float cellFrac = 1.;
        if(gMinStruct.vcl[icl].isDS) {
          // cluster is DS of the vertex
          if(dw < 0) continue;
          if(dw == 0) cellFrac = 1. - par[0] + w;
        } else {
          // cluster is US of the vertex
          if(dw > 0) continue;
          if(dw == 0) cellFrac = par[0] - w;
        }
        float ptime = clT0 + slope * (w - clW0);
        float amp = cellFrac * par[2 + icl];
        // find the +/- 3.5 sigma time range and ensure that it is within
        // the vector bounds
        float ptchk = ptime - 3.5 * clRMS;
        unsigned short ptimeLo = 0;
        if(ptchk > 0.) ptimeLo = ptchk;
        ptchk = ptime + 3.5 * clRMS;
        unsigned short ptimeHi = tsize - 1;
        if(ptchk < ptimeHi) ptimeHi = ptchk;
        for(unsigned short t = ptimeLo; t < ptimeHi; ++t) {
          float arg = (t - ptime) / clRMS;
          gMinStruct.HitSignals[w][t] += amp * exp( -0.5 * arg * arg);
        } // t
      } // w
    } // icl
    
    // now find the difference
    fval = 0.;
    float wsum = 0;
    for(unsigned short w = 0; w < wsize; ++w) {
      // wire weight defined in FitVtxPos
      float wght = gMinStruct.WireWght[w];
      // dead wire?
      if(wght < 0) continue;
      for(unsigned short t = 0; t < tsize; ++t) {
        float arg = (gMinStruct.WireSignals[w][t] 
                   - gMinStruct.HitSignals[w][t]);
        fval += wght * arg * arg;
        wsum += wght;
      } // t
    } // w
    
    fval /= wsum;
    
    // stash the chisq in the struct
    gMinStruct.fcnVal = fval;

  } // fcnA

*/

namespace hit {

//------------------------------------------------------------------------------
  CCHitRefinerAlg::CCHitRefinerAlg(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  void CCHitRefinerAlg::reconfigure(fhicl::ParameterSet const& pset)
  {
    fRefineHits = pset.get<bool >("RefineHits");
    fBEChgRat   = pset.get<float>("BEChgRat");
  }

//------------------------------------------------------------------------------
  CCHitRefinerAlg::~CCHitRefinerAlg()
  {
  }

  
  void CCHitRefinerAlg::RunCCHitRefiner(
      std::vector<recob::Hit>& allhits,
      CCHitFinderAlg::HitCuts& /*hitcuts*/,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx, 
      ClusterCrawlerAlg& fCCAlg) 
  {
    // try to refine hits near vertices. Hits on clusters are assumed to be
    // in reverse wire order, ala ClusterCrawler, i.e. Begin = DS = large wire
    // number and End = US = small wire number.
    // This alg also defines the Begin and End of clusters.


    // copy the hit/cluster association information from CC algorithm
    HitInCluster = fCCAlg.GetHitInCluster();
    // swap Begin/End w/o refining hits?
//    if(!fRefineHits && fBEChgRat > 0.) {
      SetClusterBeginEnd(allhits, tcl);
//      return;
//    }

/* BB June 19. This code is not ready for general use. Turn it off for now.
    
    mf::LogVerbatim("ClusterCrawler")<<"CCHitRefiner ";
    
    unsigned short lastplane = 100;
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].Wght < 0) continue;
  if(iv != 5) continue;
  prt = true;
      theVtx = iv;
      plane = vtx[theVtx].CTP - vtx[theVtx].CTP / 10;
      if(plane != lastplane) {
        fCCAlg.GetHitRange(allhits, vtx[theVtx].CTP, WireHitRange, fFirstWire, fLastWire);
        lastplane = plane;
      }
      // list of clusters that are US (DS) of the vertex
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        // clusters that end DS of the vtx
        if(tcl[icl].BeginVtx == theVtx) clBeg.push_back(icl);
        // clusters that end US of the vtx
        if(tcl[icl].EndVtx == theVtx) clEnd.push_back(icl);
      }
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"Vtx "
    <<theVtx<<" P:W:T "<<plane<<":"<<vtx[iv].Wire<<":"<<(int)vtx[iv].Time;
      for(unsigned int ii = 0; ii <clBeg.size(); ++ii) {
        mf::LogVerbatim("ClusterCrawler")<<" Begin cls "<<tcl[clBeg[ii]].ID;
      }
      for(unsigned int ii = 0; ii <clEnd.size(); ++ii) {
        mf::LogVerbatim("ClusterCrawler")<<" End cls "<<tcl[clEnd[ii]].ID;
      }
      if(clBeg.size() == 0 && clEnd.size() == 0) {
        mf::LogError("ClusterCrawler")
          <<"CCHitRefiner: Found vertex with no associated clusters";
        continue;
      }
      // Find the RAT range = range of wires and times near the vertex
      // where there is a signal
      bool SkipIt = false;
      FindRATRange(allhits, tcl, vtx, SkipIt);
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"RAT range: wires "<<loWire<<" "<<hiWire
    <<" times "<<loTime<<" "<<hiTime<<" SkipIt "<<SkipIt;
      if(SkipIt) continue;
      FillWireSignals(allhits);
      // Fit the cluster parameters at the RAT range boundary
      FillVcl(allhits, tcl, vtx);
      // fit the hit signals to the wire signals
      if(gMinStruct.vcl.size() == 0) continue;
      FitVtxPos(vtx);
      Printvcl(vtx);
      FitHitAmplitudes(vtx);
      Printvcl(vtx);
      // refine existing hits and create new ones
      RefineHits(allhits, tcl, vtx);

    } // iv

    if(fBEChgRat > 0.) SetClusterBeginEnd(allhits, tcl);
    
    gMinStruct.WireSignals.clear();
    gMinStruct.HitSignals.clear();
    gMinStruct.vcl.clear();
*/
  } //RunCCHitFinder

/////////////////////////////////////////
    void CCHitRefinerAlg::SetClusterBeginEnd(
        std::vector<recob::Hit>& /*allhits*/,
        std::vector<ClusterCrawlerAlg::ClusterStore>& tcl)
    {
      // This routine prepares the clusters in tcl for stuffing into
      // recob::cluster. The start and end wire, time and slope are
      // defined based on the ratio of start and end charge. Tracks are
      // assumed to be going "downstream" => increasing wire number =>
      // from the end to the beginning of the tcl hit vector, unless the
      // charge ratio is significantly different
      
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        // ignore deleted clusters
        if(tcl[ii].ID < 0) continue;
        bool GoingDS = true;
        if(tcl[ii].EndChg > fBEChgRat * tcl[ii].BeginChg) GoingDS = false;
        // need to swap the beginning/end?
        if(GoingDS) {
          // slope
          float tmp = tcl[ii].BeginSlp;
          tcl[ii].BeginSlp = tcl[ii].EndSlp;
          tcl[ii].EndSlp = tmp;
          // wire
          unsigned short itmp = tcl[ii].BeginWir;
          tcl[ii].BeginWir = tcl[ii].EndWir;
          tcl[ii].EndWir = itmp;
          // time
          tmp = tcl[ii].BeginTim;
          tcl[ii].BeginTim = tcl[ii].EndTim;
          tcl[ii].EndTim = tmp;
          // charge
          tmp = tcl[ii].BeginChg;
          tcl[ii].BeginChg = tcl[ii].EndChg;
          tcl[ii].EndChg = tmp;
          // vertex
          short jtmp = tcl[ii].BeginVtx;
          tcl[ii].BeginVtx = tcl[ii].EndVtx;
          tcl[ii].EndVtx = jtmp;
        }
      }
    } // SetClusterBeginEnd

/*

/////////////////////////////////////////
  void CCHitRefinerAlg::RefineHits(
    std::vector<recob::Hit>& allhits,
    std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
    std::vector<ClusterCrawlerAlg::VtxStore>& )
  {

    // get the amplitude -> charge normalization from a hit on a cluster
    // in this plane
    unsigned short icl = gMinStruct.vcl[0].tclID;
    unsigned short iht = tcl[icl].tclhits[0];
    float ChgNorm = 2.507 * allhits[iht].PeakAmplitude() * allhits[iht].RMS()
        / allhits[iht].Integral();

    unsigned short wsize = hiWire - loWire + 1;
    for(unsigned short ivcl = 0; ivcl < gMinStruct.vcl.size(); ++ivcl) {
      float rms = gMinStruct.vcl[ivcl].RMS;
      // ID of the cluster in the tcl struct
      short tclID = gMinStruct.vcl[ivcl].tclID;
      if(tclID < 0) {
        // need to make a new cluster
        mf::LogError("ClusterCrawler")<<"Need code to make a new cluster";
        continue;
      }
      bool reMakeCluster = false;
      // vector of tcl hits on every wire within the RAT range
      std::vector<short> hIndx(wsize, -1);
      for(unsigned short ii = 0; ii < tcl[tclID].tclhits.size(); ++ii) {
        unsigned short hit = tcl[tclID].tclhits[ii];
        unsigned short wire = allhits[hit].WireID().Wire;
        if(wire < loWire || wire > hiWire) continue;
        unsigned short w = wire - loWire;
        hIndx[w] = hit;
      } // ii
      // compare hIndx with the vcl struct and decide whether to create,
      // refine or delete hits.
      // vector of IDs of hits within the RAT range. This will be appended
      // to the beginning or end of the tclhits vector
      std::vector<unsigned short> fclhits;
      // start at the DS end of the RAT range to keep cluster hits in order
      short w = wsize - 1;
      while(w >= 0) {
        if(hIndx[w] >= 0 && gMinStruct.vcl[ivcl].Time[w] > 0) {
          // A hit exists on the tcl cluster and one in the vcl as well.
          // Update the hit parameters
          unsigned short theHit = hIndx[w];
          float amp = gMinStruct.vcl[ivcl].Amp[w];
          float amperr = gMinStruct.vcl[ivcl].AmpErr[w];
          recob::Hit const& hit = allhits[theHit];
          allhits[theHit] = recob::Hit(
            hit.Channel(),
            hit.StartTick(),
            hit.EndTick(),
            hit.StartTick() + gMinStruct.vcl[ivcl].Time[w] // peak_time
            0.,
            rms,
            amp,                                           // peak_amplitude
            amperr,                                        // sigma_peak_amplitude
            hit.SummedADC(),
            .507 * amp * rms / ChgNorm,                    // hit_integral
            2.507 * amperr * rms / ChgNorm,                // hit_sigma_integral
            gMinStruct.vcl.size(),                         // multiplicity
            w,                                             // local_index
            gMinStruct.vcl[ivcl].HitChiDOF[w],             // goodness_of_fit
            hit.DegreesOfFreedom()                         // dof FIXME
            hit.View(),
            hit.SignalType(),
            hit.WireID()
            );
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Update hit: ivcl "<<ivcl<<" w "<<w<<" theHit "<<theHit;
          fclhits.push_back(theHit);
        } // updated an existing hit
        else if(hIndx[w] < 0 && gMinStruct.vcl[ivcl].Time[w] > 0) {
          // No hit exists on the tcl cluster. One needs to be added
          unsigned short wire = w + loWire;
          float amp = gMinStruct.vcl[ivcl].Amp[w];
          float amperr = gMinStruct.vcl[ivcl].AmpErr[w];
          recob::Hit newhit; // fix with a constructor
          newhit.Charge = 2.507 * amp * rms / ChgNorm;
          newhit.Amplitude = amp;
          newhit.AmplitudeErr = amperr;
          newhit.Time = loTime + gMinStruct.vcl[ivcl].Time[w];
          newhit.TimeErr = 0.;
          newhit.RMS = rms;
          newhit.RMSErr = 0.;
          newhit.ChiDOF = gMinStruct.vcl[ivcl].HitChiDOF[w];
          // define the wire ID using a hit on this wire
          unsigned short aHit = WireHitRange[wire - fFirstWire].first;
          newhit.Wire = allhits[aHit].Wire;
          newhit.WireNum = wire;
          newhit.numHits = gMinStruct.vcl.size();
          newhit.LoHitID = allhits.size();
          newhit.LoTime = loTime;
          newhit.HiTime = hiTime;
          allhits.push_back(newhit);
          unsigned short theHit = allhits.size() - 1;
          fclhits.push_back(theHit);
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"New hit: ivcl "<<ivcl<<" w "<<w<<" theHit "<<theHit
    <<" Time "<<newhit.Time<<" Amp "<<newhit.Amplitude;
          reMakeCluster = true;
        } // created a new hit
        else if(hIndx[w] >= 0 && gMinStruct.vcl[ivcl].Time[w] < 0) {
          // delete a hit on the cluster
          unsigned short theHit = hIndx[w];
          HitInCluster.makeObsolete(theHit);
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Delete hit: ivcl "<<ivcl<<" w "<<w<<" theHit "<<theHit;
          reMakeCluster = true;
        } // deleted a hit
        --w;
      } // w
      if(reMakeCluster) {
        // hits were added or removed from the cluster in tcl so it needs
        // to be re-made
        cluster::ClusterCrawlerAlg::ClusterStore ncl;
        // Create/copy the stuff that doesn't change
        ncl.ID = tcl.size() + 1;
        // flag the processor code
        ncl.ProcCode = 10000 + tcl[tclID].ProcCode;
        ncl.StopCode = tcl[tclID].StopCode;
        ncl.CTP = tcl[tclID].CTP;
        ncl.BeginVtx = tcl[tclID].BeginVtx;
        ncl.EndVtx = tcl[tclID].EndVtx;
        // vector of new cluster hits
        std::vector<unsigned short> nclhits;
        if(gMinStruct.vcl[ivcl].isDS) {
          ncl.BeginSlp = tcl[tclID].BeginSlp;
          ncl.BeginSlpErr = tcl[tclID].BeginSlpErr;
          ncl.BeginWir = tcl[tclID].BeginWir;
          ncl.BeginTim = tcl[tclID].BeginTim;
          ncl.BeginChg = tcl[tclID].BeginChg;
          ncl.BeginSlp = tcl[tclID].BeginSlp;
          ncl.BeginSlpErr = tcl[tclID].BeginSlpErr;
          ncl.EndSlp = gMinStruct.vcl[ivcl].Slope;
          ncl.EndSlpErr = gMinStruct.vcl[ivcl].SlopeErr;
          // fill with hit IDs DS of the RAT range
          for(unsigned short ii = 0; ii < tcl[tclID].tclhits.size(); ++ii) {
            unsigned short iht = tcl[tclID].tclhits[ii];
            if(allhits[iht].WireID().Wire > hiWire) {
              nclhits.push_back(iht);
              HitInCluster.setCluster(iht, ncl.ID);
            }
          } // ii
          // now append hits inside the RAT range
          for(unsigned short ii = 0; ii < fclhits.size(); ++ii) {
            unsigned short iht = fclhits[ii];
            nclhits.push_back(iht);
            HitInCluster.setCluster(iht, ncl.ID);
          } // ii
  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<" nclhits ";
    for(unsigned short ii = 0; ii < nclhits.size(); ++ii) {
      myprt<<nclhits[ii]<<" ";
    }
    myprt<<"\n";
  }
          unsigned short lastHit = fclhits[fclhits.size() - 1];
          ncl.EndWir = allhits[lastHit].WireID().Wire;
          ncl.EndTim = allhits[lastHit].PeakTime();
        } // cluster is DS
        else {
          ncl.BeginSlp = gMinStruct.vcl[ivcl].Slope;
          ncl.BeginSlpErr = gMinStruct.vcl[ivcl].SlopeErr;
          unsigned short firstHit = fclhits[0];
          ncl.BeginWir = allhits[firstHit].WireID().Wire;
          ncl.BeginTim = allhits[firstHit].PeakTime();
          ncl.BeginChg = allhits[firstHit].Integral();
          ncl.EndSlp = tcl[tclID].EndSlp;
          ncl.EndSlpErr = tcl[tclID].EndSlpErr;
          ncl.EndChg = tcl[tclID].EndChg;
          // put the new/modified hits in the RAT range into nclhits
          for(unsigned short ii = 0; ii < fclhits.size(); ++ii) {
            unsigned short iht = fclhits[ii];
            nclhits.push_back(iht);
            HitInCluster.setCluster(iht, ncl.ID);
          }
          // append the existing tcl hits
          for(unsigned short ii = 0; ii < tcl[tclID].tclhits.size(); ++ii) {
            unsigned short iht = tcl[tclID].tclhits[ii];
            if(allhits[iht].WireID().Wire < loWire) {
              nclhits.push_back(iht);
              HitInCluster.setCluster(iht, ncl.ID);
            }
          } // ii
        } // cluster is not DS
        ncl.tclhits = nclhits;
        tcl.push_back(ncl);
        // declare the old cluster obsolete
        tcl[tclID].ID = -1;
      } // reMakeCluster
    } // ivcl
  } // RefineHits

/////////////////////////////////////////
  void CCHitRefinerAlg::FitHitAmplitudes(
    std::vector<ClusterCrawlerAlg::VtxStore>&)
  {
    // Performs a fit to the wire signal on each wire in the RAT range
    // to find the amplitudes and time offsets of all of the hits

    // max number of parameters is 2x the number of clusters, fitting to
    // the amplitude and a small time offset relative to the line between
    //  vertex (W,T) the the cluster (W0,T0)
    unsigned short mxpar = 2 * gMinStruct.vcl.size();
    if(mxpar == 0) return;

  if(prt) mf::LogVerbatim("ClusterCrawler")<<"FitHitAmplitudes ";
    
    TMinuit *gMin = new TMinuit(mxpar);
    // use the Wire version fcn function
    gMin->SetFCN(fcnW);
    Int_t errFlag = 0;
    
    Double_t arglist[10];
    // turn off printing
    arglist[0] = -1.;
    gMin->mnexcm("SET PRINT", arglist, 1, errFlag);
    // set strategy 0 for faster operation
    arglist[0] = 0.;
    gMin->mnexcm("SET STRATEGY", arglist, 1, errFlag);

    // define the fit vectors
    std::vector<double> par(mxpar);
    std::vector<double> stp(mxpar);
    std::vector<double> parerr(mxpar);
    
    unsigned short wsize = gMinStruct.WireSignals.size();
    for(unsigned short w = 0; w < wsize; ++w) {
      // ignore dead wires
      if(WireHitRange[loWire + w - fFirstWire].first < 0) continue;
      short ip = -1;
      // don't fit if the time separation between hits is too small
      float minSep = 5.;
      for(unsigned short icl = 0; icl < gMinStruct.vcl.size() - 1; ++icl) {
        if(gMinStruct.vcl[icl].Time[w] < 0) continue;
        for(unsigned short jcl = icl + 1; jcl < gMinStruct.vcl.size(); ++jcl) {
          if(gMinStruct.vcl[jcl].Time[w] < 0) continue;
          float Sep = fabs(gMinStruct.vcl[icl].Time[w] - 
                           gMinStruct.vcl[jcl].Time[w]);
          if(Sep < minSep) {
            minSep = Sep;
            break;
          }
        } // jcl
        if(minSep < 5.) break;
      } // icl
      if(minSep < 5.) continue;
      for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
        // the hit times are filled in FitVtxPos. Check for
        // "no hit expected on wire" condition
        if(gMinStruct.vcl[icl].Time[w] < 0) continue;
        ++ip;
        // set the starting amplitude
        par[ip] = gMinStruct.vcl[icl].Amp[w];
        // expect a 30% fluctuation in amplitude wire-to-wire
        stp[ip] = 0.3 * par[ip];
        gMin->mnparm(ip,"", par[ip], stp[ip], 1., 3 * par[ip], errFlag);
        // now set the starting time offset
        ++ip;
        par[ip] = 0.;
        stp[ip] = 1.;
        gMin->mnparm(ip,"", par[ip], stp[ip], -10., 10., errFlag);
  mf::LogVerbatim("ClusterCrawler")<<"mnparm: wire "<<w+loWire<<" icl "<<icl
    <<" time "<<(int)gMinStruct.vcl[icl].Time[w]
    <<" amp "<<(int)gMinStruct.vcl[icl].Amp[w];
      } // icl
      if(ip < 1) continue;
      // pass the wire index to fcnW so it uses the appropriate elements of
      // WireSignals
      gMinStruct.WireIndex = w;
  if(prt) {
    // call fcn with starting parameters
    arglist[0] = 1;
    gMin->mnexcm("CAL1", arglist, 1, errFlag);
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"Wire "<<w<<" fcnW";
    for(unsigned short ii = 0; ii < par.size(); ++ii) myprt<<" par"<<ii
      <<" "<<(int)par[ii];
    myprt<<" chisq "<<gMinStruct.fcnVal;
  }
      // Minuit argument list for Migrad, max calls, tolerance
      arglist[0] = 500; // max calls
      arglist[1] = 1.; // tolerance on fval in fcn
      gMin->mnexcm("MIGRAD", arglist, 2, errFlag);
      if(errFlag != 0) mf::LogError("ClusterCrawler")<<"Bad fit "<<errFlag;
    // get the parameters and stuff them back into the vcl struct
      ip = -1;
      for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
        if(gMinStruct.vcl[icl].Time[w] < 0) continue;
        ++ip;
        gMin->GetParameter(ip, par[ip], parerr[ip]);
        // stuff the amplitude into the vcl struct
        if(par[ip] < 0.) par[ip] = 0.;
        gMinStruct.vcl[icl].Amp[w] = par[ip];
        gMinStruct.vcl[icl].AmpErr[w] = parerr[ip];
        // now adjust the hit time
        ++ip;
        gMin->GetParameter(ip, par[ip], parerr[ip]);
        gMinStruct.vcl[icl].Time[w] += par[ip];
        gMinStruct.vcl[icl].HitChiDOF[w] = gMinStruct.fcnVal;
      } // icl
      if(prt) {
        mf::LogVerbatim myprt("ClusterCrawler");
        myprt<<"Fit done: wire "<<loWire + w;
        ip = -1;
        for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
          if(gMinStruct.vcl[icl].Time[w] < 0) continue;
          ++ip;
          myprt<<" Amp "<<(int)par[ip]<<"+/-"<<(int)parerr[ip];
          ++ip;
          myprt<<" TOff "<<std::fixed<<std::setprecision(1)<<par[ip];
          // zero the time offset since we made the correction above.
          // Only need to do this if calling fcnW to print out Wire and
          // hit signals
          gMin->mnparm(ip,"", 0., 0., -10., 10., errFlag);
        }
        myprt<<" chisq "<<gMinStruct.fcnVal;
        // print out signals after fitting
//        arglist[0] = 3;
//        gMin->mnexcm("CAL1", arglist, 1, errFlag);
      } // prt

    } // w

    delete gMin;

    
  } // FitHitAmp

/////////////////////////////////////////
  void CCHitRefinerAlg::FitVtxPos(
    std::vector<ClusterCrawlerAlg::VtxStore>& vtx)
  {
    // Fit the hit signals to the wire signals by varying the vertex
    // position and hit amplitudes
    
    // Number of fit parameters = 2 for the vertex wire and time plus an
    // amplitude for each wire
    Int_t npar = 2 + gMinStruct.vcl.size();
    // define the starting parameters
    std::vector<double> par(npar);
    std::vector<double> stp(npar);
    std::vector<double> parerr(npar);

    unsigned short wsize = gMinStruct.WireSignals.size();
    unsigned short tsize = gMinStruct.WireSignals[0].size();

    TMinuit *gMin = new TMinuit(npar);
    gMin->SetFCN(fcnA);
    Int_t errFlag = 0;
    
    Double_t arglist[10];
    // turn off printing
    arglist[0] = -1.;
    gMin->mnexcm("SET PRINT", arglist, 1, errFlag);

    // vertex parameters and starting step sizes
    // put the vertex halfway into the vertex wire cell
    par[0] = vtx[theVtx].Wire - loWire + 0.5;
    stp[0] = 0.4;
    gMin->mnparm(0,"", par[0], stp[0], 0., (double)wsize, errFlag);
    par[1] = vtx[theVtx].Time - loTime;
    stp[1] = 3.;
    gMin->mnparm(1,"", par[1], stp[1], 0., (double)tsize, errFlag);
    for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
      unsigned short ip = 2 + icl;
      // find a hit on this cluster to get the amplitude
      float amp = 0.;
      for(unsigned short w = 0; w < wsize; ++w) {
        if(gMinStruct.vcl[icl].Amp[w] > 0.) {
          amp = gMinStruct.vcl[icl].Amp[w];
          break;
        }
      } // w
      par[ip] = amp;
      stp[ip] = 0.3 * par[ip];
      gMin->mnparm(ip,"", par[ip], stp[ip], 0, 0, errFlag);
    }
    
    // define a weight that will be applied to wire/hit signal differences.
    // This scheme gives a higher weight to wire signals that are different
    // from the wire signals on the next DS wire
    gMinStruct.WireWght.resize(wsize);
    std::vector<float> wsum(wsize, 0.);
    for(unsigned short w = 0; w < wsize; ++w) {
      // set WireWght to -1 to flag dead wires
      if(WireHitRange[loWire + w - fFirstWire].first == -1) {
        gMinStruct.WireWght[w] = -1;
      } else {
        // set to -2 to indicate it hasn't been defined
        gMinStruct.WireWght[w] = -2;
      }
      // sum up the ADCs on this wire
      for(unsigned short t = 0; t < tsize; ++t) {
        wsum[w] += fabs(gMinStruct.WireSignals[w][t]);
      }
    } // ww
    // find the normalized difference
    unsigned short w = wsize - 2;
    float minwght = 9999.;
    float maxwght = 0.;
    unsigned short imbig = 0;
    while(w != 0) {
      float wght = fabs(wsum[w] - wsum[w + 1]) / (wsum[w] + wsum[w + 1]);
      wght *= wght;
      gMinStruct.WireWght[w] = wght;
      if(wght < minwght) minwght = wght;
      if(wght > maxwght) {
        maxwght = wght;
        imbig = w;
      }
      --w;
    } // w != 0

    // set the weight to the max value on the wire US of the
    // wire with the max value so that the vertex fit is not unduly biased
    // if it moves one wire US
    short adjBin = imbig - 1;
    if(adjBin >= 0) gMinStruct.WireWght[adjBin] = maxwght;
    // fill in missing weights and scale by the expected wire by wire 
    // dE/dx fluctuations
    for(unsigned short w = 0; w < wsize; ++w) {
      if(gMinStruct.WireWght[w] == -2) gMinStruct.WireWght[w] = minwght;
      float dedxw = 0.;
      if(wsum[w] > 0.) dedxw = 1/(0.3 * wsum[w]);
      gMinStruct.WireWght[w] *= dedxw;
    } 

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"WireWght ";
    for(unsigned short w = 0; w < wsize; ++w) {
      myprt<<std::scientific<<std::setprecision(3)<<gMinStruct.WireWght[w]<<" ";
    }
    myprt<<"\n";
    // call fcn with starting parameters
    arglist[0] = 1;
    gMin->mnexcm("CAL1", arglist, 1, errFlag);
    myprt<<"Starting: par ";
    for(unsigned short ip = 0; ip < par.size(); ++ip) {
      myprt<<" "<<std::fixed<<std::setprecision(2)<<par[ip];
    }
    myprt<<" chisq "<<gMinStruct.fcnVal;
  }

    // set strategy 0 for faster Minuit fitting
    arglist[0] = 0.;
    gMin->mnexcm("SET STRATEGY", arglist, 1, errFlag);

    // execute Minuit command: Migrad, max calls, tolerance
    arglist[0] = 500; // max calls
    arglist[1] = 1.; // tolerance on fval in fcn

    // Ref: http://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node18.html
    gMin->mnexcm("MIGRAD", arglist, 2, errFlag);
    
    // get the parameters
    for(unsigned short ip = 0; ip < par.size(); ++ip) {
      gMin->GetParameter(ip, par[ip], parerr[ip]);
    }

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"Fit done: par ";
    for(unsigned short ip = 0; ip < par.size(); ++ip) {
      myprt<<" "<<std::fixed<<std::setprecision(2)<<par[ip];
    }
    myprt<<" chisq "<<gMinStruct.fcnVal;
  }
    
    // should make a chisq decision here
    
    // update the vertex position
    vtx[theVtx].Wire = loWire + par[0];
    vtx[theVtx].Time = loTime + par[1];
    for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
      // update the cluster slope
      float clW0 = gMinStruct.vcl[icl].Wire0;
      float clT0  = gMinStruct.vcl[icl].Time0;
      float slope = (par[1] - clT0) / (par[0] - clW0);
      gMinStruct.vcl[icl].Slope = slope;
      // stuff the average Amp back into the vcl array
      for(unsigned short w = 0; w < wsize; ++w) {
        // Initialize the Time vector < 0 to flag the "no hit expected"
        // condition
        gMinStruct.vcl[icl].Time[w] = -1;
        gMinStruct.vcl[icl].Amp[w] = 0;
        short dw = w - (short)par[0];
        float cellFrac = 1.;
        if(gMinStruct.vcl[icl].isDS) {
          // vtx is US of the wire
          if(dw < 0) continue;
          // vtx is in the wire cell. Find out if it is <50% of the way
          // through it. If so, don't make a hit on the vertex wire
          if(dw == 0) {
            cellFrac = 1. - par[0] + w;
            if(cellFrac < 0.5) continue;
          }
        } else {
          if(dw > 0) continue;
          if(dw == 0) {
            cellFrac = par[0] - w;
            if(cellFrac < 0.5) continue;
          }
        }
        gMinStruct.vcl[icl].Amp[w] = cellFrac * par[2 + icl];
        // set the correct hit time
        gMinStruct.vcl[icl].Time[w] = clT0 + (w - clW0) * slope;
      } // w
    } // icl

    delete gMin;

  } // FitVtxPos

/////////////////////////////////////////
  void CCHitRefinerAlg::FillVcl(
    std::vector<recob::Hit>& allhits,
    std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
    std::vector<ClusterCrawlerAlg::VtxStore>&)
  {
    // refit the clusters associated with theVtx. The fit is done
    // at the boundary of the RAT range. 
    
    gMinStruct.vcl.clear();
    
    unsigned short wsize = gMinStruct.WireSignals.size();
    // clusters that End at the vertex
    for(unsigned short ii = 0; ii < clEnd.size(); ++ii) {
      unsigned short icl = clEnd[ii];
      // the hit on the boundary of the RAT range
      unsigned short hit0 = 0;
      if(tcl[icl].BeginWir >= hiWire) {
        for(unsigned short jj = tcl[icl].tclhits.size() - 1; jj > 0; --jj) {
          unsigned short hit = tcl[icl].tclhits[jj];
          if(allhits[hit].WireID().Wire >= hiWire) {
            hit0 = hit;
            break;
          }
        } // jj
        MinuitStruct::VtxCluster rcl;
        rcl.tclID = icl;
        // offset the wire and time relative to loWire and hiWire to
        // facilitate fitting
        rcl.Wire0 = allhits[hit0].WireID().Wire - loWire;
        rcl.Time0 = allhits[hit0].PeakTime() - loTime;
        rcl.Slope = tcl[icl].EndSlp;
        rcl.SlopeErr = tcl[icl].EndSlpErr;
        if(rcl.SlopeErr <= 0.) rcl.SlopeErr = 0.1 * fabs(rcl.Slope);
        rcl.Time.resize(wsize);
        rcl.Amp.resize(wsize);
        rcl.AmpErr.resize(wsize);
        rcl.HitChiDOF.resize(wsize);
        for(unsigned short w = 0; w < wsize; ++w) {
          rcl.Amp[w] = allhits[hit0].PeakAmplitude();
          rcl.AmpErr[w] = 0.;
        }
        rcl.RMS   = allhits[hit0].RMS();
        // true for clusters that End at the vertex
        rcl.isDS  = true;
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"End rcl W0 "<<rcl.Wire0<<" T0 "<<(int)rcl.Time0<<" slp "<<rcl.Slope
      <<" RMS "<<rcl.RMS<<" Amp[0] "<<(int)rcl.Amp[0]<<" DS? "<<rcl.isDS;
  } // prt
        gMinStruct.vcl.push_back(rcl);
      } // tcl[icl].BeginWir > hiWire
      else {
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Cluster ends inside the RAT. Deal with it";
      } // tcl[icl].BeginWir < hiWire
    } // ii
    
    // clusters that Begin at the vertex
    for(unsigned short ii = 0; ii < clBeg.size(); ++ii) {
      unsigned short icl = clBeg[ii];
      // the hit on the boundary of the RAT range
      unsigned short hit0 = 0;
      if(tcl[icl].EndWir <= loWire) {
        for(unsigned short jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
          unsigned short hit = tcl[icl].tclhits[jj];
          if(allhits[hit].WireID().Wire <= loWire) {
            hit0 = hit;
            break;
          }
        } // jj
        MinuitStruct::VtxCluster rcl;
        rcl.tclID = icl;
        // offset the wire and time relative to loWire and hiWire to
        // facilitate fitting
        rcl.Wire0 = allhits[hit0].WireID().Wire - loWire;
        rcl.Time0 = allhits[hit0].PeakTime() - loTime;
        // save the slope in case we want to check the deviation between
        // this slope and the one calculated using the vertex fit position
        rcl.Slope = tcl[icl].BeginSlp;
        rcl.SlopeErr = tcl[icl].BeginSlpErr;
        if(rcl.SlopeErr <= 0.) rcl.SlopeErr = 0.1 * fabs(rcl.Slope);
        rcl.Time.resize(wsize);
        rcl.Amp.resize(wsize);
        rcl.AmpErr.resize(wsize);
        rcl.HitChiDOF.resize(wsize);
        for(unsigned short w = 0; w < wsize; ++w) {
          rcl.Amp[w] = allhits[hit0].PeakAmplitude();
          rcl.AmpErr[w] = 0.;
        }
        rcl.RMS   = allhits[hit0].RMS();
        // false for clusters that Begin at the vertex
        rcl.isDS  = false;
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"Beg rcl W0 "<<rcl.Wire0<<" T0 "<<(int)rcl.Time0<<" slp "<<rcl.Slope
      <<" RMS "<<rcl.RMS<<" Amp[0] "<<(int)rcl.Amp[0]<<" DS? "<<rcl.isDS;
  } // prt
        gMinStruct.vcl.push_back(rcl);
      } // tcl[icl].BeginWir > hiWire
      else {
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Cluster ends inside the RAT. Deal with it";
      } // tcl[icl].BeginWir < hiWire
    } // ii
  
  } // FillVcl


/////////////////////////////////////////
  void CCHitRefinerAlg::FillWireSignals(
    std::vector<recob::Hit>& allhits)
  {
    // Fill the wire signals vector and initialize the fitted hit
    // signals vector
    gMinStruct.WireSignals.clear();
    gMinStruct.HitSignals.clear();
    
    // define a zero vector
    std::vector<float> zero(hiTime - loTime + 1);
    
    for(unsigned short wire = loWire; wire <= hiWire; ++wire) {
      // dead wire or no hits. pushback the zero vector
      if(WireHitRange[wire - fFirstWire].first < 0) {
        gMinStruct.WireSignals.push_back(zero);
        continue;
      }
      // find a hit on the wire to get the wire object
      unsigned short iht = WireHitRange[wire - fFirstWire].first;
      art::Ptr< recob::Wire> oWire = allhits[iht].Wire; // TODO ask associations
      // would be nice to just get the wire signal where we need it...
      std::vector<float> signal( oWire->Signal() );
      // temporary RAT vector
      std::vector<float> rat;
      for(unsigned short time = loTime; time <= hiTime; ++time) {
        rat.push_back(signal[time]);
      } // time
      gMinStruct.WireSignals.push_back(rat);
    } // wire
    
    // define the HitSignals vector here as well
    for(unsigned short wire = loWire; wire <= hiWire; ++wire) {
      gMinStruct.HitSignals.push_back(zero);
    }
    
  } // FillWireSignals

/////////////////////////////////////////
  void CCHitRefinerAlg::PrintSignals()
  {
    for(unsigned short wire = loWire; wire <= hiWire; ++wire) {
      unsigned short wIndx = wire - loWire;
      for(unsigned short time = loTime; time <= hiTime; ++time) {
        unsigned short tIndx = time - loTime;
          mf::LogVerbatim("ClusterCrawler")<<"SIG "
            <<wire<<" "<<time<<" "<<(int)gMinStruct.WireSignals[wIndx][tIndx]
            <<" "<<(int)gMinStruct.HitSignals[wIndx][tIndx];
      } // time
    } // wire
  } // PrintSignals


/////////////////////////////////////////
  void CCHitRefinerAlg::FindRATRange(
      std::vector<recob::Hit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
      std::vector<ClusterCrawlerAlg::VtxStore>& ,
      bool& SkipIt)
  {
    // gets the range of wires and times of the RAT surrounding theVtx.
    // SkipIt is set true if there is nothing to do.
    loWire = 9999;
    loTime = 9999;
    hiWire = 0;
    hiTime = 0;
    
    unsigned short nMultgt1 = 0;
    
    // start with clusters that End at the vertex
    for(unsigned short ii = 0; ii < clEnd.size(); ++ii) {
      unsigned short icl = clEnd[ii];
      for(unsigned short jj = tcl[icl].tclhits.size() - 1; jj > 0; --jj) {
        unsigned short hit = tcl[icl].tclhits[jj];
        if(allhits[hit].WireID().Wire < loWire) loWire = allhits[hit].WireID().Wire;
        if(allhits[hit].StartTick() < loTime) loTime = allhits[hit].StartTick();
        if(allhits[hit].WireID().Wire > hiWire) hiWire = allhits[hit].WireID().Wire;
        if(allhits[hit].EndTick() > hiTime) hiTime = allhits[hit].EndTick();
        // stop looking if the multiplicity = 1 
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"End chk "<<tcl[icl].ID<<" "<<allhits[hit].WireID().Wire
    <<":"<<(int)allhits[hit].PeakTime()<<" mult "<<allhits[hit].Multiplicity();
        if(allhits[hit].Multiplicity() > 1) {
          ++nMultgt1;
        } else {
          break;
        }
        // stop if it looks squirrely
        if(jj < tcl[icl].tclhits.size() - 10) break;
      } // jj
    } // ii
    
    // now check clusters that Begin at the vertex
    for(unsigned short ii = 0; ii < clBeg.size(); ++ii) {
      unsigned short icl = clBeg[ii];
      for(unsigned short jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
        unsigned short hit = tcl[icl].tclhits[jj];
        if(allhits[hit].WireID().Wire < loWire) loWire = allhits[hit].WireID().Wire;
        if(allhits[hit].StartTick() < loTime) loTime = allhits[hit].StartTick();
        if(allhits[hit].WireID().Wire > hiWire) hiWire = allhits[hit].WireID().Wire;
        if(allhits[hit].EndTick() > hiTime) hiTime = allhits[hit].EndTick();
        // stop looking if the multiplicity = 1 
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Begin Chk "<<tcl[icl].ID<<" "<<allhits[hit].WireID().Wire
    <<":"<<(int)allhits[hit].PeakTime()<<" mult "<<allhits[hit].Multiplicity();
        if(allhits[hit].Multiplicity() > 1) {
          ++nMultgt1;
        } else {
          break;
        }
        // stop if it looks squirrely
        if(jj > 10) break;
      } // jj
    } // ii
    
    // no sense continuing if there aren't any multiplicity > 1 hits
    // to refine
    if(nMultgt1 == 0 || hiWire == loWire + 1) {
      SkipIt = true;
      return;
    }
    SkipIt = false;
    
    // Expand the range if there is an unused hit just outside the RAT range
    // low end first
    unsigned short low = loWire - 1;
    for(unsigned short wire = low; wire > low - 3; --wire) {
      unsigned short index = wire - fFirstWire;
      // no hits on the wire
      if(WireHitRange[index].first == -2) break;
      // dead wire
      if(WireHitRange[index].first == -1) continue;
      unsigned short firsthit = WireHitRange[index].first;
      unsigned short lasthit = WireHitRange[index].second;
      // require a hit within the time range
      bool gotone = false;
      for(unsigned short hit = firsthit; hit < lasthit; ++hit) {
        // ignore obsolete hit
        if(allhits[hit].Integral() < 0) continue;
        // ignore used hits
        if(HitInCluster.isInCluster(hit)) continue;
        if(allhits[hit].PeakTime() > loTime && allhits[hit].PeakTime() < hiTime) {
          --loWire;
          gotone = true;
          break;
        }
      } // hit
      if(!gotone) break;
    } // wire
    
    // high end next
    unsigned short hiw = hiWire + 1;
    for(unsigned short wire = hiw; wire < hiw + 3; ++wire) {
      unsigned short index = wire - fFirstWire;
      // no hits on the wire
      if(WireHitRange[index].first == -2) break;
      // dead wire
      if(WireHitRange[index].first == -1) continue;
      unsigned short firsthit = WireHitRange[index].first;
      unsigned short lasthit = WireHitRange[index].second;
      // require a hit within the time range
      bool gotone = false;
      for(unsigned short hit = firsthit; hit < lasthit; ++hit) {
        // ignore obsolete hit
        if(allhits[hit].Integral() < 0) continue;
        // ignore used hits
        if(HitInCluster.isInCluster(hit)) continue;
        if(allhits[hit].PeakTime() > loTime && allhits[hit].PeakTime() < hiTime) {
          ++hiWire;
          gotone = true;
          break;
        }
      } // hit
      if(!gotone) break;
    } // wire
    
    // pad the times by some amount
    if(loTime > 5) loTime -= 5;
    if(hiTime < detprop->NumberTimeSamples() - 5) hiTime += 5;

  } // FindRATRange

/////////////////////////////////////////

  void CCHitRefinerAlg::Printvcl(
    std::vector<ClusterCrawlerAlg::VtxStore>& vtx)
  {
    
    unsigned short wsize = gMinStruct.WireSignals.size();
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"Printvcl: vtx: Wire "
      <<std::fixed<<std::setprecision(1)<<vtx[theVtx].Wire
      <<" Time "<<vtx[theVtx].Time<<"\n";
    
    myprt<<"icl Times\n";
    myprt<<"Wire";
    for(unsigned short w = 0; w < wsize; ++w) {
      myprt<<std::setw(8)<<loWire+w;
    }
    myprt<<"\n";
    for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
      myprt<<std::setw(4)<<icl;
      for(unsigned short w = 0; w < wsize; ++w) {
        float time = gMinStruct.vcl[icl].Time[w];
        if(time > 0) time += loTime;
        myprt<<std::setw(8)<<std::fixed<<std::setprecision(1)<<time;
      } // w
      myprt<<"\n";
    } // icl
    
    myprt<<"icl Amps:\n";
    myprt<<"Wire";
    for(unsigned short w = 0; w < wsize; ++w) {
      myprt<<std::setw(8)<<loWire+w;
    }
    myprt<<"\n";
    for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
      myprt<<std::setw(4)<<icl;
      for(unsigned short w = 0; w < wsize; ++w) {
        myprt<<std::setw(8)<<std::fixed<<std::setprecision(1)
          <<gMinStruct.vcl[icl].Amp[w];
      } // w
      myprt<<"\n";
    } // icl
    
    myprt<<"icl Amp Errors:\n";
    myprt<<"Wire";
    for(unsigned short w = 0; w < wsize; ++w) {
      myprt<<std::setw(8)<<loWire+w;
    }
    myprt<<"\n";
    for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
      myprt<<std::setw(4)<<icl;
      for(unsigned short w = 0; w < wsize; ++w) {
        myprt<<std::setw(8)<<std::fixed<<std::setprecision(1)
          <<gMinStruct.vcl[icl].AmpErr[w];
      } // w
      myprt<<"\n";
    } // icl
    
    myprt<<"icl HitChiDOF\n";
    myprt<<"Wire";
    for(unsigned short w = 0; w < wsize; ++w) {
      myprt<<std::setw(8)<<loWire+w;
    }
    myprt<<"\n";
    for(unsigned short icl = 0; icl < gMinStruct.vcl.size(); ++icl) {
      myprt<<std::setw(4)<<icl;
      for(unsigned short w = 0; w < wsize; ++w) {
        myprt<<std::setw(8)<<std::fixed<<std::setprecision(1)
          <<gMinStruct.vcl[icl].HitChiDOF[w];
      } // w
      myprt<<"\n";
    } // icl
  } // Printvcl

*/
} // namespace hit

