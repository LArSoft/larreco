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


extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <stdint.h>
#include <iostream>
#include <iomanip>

// ROOT Includes
#include "TMinuit.h"

#include "RecoAlg/CCHitRefinerAlg.h"
#include "RecoAlg/CCHitFinderAlg.h"
#include "RecoAlg/ClusterCrawlerAlg.h"

MinuitStruct gMinuitStruct;

/////////////////////////////////////////
  void fcn(Int_t &npar, Double_t *gin, Double_t &fval,
      double *par, Int_t iflag)
  {
    // Minuit function for fitting the vertex position and hit signals
    // in the RAT range

    if(iflag == 1) {
      // initialize
    } // iflag == 1

    // always calculate fval
    // Fill hit signals using these parameters
    // access the gMinuitStruct.wireSignals
    unsigned short wsize = gMinuitStruct.WireSignals.size();
    if(wsize != gMinuitStruct.HitSignals.size()) {
      mf::LogError("ClusterCrawler")<<"Inconsistent wires ";
    }
    unsigned short tsize = gMinuitStruct.WireSignals[0].size();
    if(tsize != gMinuitStruct.HitSignals[0].size()) {
      mf::LogError("ClusterCrawler")<<"Inconsistent times ";
    }

    // ensure that the parameters are within the bounds of the RAT range;
//    fval = 1.E6;
//    if(par[0] < 0. || (unsigned short)par[0] > wsize - 1) return;
//    if(par[1] < 0. || (unsigned short)par[1] > tsize - 1) return;
    
    // clear the hitsignals vector
    for(unsigned short w = 0; w < wsize; ++w) {
      for(unsigned short t = 0; t< tsize; ++t) {
        gMinuitStruct.HitSignals[w][t] = 0.;
      }
    }
/*
  mf::LogVerbatim("ClusterCrawler")<<"par ";
  for(int ip = 0; ip < 4; ++ip) {
    mf::LogVerbatim("ClusterCrawler")<<par[ip]<<" ";
  }
*/
    // fill the hit signals vector using the input pars for each cluster
    for(unsigned short icl = 0; icl < gMinuitStruct.vcl.size(); ++icl) {
      float clW0  = gMinuitStruct.vcl[icl].Wire0;
      float clT0  = gMinuitStruct.vcl[icl].Time0;
      float clRMS = gMinuitStruct.vcl[icl].RMS;
      float slope = (par[1] - clT0) / (par[0] - clW0);
      for(unsigned short w = 0; w < wsize; ++w) {
        // check for cluster US/DS of the vertex
        short dw = w - (short)par[0];
        // fraction of the cell travelled by the cluster
        float cellFrac = 1.;
        if(gMinuitStruct.vcl[icl].isDS) {
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
//  mf::LogVerbatim("ClusterCrawler")<<"ptime "<<icl<<" wire "<<w<<" "<<ptime
//    <<" amp "<<amp<;
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
          gMinuitStruct.HitSignals[w][t] = amp * exp( -0.5 * arg * arg);
        } // t
      } // w
    } // icl
    
    // now calculate the difference
    fval = 0.;
    float wsum = 0;
    for(unsigned short w = 0; w < wsize; ++w) {
      for(unsigned short t = 0; t < tsize; ++t) {
        if(gMinuitStruct.WireSignals[w][t] != 0.) {
          // assume a 30% rms error on the wire signals
          float wght = 1. / (0.3 * fabs(gMinuitStruct.WireSignals[w][t]));
          float arg = (gMinuitStruct.WireSignals[w][t] 
                     - gMinuitStruct.HitSignals[w][t]);
          fval += wght * arg * arg;
          wsum += wght;
//  mf::LogVerbatim("ClusterCrawler")<<"SIG "<<w<<" "<<t<<" "<<gMinuitStruct.WireSignals[w][t]
//    <<" "<<gMinuitStruct.HitSignals[w][t];
        }
      } // t
    } // w
    
    fval /= wsum;
    // stash the chisq in the struct
    gMinuitStruct.MatchChisq = fval;
    
/*
    if(iflag == 3) {
      // last fit finished
    }
*/

  } // fcn



namespace cluster{

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
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      CCHitFinderAlg::HitCuts& hitcuts,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx, 
      ClusterCrawlerAlg& fCCAlg) 
  {
    // try to refine hits near vertices. Hits on clusters are assumed to be
    // in reverse wire order, ala ClusterCrawler, i.e. Begin = DS = large wire
    // number and End = US = small wire number.
    // This alg also defines the Begin and End of clusters.


    // swap Begin/End w/o refining hits?
    if(!fRefineHits && fBEChgRat > 0.) {
      SetClusterBeginEnd(allhits, tcl);
      return;
    }

    mf::LogVerbatim("ClusterCrawler")<<"CCHitRefiner ";
    
    unsigned short lastplane = 100;
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].Wght < 0) continue;
  if(iv != 4) continue;
  prt = true;
      // move the vertex into the middle of the wire cell
      vtx[iv].Wire += 0.5;
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
      FindRATRange(allhits, tcl, vtx);
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"RAT range: wires "<<loWire<<" "<<hiWire
    <<" times "<<loTime<<" "<<hiTime;
      FillWireSignals(allhits);
      // Fit the cluster parameters at the RAT range boundary
      RefitClusters(allhits, tcl, vtx, fCCAlg);
      // fit the hit signals to the wire signals
      if(gMinuitStruct.vcl.size() == 0) continue;
      FitHitSignals(vtx);
      // refine existing hits and create new ones
      RefineHits(allhits, tcl, vtx);
    } // iv

    if(fBEChgRat > 0.) SetClusterBeginEnd(allhits, tcl);
    
    gMinuitStruct.WireSignals.clear();
    gMinuitStruct.HitSignals.clear();
    gMinuitStruct.vcl.clear();

  } //RunCCHitFinder

/////////////////////////////////////////
  void CCHitRefinerAlg::RefineHits(
    std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
    std::vector<ClusterCrawlerAlg::VtxStore>& vtx)
  {

    // get the amplitude -> charge normalization from a hit on a cluster
    // in this plane
    unsigned short icl = gMinuitStruct.vcl[0].clIndx;
    unsigned short iht = tcl[icl].tclhits[0];
    float ChgNorm = 2.507 * allhits[iht].Amplitude * allhits[iht].RMS
        / allhits[iht].Charge;

  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"vtx chk "<<vtx[theVtx].Wire<<" "<<vtx[theVtx].Time;

    // Refine existing hits and create new ones
    for(unsigned short ii = 0; ii < gMinuitStruct.vcl.size(); ++ii) {
      // cluster ID in the tcl struct
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"vcl chk "<<gMinuitStruct.vcl[ii].Wire0<<" "
    <<gMinuitStruct.vcl[ii].Time0;
      float slope = (vtx[theVtx].Time - loTime - gMinuitStruct.vcl[ii].Time0) / 
                    (vtx[theVtx].Wire - loWire - gMinuitStruct.vcl[ii].Wire0);
  if(prt) mf::LogVerbatim("ClusterCrawler")<<"slope "<<slope;
      float rms = gMinuitStruct.vcl[ii].RMS;
      float amp = gMinuitStruct.vcl[ii].Amp;
      float amperr = gMinuitStruct.vcl[ii].AmpErr;
      unsigned short icl = gMinuitStruct.vcl[ii].clIndx;
      // set true if the existing cluster needs to be replaced
      bool reMakeCluster = false;
      // index of the hit on each wire (if one exists) in the RAT range
      // -2 = no hit expected on the wire
      // -1 = hit expected on the wire but none exists
      // >= 0 hit index
      unsigned short nIndx = hiWire - loWire + 1;
      std::vector<short> hIndx(nIndx, -2);
      unsigned short vWire = vtx[theVtx].Wire - loWire;
      if(gMinuitStruct.vcl[ii].isDS) {
        // set the indices if a hit is expected on a wire
        for(unsigned short jj = vWire; jj < nIndx; ++jj) {
          hIndx[jj] = -1;
        } // jj
        // next fill in the hit indices
        for(unsigned short jj = tcl[icl].tclhits.size()-1; jj > 0; --jj) {
          unsigned short iht = tcl[icl].tclhits[jj];
          // wire index in the RAT range hIndx vector
          short wIndx = allhits[iht].WireNum - loWire;
          if(wIndx < 0) continue;
          // delete any hits that are US of the vertex wire and flag this
          // cluster for re-making
          if(wIndx < vWire) {
            allhits[iht].InClus = -1;
            reMakeCluster = true;
            continue;
          }
          if(wIndx == nIndx) break;
          // a hit exists
          hIndx[wIndx] = iht;
        } // jj
      } else {
        // cluster is US of the vertex
        // Note: assume no hit on the vertex wire for US clusters...
        // set the "cluster expected on wire" condition
        for(unsigned short jj = 0; jj < vWire; ++jj) {
          hIndx[jj] = -1;
        } // jj
        // next fill in the hit indices
        for(unsigned short jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
          unsigned short iht = tcl[icl].tclhits[jj];
          short wIndx = allhits[iht].WireNum - loWire;
          if(wIndx < 0) break;
          // check for hits DS of the vertex
          if(allhits[iht].WireNum - loWire >= vWire) {
            allhits[iht].InClus = -1;
            reMakeCluster = true;
            continue;
          }
          if(wIndx == nIndx) continue;
          hIndx[wIndx] = iht;
        } // jj
      } // !isDS

  if(prt) {
    mf::LogVerbatim myprt("ClusterCrawler");
    myprt<<"hIndx ";
    for(unsigned short w = 0; w < hIndx.size(); ++w) {
    myprt<<hIndx[w]<<" ";
    }
  }
      // now refine/create hits
      for(unsigned short w = 0; w < hIndx.size(); ++w) {
        // ensure that this isn't a dead wire
        if(hIndx[w] == -2) continue;
        unsigned short wire = loWire + w;
        float dwire = w - gMinuitStruct.vcl[ii].Wire0;
        float hTime = loTime + gMinuitStruct.vcl[ii].Time0 + slope * dwire;
        float cellFrac = 1.;
        // scale amplitude by the cell fraction
        if(w == vWire) {
          if(gMinuitStruct.vcl[ii].isDS) {
            cellFrac = 1. - vtx[theVtx].Wire + int(vtx[theVtx].Wire);
          } else {
            cellFrac = vtx[theVtx].Wire - int(vtx[theVtx].Wire);
          }
        }
        unsigned short theHit = hIndx[w];
        if(hIndx[w] == -1) {
          // make a new hit
          cluster::CCHitFinderAlg::CCHit newhit;
          theHit = allhits.size();
          // initialize the struct
          newhit.Charge = 0.; newhit.ChargeErr = 0.;
          newhit.Amplitude = 0.; newhit.AmplitudeErr = 0.;
          newhit.Time = 0.; newhit.TimeErr = 0.;
          newhit.RMS = 0.; newhit.RMSErr = 0.;
          newhit.ChiDOF = 0.; newhit.numHits = 0;
          // define the wire ID using a hit on this wire
          unsigned short aHit = WireHitRange[wire - fFirstWire].first;
          newhit.Wire = allhits[aHit].Wire;
          newhit.WireNum = wire;
          newhit.LoHitID = theHit; newhit.LoTime = 0; newhit.HiTime = 0;
          newhit.InClus = 0;
          allhits.push_back(newhit);
          hIndx[w] = theHit;
          // Need to re-make the cluster since we are adding hits to it
          reMakeCluster = true;
        } // hIndx[w] == -1
        // theHit now points to a valid hit. Update the parameters
        allhits[theHit].Charge = 2.507 * amp * rms / ChgNorm;
        allhits[theHit].ChargeErr = 2.507 * amperr * rms / ChgNorm;
        allhits[theHit].Amplitude = cellFrac * amp;
        allhits[theHit].AmplitudeErr = amperr;
        allhits[theHit].Time = hTime;
        allhits[theHit].TimeErr = 0.;
        allhits[theHit].RMS = rms;
        allhits[theHit].RMSErr = 0.;
        allhits[theHit].ChiDOF = 0.;
        allhits[theHit].numHits = 1;
        allhits[theHit].LoHitID = theHit;
        allhits[theHit].LoTime = hTime - 3 * rms;
        allhits[theHit].HiTime = hTime + 3 * rms;
        allhits[theHit].InClus = tcl[icl].ID;
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"RAT range hit "<<hIndx[w]<<" wire "<<wire;
      } // w
      if(reMakeCluster) {
        // hits were added or removed from the cluster in tcl so it needs
        // to be re-made. 
        std::vector<unsigned short> fcl2hits;
        // Start by copying the hits DS of the RAT range
        // add a slot for a new cluster
        cluster::ClusterCrawlerAlg::ClusterStore ncl;
        unsigned short newClID = tcl.size() + 1;
        ncl.ID = newClID;
        ncl.ProcCode = 10000 + tcl[icl].ProcCode;
        ncl.Assn = -1;
        ncl.StopCode = tcl[icl].StopCode;
        ncl.CTP = tcl[icl].CTP;
        if(gMinuitStruct.vcl[ii].isDS) {
          ncl.BeginSlp = tcl[icl].BeginSlp;
          ncl.BeginSlpErr = tcl[icl].BeginSlpErr;
          ncl.BeginWir = tcl[icl].BeginWir;
          ncl.BeginTim = tcl[icl].BeginTim;
          ncl.BeginChg = tcl[icl].BeginChg;
          ncl.EndSlp = gMinuitStruct.vcl[ii].Slope;
          ncl.EndSlpErr = 0.;
          ncl.EndWir = vtx[theVtx].Wire;
          ncl.EndTim = vtx[theVtx].Time;
          ncl.EndChg = 2.507 * amp * rms / ChgNorm;
          // fill the fcl2hits vector in the proper order.
          // start with the DS hits
          for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
            unsigned short iht = tcl[icl].tclhits[ii];
            if(allhits[iht].WireNum > hiWire) {
              fcl2hits.push_back(tcl[icl].tclhits[ii]);
              allhits[iht].InClus = newClID;
            }
          } // ii
          // now stuff in the new hits in the RAT range
          unsigned short wire = hiWire;
          while(wire >= loWire) {
            unsigned short w = wire - loWire;
            if(hIndx[w] >= 0) {
              fcl2hits.push_back(hIndx[w]);
              allhits[hIndx[w]].InClus = newClID;
            }
            --wire;
          } // wire >= loWire
        } else {
          ncl.BeginSlp = gMinuitStruct.vcl[ii].Slope;
          ncl.BeginSlpErr = 0.;
          ncl.BeginWir = vtx[theVtx].Wire;
          ncl.BeginTim = vtx[theVtx].Time;
          ncl.BeginChg = 2.507 * amp * rms / ChgNorm;
          ncl.EndSlp = tcl[icl].EndSlp;
          ncl.EndSlpErr = tcl[icl].EndSlpErr;
          ncl.EndWir = tcl[icl].EndWir;
          ncl.EndTim = tcl[icl].EndTim;
          ncl.EndChg = tcl[icl].EndChg;
          // Stuff in the new hits in the RAT range
          unsigned short wire = hiWire;
          while(wire >= loWire) {
            unsigned short w = wire - loWire;
            if(hIndx[w] >= 0) {
              fcl2hits.push_back(hIndx[w]);
              allhits[hIndx[w]].InClus = newClID;
            }
            --wire;
          } // wire >= loWire
          for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
            unsigned short iht = tcl[icl].tclhits[ii];
            if(allhits[iht].WireNum < loWire) {
              fcl2hits.push_back(tcl[icl].tclhits[ii]);
              allhits[iht].InClus = newClID;
            } else {
              // declare the hit obsolete in the RAT range
              allhits[iht].InClus = -1;
            } // allhits[iht].WireNum >
          } // ii
        }
        ncl.BeginVtx = tcl[icl].BeginVtx;
        ncl.EndVtx = tcl[icl].EndVtx;
        ncl.tclhits = fcl2hits;
/*
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")<<"fcl2hits ";
    for(unsigned short jj = 0; jj < fcl2hits.size(); ++jj) {
      mf::LogVerbatim("ClusterCrawler")<<fcl2hits[jj]<<" "<<allhits[fcl2hits[jj]].WireNum;
    }
*/
        tcl.push_back(ncl);
        // declare the old cluster obsolete
        tcl[icl].ID = -1;
      } // reMakeCluster
    } // ii (cluster)
  } // RefineHits


/////////////////////////////////////////
  void CCHitRefinerAlg::FitHitSignals(
    std::vector<ClusterCrawlerAlg::VtxStore>& vtx)
  {
    // Fit the hit signals to the wire signals by varying the vertex
    // position and hit amplitudes
    
    // Number of fit parameters = 2 for the vertex wire and time plus an
    // amplitude for each wire
    Int_t npar = 2 + gMinuitStruct.vcl.size();
    // define the starting parameters
    std::vector<double> par(npar);
    std::vector<double> stp(npar);
    std::vector<double> parerr(npar);

    TMinuit *gMin = new TMinuit(npar);
    gMin->SetFCN(fcn);
    Int_t errFlag = 0;
    
    Double_t arglist[10];
    // turn off printing
    arglist[0] = -1.;
    gMin->mnexcm("SET PRINT", arglist, 1, errFlag);

    // vertex parameters and starting step sizes
    par[0] = vtx[theVtx].Wire - loWire;
    stp[0] = 0.4;
    gMin->mnparm(0,"", par[0], stp[0], 0, 0, errFlag);
    par[1] = vtx[theVtx].Time - loTime;
    stp[1] = 3.;
    gMin->mnparm(1,"", par[1], stp[1], 0, 0, errFlag);
    for(unsigned short icl = 0; icl < gMinuitStruct.vcl.size(); ++icl) {
      unsigned short ip = 2 + icl;
      par[ip] = gMinuitStruct.vcl[icl].Amp;
      stp[ip] = 0.3 * par[ip];
      gMin->mnparm(ip,"", par[ip], stp[ip], 0, 0, errFlag);
    }
    
    // set strategy 0 for faster operation
    arglist[0] = 0.;
    gMin->mnexcm("SET STRATEGY", arglist, 1, errFlag);

    // execute Minuit command: Migrad, max calls, tolerance
    arglist[0] = 500; // max calls
    arglist[1] = 1.; // tolerance on fval in fcn

    // Ref: http://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node18.html
    gMin->mnexcm("MIGRAD", arglist, 2, errFlag);
    
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Chisq "<<gMinuitStruct.MatchChisq;
    // get the parameters
    for(unsigned short ip = 0; ip < par.size(); ++ip) {
      gMin->GetParameter(ip, par[ip], parerr[ip]);
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Par "<<ip<<" "<<par[ip];
    }
    
    // should make a chisq decision here
    
    // update the vertex position
    vtx[theVtx].Wire = loWire + par[0];
    vtx[theVtx].Time = loTime + par[1];
    // stuff the average Amp back into the vcl array
    for(unsigned short icl = 0; icl < gMinuitStruct.vcl.size(); ++icl) {
      gMinuitStruct.vcl[icl].Amp = par[2 + icl];
      gMinuitStruct.vcl[icl].AmpErr = parerr[2 + icl];
    }

    delete gMin;

  } // FitHitSignals

/////////////////////////////////////////
  void CCHitRefinerAlg::RefitClusters(
    std::vector<CCHitFinderAlg::CCHit>& allhits,
    std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
    std::vector<ClusterCrawlerAlg::VtxStore>& vtx,
    ClusterCrawlerAlg& fCCAlg)
  {
    // refit the clusters associated with theVtx. The fit is done
    // at the boundary of the RAT range. This routine is only called during
    // setup and may be replaced completely by fcn
    
    gMinuitStruct.vcl.clear();
    
    for(unsigned short ii = 0; ii < clEnd.size(); ++ii) {
      unsigned short icl = clEnd[ii];
      // the hit on the boundary of the RAT range
      unsigned short hit0 = 0;
      if(tcl[icl].BeginWir >= hiWire) {
        for(unsigned short jj = tcl[icl].tclhits.size() - 1; jj > 0; --jj) {
          unsigned short hit = tcl[icl].tclhits[jj];
          if(allhits[hit].WireNum >= hiWire) {
            hit0 = hit;
            break;
          }
        } // jj
        // fit the 3 points DS of the RAT range. Also get the average
        // width and amplitude
        fCCAlg.FitClusterMid(allhits, tcl, icl, hit0, -3);
        MinuitStruct::VtxCluster rcl;
        rcl.clIndx = icl;
        // offset the wire and time relative to loWire and hiWire to
        // facilitate fitting
        rcl.Wire0 = allhits[hit0].WireNum - loWire;
        rcl.Time0 = fCCAlg.clpar[0] - loTime;
        // save the slope in case we want to check the deviation between
        // this slope and the one calculated using the vertex fit position
        rcl.Slope = fCCAlg.clpar[1];
        rcl.RMS   = fCCAlg.fAveRMS;
        // true for clusters that End at the vertex
        rcl.isDS  = true;
        rcl.Amp   = fCCAlg.fAveAmp;
  if(prt) {
    mf::LogVerbatim("ClusterCrawler")
      <<"rcl "<<rcl.Wire0<<" T0 "<<(int)rcl.Time0<<" slp "<<rcl.Slope
      <<" RMS "<<rcl.RMS<<" Amp "<<(int)rcl.Amp<<" DS? "<<rcl.isDS;
  } // prt
        gMinuitStruct.vcl.push_back(rcl);
      } // tcl[icl].BeginWir > hiWire
      else {
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Cluster ends inside the RAT. Deal with it";
      } // tcl[icl].BeginWir < hiWire
    } // ii
  
  } // RefitClusters


/////////////////////////////////////////
  void CCHitRefinerAlg::FillWireSignals(
    std::vector<CCHitFinderAlg::CCHit>& allhits)
  {
    // Fill the wire signals vector and initialize the fitted hit
    // signals vector
    gMinuitStruct.WireSignals.clear();
    gMinuitStruct.HitSignals.clear();
    
    // define a zero vector
    std::vector<float> zero(hiTime - loTime + 1);
    
    for(unsigned short wire = loWire; wire <= hiWire; ++wire) {
      // dead wire or no hits. pushback the zero vector
      if(WireHitRange[wire - fFirstWire].first < 0) {
        gMinuitStruct.WireSignals.push_back(zero);
        continue;
      }
      // find a hit on the wire to get the wire object
      unsigned short iht = WireHitRange[wire - fFirstWire].first;
      art::Ptr< recob::Wire> oWire = allhits[iht].Wire;
      // would be nice to just get the wire signal where we need it...
      std::vector<float> signal( oWire->Signal() );
      // temporary RAT vector
      std::vector<float> rat;
      for(unsigned short time = loTime; time <= hiTime; ++time) {
        rat.push_back(signal[time]);
      } // time
      gMinuitStruct.WireSignals.push_back(rat);
    } // wire
    
    // define the HitSignals vector here as well
    for(unsigned short wire = loWire; wire <= hiWire; ++wire) {
      gMinuitStruct.HitSignals.push_back(zero);
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
            <<wire<<" "<<time<<" "<<(int)gMinuitStruct.WireSignals[wIndx][tIndx]
            <<" "<<(int)gMinuitStruct.HitSignals[wIndx][tIndx];
      } // time
    } // wire
  } // PrintSignals


/////////////////////////////////////////
  void CCHitRefinerAlg::FindRATRange(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx)
  {
    // gets the range of wires and times of the RAT surrounding theVtx
    loWire = 9999;
    loTime = 9999;
    hiWire = 0;
    hiTime = 0;
    
    // start with clusters that End at the vertex
    for(unsigned short ii = 0; ii < clEnd.size(); ++ii) {
      unsigned short icl = clEnd[ii];
      for(unsigned short jj = tcl[icl].tclhits.size() - 1; jj > 0; --jj) {
        unsigned short hit = tcl[icl].tclhits[jj];
        if(allhits[hit].WireNum < loWire) loWire = allhits[hit].WireNum;
        if(allhits[hit].LoTime  < loTime) loTime = allhits[hit].LoTime;
        if(allhits[hit].WireNum > hiWire) hiWire = allhits[hit].WireNum;
        if(allhits[hit].HiTime  > hiTime) hiTime = allhits[hit].HiTime;
        // stop looking if the multiplicity = 1 and the hit fit is good
  if(prt) mf::LogVerbatim("ClusterCrawler")
    <<"Chk "<<tcl[icl].ID<<" "<<allhits[hit].WireNum
    <<":"<<(int)allhits[hit].Time<<" mult "<<allhits[hit].numHits;
        if(allhits[hit].numHits == 1) break;
        // stop if it looks squirrely
        if(jj < tcl[icl].tclhits.size() - 10) break;
      } // jj
    } // ii
    
    // now check clusters that Begin at the vertex
    for(unsigned short ii = 0; ii < clBeg.size(); ++ii) {
      unsigned short icl = clBeg[ii];
      for(unsigned short jj = 0; jj < tcl[icl].tclhits.size(); ++jj) {
        unsigned short hit = tcl[icl].tclhits[jj];
        if(allhits[hit].WireNum < loWire) loWire = allhits[hit].WireNum;
        if(allhits[hit].LoTime  < loTime) loTime = allhits[hit].LoTime;
        if(allhits[hit].WireNum > hiWire) hiWire = allhits[hit].WireNum;
        if(allhits[hit].HiTime  > hiTime) hiTime = allhits[hit].HiTime;
        // stop looking if the multiplicity = 1 and the hit fit is good
        if(allhits[hit].numHits == 1) break;
        // stop if it looks squirrely
        if(jj > 10) break;
      } // jj
    } // ii
    
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
        if(allhits[hit].Charge < 0) continue;
        // ignore used hits
        if(allhits[hit].InClus > 0) continue;
        if(allhits[hit].Time > loTime && allhits[hit].Time < hiTime) {
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
        if(allhits[hit].Charge < 0) continue;
        // ignore used hits
        if(allhits[hit].InClus > 0) continue;
        if(allhits[hit].Time > loTime && allhits[hit].Time < hiTime) {
          ++hiWire;
          gotone = true;
          break;
        }
      } // hit
      if(!gotone) break;
    } // wire
    
    // pad the times by some amount
    if(loTime > 10) loTime -= 10;
    if(hiTime < detprop->NumberTimeSamples() - 10) hiTime += 10;

  } // GetRATRange

/////////////////////////////////////////
    void CCHitRefinerAlg::SetClusterBeginEnd(
        std::vector<CCHitFinderAlg::CCHit>& allhits,
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



} // namespace cluster

