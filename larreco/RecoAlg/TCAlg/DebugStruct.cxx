#include "larreco/RecoAlg/TCAlg/DebugStruct.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

#include <boost/algorithm/string/classification.hpp>   // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp>            // Include for boost::split

#include "lardataobj/RecoBase/Hit.h" // for Hit
#include "larreco/RecoAlg/TCAlg/PFPUtils.h"
#include "larreco/RecoAlg/TCAlg/StepUtils.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h" // for tcc


namespace tca {
  DebugStuff debug;

  ////////////////////////////////////////////////
  bool
  DecodeDebugString(std::string strng)
  {
    // try to unpack the string as Cryostat:TPC:Plane:Wire:Tick or something
    // like Slice:<slice index>

    if (strng == "instruct") {
      std::cout << "****** Unrecognized DebugConfig. Here are your options\n";
      std::cout << " 'C:T:P:W:Tick' where C = cryostat, T = TPC, W = wire, Tick (+/-5) to debug "
                   "stepping (DUNE)\n";
      std::cout << " 'P:W:Tick' for single cryostat/TPC detectors (uB, LArIAT, etc)\n";
      std::cout << " 'WorkID <id> <slice index>' where <id> is a tj work ID (< 0) in slice <slice "
                   "index> (default = 0)\n";
      std::cout << " 'CTP <CTP>' to restrict debug output to CTP\n";
      std::cout << " 'Merge <CTP>' to debug trajectory merging\n";
      std::cout << " '2V <CTP>' to debug 2D vertex finding\n";
      std::cout << " '3V' to debug 3D vertex finding\n";
      std::cout << " 'VxMerge' to debug 2D vertex merging\n";
      std::cout << " 'JunkVx' to debug 2D junk vertex finder\n";
      std::cout << " 'PFP' to debug 3D matching and PFParticles\n";
      std::cout << " 'MVI <MVI> <MVI Iteration>' for detailed debugging of one PFP MatchVecIndex\n";
      std::cout << " 'DeltaRay' to debug delta ray tagging\n";
      std::cout << " 'Muon' to debug muon tagging\n";
      std::cout << " '2S <CTP>' to debug a 2D shower in CTP\n";
      std::cout << " 'Reco TPC <TPC>' to only reconstruct hits in the specified TPC\n";
      std::cout << " 'Reco Slice <ID>' to reconstruct all sub-slices in the recob::Slice with the "
                   "specified ID\n";
      std::cout << " 'SubSlice <sub-slice index>' where <slice index> restricts output to the "
                   "specified sub-slice index\n";
      std::cout << " 'Stitch' to debug PFParticle stitching between TPCs\n";
      std::cout << " 'Sum' or 'Summary' to print a debug summary report\n";
      std::cout << " 'Dump <WorkID>' or 'Dump <UID>' to print all TPs in the trajectory to "
                   "tcdump<UID>.csv\n";
      std::cout << " Note: Algs with debug printing include HamVx, HamVx2, SplitTjCVx, Comp3DVx, "
                   "Comp3DVxIG, VtxHitsSwap\n";
      std::cout << " Set SkipAlgs: [\"bogusText\"] to print a list of algorithm names\n";
      return false;
    } // instruct

    // handle the simple cases that don't need decoding
    for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
      if(strng.find(AlgBitNames[ib]) == 0) {
        tcc.dbgAlg[ib] = true;
        return true;
      }
    } // ib
    if (strng.find("3V") == 0) {
      tcc.dbg3V = true;
      tcc.modes[kModeDebug] = true;
      return true;
    }
    if (strng.find("3S") != std::string::npos) {
      tcc.dbg3S = true;
      tcc.modes[kModeDebug] = true;
      return true;
    }
    if (strng.find("Stitch") != std::string::npos) {
      tcc.dbgStitch = true;
      tcc.modes[kModeDebug] = true;
      return true;
    }
    if (strng.find("Sum") != std::string::npos) {
      tcc.dbgSummary = true;
      tcc.modes[kModeDebug] = true;
      return true;
    }

    std::vector<std::string> words;
    boost::split(words, strng, boost::is_any_of(" :"), boost::token_compress_on);
    if (words.size() == 5) {
      // configure for DUNE
      debug.Cryostat = std::stoi(words[0]);
      debug.TPC = std::stoi(words[1]);
      debug.Plane = std::stoi(words[2]);
      debug.Wire = std::stoi(words[3]);
      debug.Tick = std::stoi(words[4]);
      debug.CTP = EncodeCTP(debug.Cryostat, debug.TPC, debug.Plane);
      tcc.modes[kModeDebug] = true;
      tcc.dbgStp = true;
      // also dump this tj
      tcc.dbgDump = true;
      return true;
    } // nums.size() == 5
    if (words[0] == "PFP" || words[0] == "MVI") {
      tcc.dbgPFP = true;
      tcc.modes[kModeDebug] = true;
      // Use debug.Hit to identify the matchVec index
      if (words.size() > 2) {
        debug.MVI = std::stoi(words[1]);
        if (words.size() == 3) debug.MVI_Iter = std::stoi(words[2]);
      }
      return true;
    } // PFP
    if (words.size() == 2 && words[0] == "CTP") {
      debug.CTP = std::stoi(words[1]);
      auto plnID = DecodeCTP(debug.CTP);
      debug.Cryostat = plnID.Cryostat;
      debug.TPC = plnID.TPC;
      debug.Plane = plnID.Plane;
      return true;
    } // CTP
    if (words.size() == 2 && words[0] == "Dump") {
      debug.WorkID = std::stoi(words[1]);
      debug.Slice = 0;
      tcc.modes[kModeDebug] = true;
      tcc.dbgDump = true;
      return true;
    }
    if (words.size() > 1 && words[0] == "WorkID") {
      debug.WorkID = std::stoi(words[1]);
      if (debug.WorkID >= 0) return false;
      // default to sub-slice index 0
      debug.Slice = 0;
      if (words.size() > 2) debug.Slice = std::stoi(words[2]);
      tcc.modes[kModeDebug] = true;
      // dbgStp is set true after debug.WorkID is found
      tcc.dbgStp = false;
      return true;
    } // words.size() == 3 && words[0] == "WorkID"
    if (words.size() == 3 && words[0] == "Reco" && words[1] == "TPC") {
      tcc.recoTPC = std::stoi(words[2]);
      tcc.modes[kModeDebug] = true;
      std::cout << "Reconstructing only in TPC " << tcc.recoTPC << "\n";
      return true;
    }
    if(words.size() == 3 && words[0] == "Reco" && words[1] == "Slice") {
      tcc.recoSlice = std::stoi(words[2]);
      std::cout<<"Reconstructing Slice "<<tcc.recoSlice<<"\n";
      return true;
    }
    if (words.size() == 3) {
      // configure for uB, LArIAT, etc
      debug.Cryostat = 0;
      debug.TPC = 0;
      debug.Plane = std::stoi(words[0]);
      debug.Wire = std::stoi(words[1]);
      debug.Tick = std::stoi(words[2]);
      debug.CTP = EncodeCTP(debug.Cryostat, debug.TPC, debug.Plane);
      tcc.modes[kModeDebug] = true;
      tcc.dbgStp = true;
      return true;
    }
    if (words.size() == 2 && words[0] == "Merge") {
      debug.CTP = std::stoi(words[1]);
      tcc.dbgMrg = true;
      tcc.modes[kModeDebug] = true;
      return true;
    }
    if (words.size() == 2 && words[0] == "2V") {
      debug.CTP = std::stoi(words[1]);
      auto plnID = DecodeCTP(debug.CTP);
      debug.Cryostat = plnID.Cryostat;
      debug.TPC = plnID.TPC;
      debug.Plane = plnID.Plane;
      tcc.dbg2V = true;
      tcc.modes[kModeDebug] = true;
      return true;
    }
    if (words.size() == 2 && words[0] == "2S") {
      debug.CTP = std::stoi(words[1]);
      tcc.dbg2S = true;
      tcc.modes[kModeDebug] = true;
      return true;
    }
    // Slice could apply to several debug options.
    if (words.size() == 2 && words[0] == "SubSlice") {
      debug.Slice = std::stoi(words[1]);
      return true;
    }
    return false;
  } // DecodeDebugString

  void
  DumpTj()
  {
    // Dump all of the points in a trajectory to the output in a form that can
    // be imported by another application, e.g. Excel
    // Search for the trajectory with the specified WorkID or Unique ID

    for (auto& slc : slices) {
      for (auto& tj : slc.tjs) {
        if (tj.WorkID != debug.WorkID && tj.UID != debug.WorkID) continue;
        // print a header
        std::ofstream outfile;
        std::string fname = "tcdump" + std::to_string(tj.UID) + ".csv";
        outfile.open(fname, std::ios::out | std::ios::trunc);
        outfile << "Dump trajectory T" << tj.UID << " WorkID " << tj.WorkID;
        outfile << " ChgRMS " << std::setprecision(2) << tj.ChgRMS;
        outfile << "\n";
        outfile << "Wire, Chg T" << tj.UID
                << ", totChg, Tick, Delta, NTPsFit, Ang, ChiDOF, KinkSig, HitPosErr\n";
        for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          auto& tp = tj.Pts[ipt];
          outfile << std::fixed;
          outfile << std::setprecision(0) << std::nearbyint(tp.Pos[0]);
          outfile << "," << (int)tp.Chg;
          // total charge near the TP
          float totChg = 0;
          for (auto iht : tp.Hits) {
            auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
            totChg += hit.Integral();
          }
          outfile << "," << (int)totChg;
          outfile << "," << std::setprecision(0) << std::nearbyint(tp.Pos[1] / tcc.unitsPerTick);
          outfile << "," << std::setprecision(2) << tp.Delta;
          outfile << "," << tp.NTPsFit;
          outfile << "," << std::setprecision(3) << tp.Ang;
          outfile << "," << std::setprecision(2) << tp.FitChi;
          outfile << "," << std::setprecision(2) << tp.KinkSig;
          outfile << "," << std::setprecision(2) << sqrt(tp.HitPosErr2);
          outfile << "\n";
        } // ipt
        outfile.close();
        std::cout<<"Points on T"<<tj.UID<<" dumped to "<<fname<<"\n";
//        tcc.dbgDump = false;
//        return;
      } // tj
    }   // slc

  } // DumpTj

  ////////////////////////////////////////////////
  void
  PrintDebugMode()
  {
    // print the debug mode configuration to the screen
    std::cout << "*** TrajCluster debug mode configuration in";
    std::cout << " CTP=";
    if (debug.CTP == UINT_MAX) { std::cout << "NA"; }
    else {
      std::cout << debug.CTP;
      auto plnID = DecodeCTP(debug.CTP);
      debug.Cryostat = plnID.Cryostat;
      debug.TPC = plnID.TPC;
      debug.Plane = plnID.Plane;
    }
    std::cout << " Cryostat=" << debug.Cryostat;
    std::cout << " TPC=" << debug.TPC;
    std::cout << " Plane=" << debug.Plane;
    std::cout << " Wire=" << debug.Wire;
    std::cout << " Tick=" << debug.Tick;
    std::cout << " Hit=";
    if (debug.Hit == UINT_MAX) { std::cout << "NA"; }
    else {
      std::cout << debug.Hit;
    }
    std::cout << " WorkID=";
    if (debug.WorkID == 0) { std::cout << "NA"; }
    else {
      std::cout << debug.WorkID;
    }
    std::cout << " Slice=";
    if (debug.Slice == -1) { std::cout << "All"; }
    else {
      std::cout << debug.Slice;
    }
    std::cout << "\n";
    std::cout << "*** tcc.dbg modes:";
    if (tcc.dbgSlc) std::cout << " dbgSlc";
    if (tcc.dbgStp) std::cout << " dbgStp";
    if (tcc.dbgMrg) std::cout << " dbgMrg";
    if (tcc.dbg2V) std::cout << " dbg2V";
    if (tcc.dbg2S) std::cout << " dbg2S";
    if (tcc.dbg3V) std::cout << " dbg3V";
    if (tcc.dbgPFP) std::cout << " dbgPFP";
    if (tcc.dbgStitch) std::cout << " dbgStitch";
    if (tcc.dbgSummary) std::cout << " dbgSummary";
    if (tcc.dbgDump) std::cout << " dbgDump";
    for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
      if(tcc.dbgAlg[ib]) std::cout<<" "<<AlgBitNames[ib];
    } // ib
    std::cout << "\n";
    std::cout << "*** Using algs:";
    unsigned short cnt = 0;
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
      if (tcc.useAlg[ib] && ib != kKilled) {
        ++cnt;
        if (cnt % 10 == 0) std::cout << "\n   ";
        std::cout << " " << AlgBitNames[ib];
      }
    }
    std::cout << "\n";
    std::cout << "*** Skipping algs:";
    cnt = 0;
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
      if (!tcc.useAlg[ib] && ib != kKilled) {
        ++cnt;
        if (cnt % 10 == 0) std::cout << "\n   ";
        std::cout << " " << AlgBitNames[ib];
      }
    }
    std::cout << "\n";
    std::cout << "TrajCluster configuration modes:";
    if(tcc.modes[kModeStepPos]) {
      std::cout << " StepPos"; 
    } else {
       std::cout << " StepNeg";
    }
    if(tcc.modes[kModeTestBeam]) {
      std::cout<<", TestBeam";
    } else if(tcc.modes[kModeLEPhysics]) {
      std::cout<<", LEPhysics";
    } else {
      std::cout<<", Neutrino (default)";
    }
    std::cout << "\n";
  } // PrintDebugMode

   void PrintAll(detinfo::DetectorPropertiesData const& detProp, std::string someText)
  {
    // print everything in all slices
    bool prt3V = false;
    bool prt2V = false;
    bool prtT = false;
    bool prtP = false;
    for (size_t isl = 0; isl < slices.size(); ++isl) {
      if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
      auto& slc = slices[isl];
      if (!slc.vtx3s.empty()) prt3V = true;
      if (!slc.vtxs.empty()) prt2V = true;
      if (!slc.tjs.empty()) prtT = true;
      if (!slc.pfps.empty()) prtP = true;
    } // slc
    mf::LogVerbatim myprt("TC");
    myprt << "Debug report from caller " << someText << "\n";
    someText = "";
    myprt << " 'prodID' = <sliceID>:<subSliceIndex>:<productID>/<productUID>\n";
    if (prtP) {
      bool printHeader = true;
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.pfps.empty()) continue;
        for (auto& pfp : slc.pfps)
          PrintP(someText, myprt, pfp, printHeader);
      } // slc
    }   // prtP
    if (prt3V) {
      bool printHeader = true;
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.vtx3s.empty()) continue;
        for (auto& vx3 : slc.vtx3s)
          Print3V(detProp, someText, myprt, vx3, printHeader);
      } // slc
    }   // prt3V
    if (prt2V) {
      bool printHeader = true;
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.vtxs.empty()) continue;
        for (auto& vx2 : slc.vtxs)
          Print2V(someText, myprt, vx2, printHeader);
      } // slc
    }   // prt2V
    if (prtT) {
      bool printHeader = true;
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.tjs.empty()) continue;
        for (auto& tj : slc.tjs)
          PrintT(someText, myprt, tj, printHeader);
      } // slc
    }   // prtT
  }     // PrintAll

  ////////////////////////////////////////////////
  void
  PrintP(std::string someText, mf::LogVerbatim& myprt, PFPStruct& pfp, bool& printHeader)
  {
    if (pfp.ID <= 0) return;
    if (printHeader) {
      myprt << "************ PFParticles ************\n";
      myprt << "     prodID    sVx  _____sPos____ CS _______sDir______ ____sdEdx_____    eVx  "
               "_____ePos____ CS ____edEdx_____  MVI MCSMom  Len nTP3 nSec SLk? PDG  Par \n";
      printHeader = false;
    } // printHeader
    auto sIndx = GetSliceIndex("P", pfp.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    std::string str =
      std::to_string(slc.ID) + ":" + std::to_string(sIndx.first) + ":" + std::to_string(pfp.ID);
    str += "/" + std::to_string(pfp.UID);
    myprt << std::setw(12) << str;
    // start and end stuff
    for (unsigned short end = 0; end < 2; ++end) {
      str = "--";
      if (pfp.Vx3ID[end] > 0) str = "3V" + std::to_string(slc.vtx3s[pfp.Vx3ID[end] - 1].UID);
      myprt << std::setw(6) << str;
      myprt << std::fixed << std::right << std::setprecision(0);
      auto tp3d = EndTP3D(pfp, end);
      auto pos = tp3d.Pos;
      myprt << std::setw(5) << pos[0];
      myprt << std::setw(5) << pos[1];
      myprt << std::setw(5) << pos[2];
      // print character for Outside or Inside the FV
      if (InsideFV(slc, pfp, end)) { myprt << "  I"; }
      else {
        myprt << "  O";
      }
      // only print the starting direction
      if (end == 0) {
        myprt << std::fixed << std::right << std::setprecision(2);
        auto dir = tp3d.Dir;
        myprt << std::setw(6) << dir[0];
        myprt << std::setw(6) << dir[1];
        myprt << std::setw(6) << dir[2];
      } // end == 0
      for (auto& dedx : pfp.dEdx[end]) {
        if (dedx < 50) { myprt << std::setw(5) << std::setprecision(1) << dedx; }
        else {
          myprt << std::setw(5) << std::setprecision(0) << dedx;
        }
      } // dedx
      if (pfp.dEdx[end].size() < 3) {
        for (size_t i = 0; i < 3 - pfp.dEdx[end].size(); ++i) {
          myprt << std::setw(6) << ' ';
        }
      }
    } // startend
    myprt << std::setw(6) << pfp.MVI;
    // global stuff
    myprt << std::setw(7) << MCSMom(slc, pfp.TjIDs);
    float length = Length(pfp);
    if (length < 100) { myprt << std::setw(5) << std::setprecision(1) << length; }
    else {
      myprt << std::setw(5) << std::setprecision(0) << length;
    }
    myprt << std::setw(5) << pfp.TP3Ds.size();
    myprt << std::setw(5) << pfp.SectionFits.size();
    myprt << std::setw(5) << IsShowerLike(slc, pfp.TjIDs);
    myprt << std::setw(5) << pfp.PDGCode;
    myprt << std::setw(4) << pfp.ParentUID;
    if (!pfp.TjIDs.empty()) {
      if (pfp.TjUIDs.empty()) {
        // print Tjs in one TPC
        for (auto tjid : pfp.TjIDs)
          myprt << " TU" << slc.tjs[tjid - 1].UID;
      }
      else {
        // print Tjs in all TPCs (if this is called after FinishEvent)
        for (auto tjuid : pfp.TjUIDs)
          myprt << " TU" << tjuid;
      }
    } // TjIDs exist
    if (!pfp.DtrUIDs.empty()) {
      myprt << " dtrs";
      for (auto dtruid : pfp.DtrUIDs)
        myprt << " PU" << dtruid;
    } // dtr ids exist
    myprt << "\n";
  } // PrintP

  ////////////////////////////////////////////////
  void
  Print3V(detinfo::DetectorPropertiesData const& detProp,
          std::string someText,
          mf::LogVerbatim& myprt,
          Vtx3Store& vx3,
          bool& printHeader)
  {
    // print a 3D vertex on one line
    if (vx3.ID <= 0) return;
    auto sIndx = GetSliceIndex("3V", vx3.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    if (printHeader) {
      myprt
        << "****** 3D vertices ******************************************__2DVtx_UID__*******\n";
      myprt << "     prodID    Cstat TPC     X       Y       Z    pln0   pln1   pln2 Wire score "
               "Prim? Nu? nTru";
      myprt << " ___________2D_Pos____________ _____Tj UIDs________\n";
      printHeader = false;
    }
    std::string str = "3V" + std::to_string(vx3.ID) + "/3VU" + std::to_string(vx3.UID);
    myprt << std::right << std::setw(12) << std::fixed << str;
    myprt << std::setprecision(0);
    myprt << std::right << std::setw(7) << vx3.TPCID.Cryostat;
    myprt << std::right << std::setw(5) << vx3.TPCID.TPC;
    myprt << std::right << std::setw(8) << vx3.X;
    myprt << std::right << std::setw(8) << vx3.Y;
    myprt << std::right << std::setw(8) << vx3.Z;
    for (auto vx2id : vx3.Vx2ID) {
      if (vx2id > 0) {
        str = "2VU" + std::to_string(slc.vtxs[vx2id - 1].UID);
        myprt << std::right << std::setw(7) << str;
      }
      else {
        myprt << "   --";
      }
    } // vx2id
    myprt << std::right << std::setw(5) << vx3.Wire;
    unsigned short nTruMatch = 0;
    for (unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
      if (vx3.Vx2ID[ipl] == 0) continue;
      unsigned short iv2 = vx3.Vx2ID[ipl] - 1;
      if (slc.vtxs[iv2].Stat[kVxTruMatch]) ++nTruMatch;
    } // ipl
    myprt << std::right << std::setw(6) << std::setprecision(1) << vx3.Score;
    myprt << std::setw(6) << vx3.Primary;
    myprt << std::setw(4) << vx3.Neutrino;
    myprt << std::right << std::setw(5) << nTruMatch;
    Point2_t pos;
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      PosInPlane(detProp, slc, vx3, plane, pos);
      myprt << " " << PrintPos(slc, pos);
    } // plane
    if (vx3.Wire == -2) {
      // find the Tjs that are attached to it
      for (unsigned short end = 0; end < 2; ++end) {
        for (auto& pfp : slc.pfps) {
          if (pfp.Vx3ID[end] == vx3.ID) {
            for (auto tjID : pfp.TjIDs) {
              auto& tj = slc.tjs[tjID - 1];
              myprt << " T" << tj.UID;
            } // tjID
          }   // pfp.Vx3ID[0] == vx3.ID
        }     // pfp
      }       // end
    }
    else {
      auto vxtjs = GetAssns(slc, "3V", vx3.ID, "T");
      for (auto tjid : vxtjs) {
        auto& tj = slc.tjs[tjid - 1];
        myprt << " TU" << tj.UID;
      }
    } // vx3.Wire != -2
    myprt << "\n";
  } // Print3V

  ////////////////////////////////////////////////
  void
  Print2V(std::string someText, mf::LogVerbatim& myprt, VtxStore& vx2, bool& printHeader)
  {
    // print a 2D vertex on one line
    if (vx2.ID <= 0) return;
    if (debug.CTP != UINT_MAX && vx2.CTP != debug.CTP) return;
    auto sIndx = GetSliceIndex("2V", vx2.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    if (printHeader) {
      myprt << someText << "************ 2D vertices ************\n";
      myprt << someText << "     prodID    CTP    wire  err   tick   err  ChiDOF  NTj Pass  Topo"
               " ChgFrac Score  v3D Tj UIDs\n";
      printHeader = false;
    }
    std::string str = "2V" + std::to_string(vx2.ID) + "/2VU" + std::to_string(vx2.UID);
    myprt << someText;
    myprt << std::right << std::setw(12) << std::fixed << str;
    myprt << std::right << std::setw(6) << vx2.CTP;
    myprt << std::right << std::setw(8) << std::setprecision(0) << std::nearbyint(vx2.Pos[0]);
    myprt << std::right << std::setw(5) << std::setprecision(1) << vx2.PosErr[0];
    myprt << std::right << std::setw(8) << std::setprecision(0)
          << std::nearbyint(vx2.Pos[1] / tcc.unitsPerTick);
    myprt << std::right << std::setw(5) << std::setprecision(1) << vx2.PosErr[1] / tcc.unitsPerTick;
    myprt << std::right << std::setw(7) << vx2.ChiDOF;
    myprt << std::right << std::setw(5) << vx2.NTraj;
    myprt << std::right << std::setw(5) << vx2.Pass;
    myprt << std::right << std::setw(6) << vx2.Topo;
    myprt << std::right << std::setw(9) << std::setprecision(2) << vx2.TjChgFrac;
    myprt << std::right << std::setw(6) << std::setprecision(1) << vx2.Score;
    int v3id = 0;
    if (vx2.Vx3ID > 0) v3id = slc.vtx3s[vx2.Vx3ID - 1].UID;
    myprt << std::right << std::setw(5) << v3id;
    myprt << "    ";
    // display the traj IDs
    for (unsigned short ii = 0; ii < slc.tjs.size(); ++ii) {
      auto const& tj = slc.tjs[ii];
      if (tj.AlgMod[kKilled]) continue;
      for (unsigned short end = 0; end < 2; ++end) {
        if (tj.VtxID[end] != (short)vx2.ID) continue;
        std::string tid = " TU" + std::to_string(tj.UID) + "_" + std::to_string(end);
        myprt << std::right << std::setw(6) << tid;
      } // end
    }   // ii
    myprt << " Stat:";
    // Special flags. Ignore the first flag bit (0 = kVxTrjTried) which is done for every vertex
    for (unsigned short ib = 1; ib < VtxBitNames.size(); ++ib)
      if (vx2.Stat[ib]) myprt << " " << VtxBitNames[ib];
    myprt << "\n";
  } // Print2V

  ////////////////////////////////////////////////
  void
  PrintT(std::string someText, mf::LogVerbatim& myprt, Trajectory& tj, bool& printHeader)
  {
    // print a 2D vertex on one line
    if(tj.ID <= 0) return;
    if(debug.CTP != UINT_MAX && tj.CTP != debug.CTP) return;
    if(printHeader) {
      myprt<<"************ Trajectories ************\n";
      myprt<<"Tj AngleCode-EndFlag decoder (EF): <AngleCode> + <end flag>";
      myprt<<" (B=Bragg Peak, V=Vertex, A=AngleKink, C=ChargeKink, T=Trajectory, S=StartEnd)\n";
      myprt<<"     prodID    CTP  Pts     W:T      Ang EF   AveQ     W:T      Ang EF   ";
      myprt<<"AveQ Chg(k) chgRMS  Mom __Vtx__  PDG WorkID \n";
      printHeader = false;
    }
    auto sIndx = GetSliceIndex("T", tj.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    std::string str = "T";
    if(tj.AlgMod[kKilled]) str = "x";
    str = str + std::to_string(tj.ID) + "/TU" + std::to_string(tj.UID);
    myprt << std::fixed << std::setw(12) << str;
    myprt << std::setw(6) << tj.CTP;
    myprt << std::setw(5) << tj.EndPt[1] - tj.EndPt[0] + 1;
    unsigned short endPt0 = tj.EndPt[0];
    auto& tp0 = tj.Pts[endPt0];
    int itick = tp0.Pos[1]/tcc.unitsPerTick;
    if(itick < 0) itick = 0;
    myprt<<std::setw(6)<<(int)(tp0.Pos[0]+0.5)<<":"<<itick; // W:T
    if(itick < 10) { myprt<<" "; }
    if(itick < 100) { myprt<<" "; }
    if(itick < 1000) { myprt<<" "; }
    myprt<<std::setw(6)<<std::setprecision(2)<<tp0.Ang;
    myprt<<std::setw(2)<<tp0.AngleCode;
    myprt<<std::setw(3)<<std::left<<PackEndFlags(tj, 0)<<std::right;
    myprt<<std::setw(5)<<(int)tp0.AveChg;
    unsigned short endPt1 = tj.EndPt[1];
    auto& tp1 = tj.Pts[endPt1];
    itick = tp1.Pos[1]/tcc.unitsPerTick;
    myprt<<std::setw(6)<<(int)(tp1.Pos[0]+0.5)<<":"<<itick; // W:T
    if(itick < 10) { myprt<<" "; }
    if(itick < 100) { myprt<<" "; }
    if(itick < 1000) { myprt<<" "; }
    myprt<<std::setw(6)<<std::setprecision(2)<<tp1.Ang;
    myprt<<std::setw(2)<<tp1.AngleCode;
    myprt<<std::setw(3)<<std::left<<PackEndFlags(tj, 1)<<std::right;
    myprt<<std::setw(5)<<(int)tp1.AveChg;
    myprt<<std::setw(7)<<std::setprecision(1)<<tj.TotChg/1000;
    myprt<<std::setw(7)<<std::setprecision(2)<<tj.ChgRMS;
    myprt<<std::setw(5)<<tj.MCSMom;
    int vxid = 0;
    if (tj.VtxID[0] > 0) vxid = slc.vtxs[tj.VtxID[0] - 1].UID;
    myprt << std::setw(4) << vxid;
    vxid = 0;
    if (tj.VtxID[1] > 0) vxid = slc.vtxs[tj.VtxID[1] - 1].UID;
    myprt << std::setw(4) << vxid;
    myprt << std::setw(5) << tj.PDGCode;
    myprt << std::setw(7) << tj.WorkID;
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
      if (tj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
    for (unsigned short ib = 0; ib < StrategyBitNames.size(); ++ib)
      if (tj.Strategy[ib]) myprt << " " << StrategyBitNames[ib];
    myprt << "\n";
  } // PrintT

  ////////////////////////////////////////////////
  std::string
  PackEndFlags(const Trajectory& tj, unsigned short end)
  {
    // Returns a short string of trajectory EndFlags
    std::string endCode;
    if(end > 1) return endCode;
    if(end == tj.StartEnd) endCode = "S";
    if(tj.EndFlag[1][kEndBragg]) endCode = endCode + "B";
    if(tj.EndFlag[1][kEndKink]) endCode = endCode + "K";
    if(tj.EndFlag[1][kHitsAfterEnd]) endCode = endCode + "H";
    return endCode;
  } // PackEndFlags

  ////////////////////////////////////////////////
  void
  PrintAllTraj(detinfo::DetectorPropertiesData const& detProp,
               std::string someText,
               TCSlice& slc,
               unsigned short itj,
               unsigned short ipt,
               bool prtVtx)
  {

    mf::LogVerbatim myprt("TC");

    if (prtVtx) {
      if (!slc.vtx3s.empty()) {
        // print out 3D vertices
        myprt
          << someText
          << "****** 3D vertices ******************************************__2DVtx_ID__*******\n";
        myprt << someText
              << "  Vtx  Cstat  TPC     X       Y       Z    XEr  YEr  ZEr pln0 pln1 pln2 Wire "
                 "score Prim? Nu? nTru";
        myprt << " ___________2D_Pos____________ _____Tjs________\n";
        for (unsigned short iv = 0; iv < slc.vtx3s.size(); ++iv) {
          if (slc.vtx3s[iv].ID == 0) continue;
          const Vtx3Store& vx3 = slc.vtx3s[iv];
          myprt << someText;
          std::string vid = "3v" + std::to_string(vx3.ID);
          myprt << std::right << std::setw(5) << std::fixed << vid;
          myprt << std::setprecision(1);
          myprt << std::right << std::setw(7) << vx3.TPCID.Cryostat;
          myprt << std::right << std::setw(5) << vx3.TPCID.TPC;
          myprt << std::right << std::setw(8) << vx3.X;
          myprt << std::right << std::setw(8) << vx3.Y;
          myprt << std::right << std::setw(8) << vx3.Z;
          myprt << std::right << std::setw(5) << vx3.XErr;
          myprt << std::right << std::setw(5) << vx3.YErr;
          myprt << std::right << std::setw(5) << vx3.ZErr;
          myprt << std::right << std::setw(5) << vx3.Vx2ID[0];
          myprt << std::right << std::setw(5) << vx3.Vx2ID[1];
          myprt << std::right << std::setw(5) << vx3.Vx2ID[2];
          myprt << std::right << std::setw(5) << vx3.Wire;
          unsigned short nTruMatch = 0;
          for (unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
            if (vx3.Vx2ID[ipl] == 0) continue;
            unsigned short iv2 = vx3.Vx2ID[ipl] - 1;
            if (slc.vtxs[iv2].Stat[kVxTruMatch]) ++nTruMatch;
          } // ipl
          myprt << std::right << std::setw(6) << std::setprecision(1) << vx3.Score;
          myprt << std::setw(6) << vx3.Primary;
          myprt << std::setw(4) << vx3.Neutrino;
          myprt << std::right << std::setw(5) << nTruMatch;
          Point2_t pos;
          for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
            PosInPlane(detProp, slc, vx3, plane, pos);
            myprt << " " << PrintPos(slc, pos);
          } // plane
          if (vx3.Wire == -2) {
            // find the Tjs that are attached to it
            for (auto& pfp : slc.pfps) {
              if (pfp.Vx3ID[0] == slc.vtx3s[iv].ID) {
                for (auto& tjID : pfp.TjIDs)
                  myprt << " t" << tjID;
              }
              if (pfp.Vx3ID[1] == slc.vtx3s[iv].ID) {
                for (auto& tjID : pfp.TjIDs)
                  myprt << " t" << tjID;
              }
            } // ipfp
          }
          else {
            auto vxtjs = GetAssns(slc, "3V", vx3.ID, "T");
            for (auto tjid : vxtjs)
              myprt << " t" << tjid;
          }
          myprt << "\n";
        }
      } // slc.vtx3s.size
      if (!slc.vtxs.empty()) {
        bool foundOne = false;
        for (unsigned short iv = 0; iv < slc.vtxs.size(); ++iv) {
          auto& vx2 = slc.vtxs[iv];
          if (debug.Plane < 3 && debug.Plane != (int)DecodeCTP(vx2.CTP).Plane) continue;
          if (vx2.NTraj == 0) continue;
          foundOne = true;
        } // iv
        if (foundOne) {
          // print out 2D vertices
          myprt << someText << "************ 2D vertices ************\n";
          myprt << someText
                << " ID   CTP    wire  err   tick   err  ChiDOF  NTj Pass  Topo ChgFrac Score  v3D "
                   "TjIDs\n";
          for (auto& vx2 : slc.vtxs) {
            if (vx2.ID == 0) continue;
            if (debug.Plane < 3 && debug.Plane != (int)DecodeCTP(vx2.CTP).Plane) continue;
            myprt << someText;
            std::string vid = "2v" + std::to_string(vx2.ID);
            myprt << std::right << std::setw(5) << std::fixed << vid;
            myprt << std::right << std::setw(6) << vx2.CTP;
            myprt << std::right << std::setw(8) << std::setprecision(0)
                  << std::nearbyint(vx2.Pos[0]);
            myprt << std::right << std::setw(5) << std::setprecision(1) << vx2.PosErr[0];
            myprt << std::right << std::setw(8) << std::setprecision(0)
                  << std::nearbyint(vx2.Pos[1] / tcc.unitsPerTick);
            myprt << std::right << std::setw(5) << std::setprecision(1)
                  << vx2.PosErr[1] / tcc.unitsPerTick;
            myprt << std::right << std::setw(7) << vx2.ChiDOF;
            myprt << std::right << std::setw(5) << vx2.NTraj;
            myprt << std::right << std::setw(5) << vx2.Pass;
            myprt << std::right << std::setw(6) << vx2.Topo;
            myprt << std::right << std::setw(9) << std::setprecision(2) << vx2.TjChgFrac;
            myprt << std::right << std::setw(6) << std::setprecision(1) << vx2.Score;
            myprt << std::right << std::setw(5) << vx2.Vx3ID;
            myprt << "    ";
            // display the traj IDs
            for (unsigned short ii = 0; ii < slc.tjs.size(); ++ii) {
              auto const& aTj = slc.tjs[ii];
              if (debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aTj.CTP).Plane) continue;
              if (aTj.AlgMod[kKilled]) continue;
              for (unsigned short end = 0; end < 2; ++end) {
                if (aTj.VtxID[end] != (short)vx2.ID) continue;
                std::string tid = " t" + std::to_string(aTj.ID) + "_" + std::to_string(end);
                myprt << std::right << std::setw(6) << tid;
              } // end
            }   // ii
            // Special flags. Ignore the first flag bit (0 = kVxTrjTried) which is done for every vertex
            for (unsigned short ib = 1; ib < VtxBitNames.size(); ++ib)
              if (vx2.Stat[ib]) myprt << " " << VtxBitNames[ib];
            myprt << "\n";
          } // iv
        }
      } // slc.vtxs.size
    }

    if (slc.tjs.empty()) {
      mf::LogVerbatim("TC") << someText << " No allTraj trajectories to print";
      return;
    }

    // Print all trajectories in slc.tjs if itj == USHRT_MAX
    // Print a single traj (itj) and a single TP (ipt) or all TPs (USHRT_MAX)
    if (itj == USHRT_MAX) {
      // Print summary trajectory information
      myprt << "Tj AngleCode-EndFlag (EF) decoder: <AngleCode> + <reason for stopping>";
      myprt << " (B=Bragg Peak, V=Vertex, A=AngleKink, C=ChargeKink, T=Trajectory)\n";
      std::vector<unsigned int> tmp;
      myprt << someText
            << "   UID   CTP Pass  Pts     W:T      Ang EF AveQ     W:T      Ang EF AveQ Chg(k) "
               "chgRMS  Mom SDr __Vtx__  PDG  Par Pri NuPar   WorkID \n";
      for (unsigned short ii = 0; ii < slc.tjs.size(); ++ii) {
        auto& aTj = slc.tjs[ii];
        if (debug.CTP != UINT_MAX && aTj.CTP != debug.CTP) continue;
        myprt << someText << " ";
        std::string tid;
        if (aTj.AlgMod[kKilled]) { tid = "k" + std::to_string(aTj.UID); }
        else {
          tid = "t" + std::to_string(aTj.UID);
        }
        myprt << std::fixed << std::setw(5) << tid;
        myprt << std::setw(6) << aTj.CTP;
        myprt << std::setw(5) << aTj.Pass;
        myprt << std::setw(5) << aTj.EndPt[1] - aTj.EndPt[0] + 1;
        unsigned short endPt0 = aTj.EndPt[0];
        auto& tp0 = aTj.Pts[endPt0];
        int itick = tp0.Pos[1] / tcc.unitsPerTick;
        if (itick < 0) itick = 0;
        myprt << std::setw(6) << (int)(tp0.Pos[0] + 0.5) << ":" << itick; // W:T
        if (itick < 10) { myprt << " "; }
        if (itick < 100) { myprt << " "; }
        if (itick < 1000) { myprt << " "; }
        myprt << std::setw(6) << std::setprecision(2) << tp0.Ang;
        myprt << std::setw(2) << tp0.AngleCode;
        myprt << std::setw(4) << PackEndFlags(aTj, 0);
        myprt << std::setw(5) << (int)tp0.AveChg;
        unsigned short endPt1 = aTj.EndPt[1];
        auto& tp1 = aTj.Pts[endPt1];
        itick = tp1.Pos[1] / tcc.unitsPerTick;
        myprt << std::setw(6) << (int)(tp1.Pos[0] + 0.5) << ":" << itick; // W:T
        if (itick < 10) { myprt << " "; }
        if (itick < 100) { myprt << " "; }
        if (itick < 1000) { myprt << " "; }
        myprt << std::setw(6) << std::setprecision(2) << tp1.Ang;
        myprt << std::setw(2) << tp1.AngleCode;
        myprt << std::setw(4) << PackEndFlags(aTj, 1);
        myprt << std::setw(5) << (int)tp1.AveChg;
        myprt << std::setw(7) << std::setprecision(1) << aTj.TotChg / 1000;
        myprt << std::setw(7) << std::setprecision(2) << aTj.ChgRMS;
        myprt << std::setw(5) << aTj.MCSMom;
        myprt << std::setw(4) << aTj.StepDir;
        myprt << std::setw(4) << aTj.VtxID[0];
        myprt << std::setw(4) << aTj.VtxID[1];
        myprt << std::setw(5) << aTj.PDGCode;
        myprt << std::setw(5) << aTj.ParentID;
        myprt << std::setw(5) << PrimaryID(slc, aTj);
        myprt << std::setw(6) << NeutrinoPrimaryTjID(slc, aTj);
        myprt << std::setw(7) << aTj.WorkID;
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
          if (aTj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
        myprt << "\n";
      } // ii
      return;
    } // itj > slc.tjs.size()-1

    if (itj > slc.tjs.size() - 1) return;

    auto const& aTj = slc.tjs[itj];

    mf::LogVerbatim("TC") << "Print slc.tjs[" << itj << "] Vtx[0] " << aTj.VtxID[0] << " Vtx[1] "
                          << aTj.VtxID[1];
    myprt << "AlgBits";
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
      if (aTj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
    myprt << "\n";

    PrintTPHeader(someText);
    if (ipt == USHRT_MAX) {
      // print all points
      for (unsigned short ii = 0; ii < aTj.Pts.size(); ++ii)
        PrintTP(someText, slc, ii, aTj.StepDir, aTj.Pass, aTj.Pts[ii]);
    }
    else {
      // print just one
      PrintTP(someText, slc, ipt, aTj.StepDir, aTj.Pass, aTj.Pts[ipt]);
    }
  } // PrintAllTraj

  //////////////////////////////////////////
  void
  PrintTrajectory(std::string someText,
                  const TCSlice& slc,
                  const Trajectory& tj,
                  unsigned short tPoint)
  {
    // prints one or all trajectory points on tj

    if (tPoint == USHRT_MAX) {
      if (tj.ID < 0) {
        mf::LogVerbatim myprt("TC");
        myprt << someText << " ";
        myprt << "Work:   UID " << tj.UID << "    CTP " << tj.CTP << " StepDir " << tj.StepDir
              << " PDG " << tj.PDGCode << " slc.vtxs " << tj.VtxID[0] << " " << tj.VtxID[1]
              << " nPts " << tj.Pts.size() << " EndPts " << tj.EndPt[0] << " " << tj.EndPt[1];
        myprt << " MCSMom " << tj.MCSMom;
        myprt << " EndFlags " << PrintEndFlag(tj, 0) << " " << PrintEndFlag(tj, 1);
        myprt << " AlgMods:";
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
          if (tj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
      }
      else {
        mf::LogVerbatim myprt("TC");
        myprt << someText << " ";
        myprt << "slcID " << slc.ID << " T" << tj.ID << " uT" << tj.UID << " WorkID " << tj.WorkID
              << " StepDir " << tj.StepDir << " PDG " << tj.PDGCode << " VtxID " << tj.VtxID[0]
              << " " << tj.VtxID[1] << " nPts " << tj.Pts.size() << " EndPts " << tj.EndPt[0] << " "
              << tj.EndPt[1];
        myprt << " MCSMom " << tj.MCSMom;
        myprt << " EndFlags " << PrintEndFlag(tj, 0) << " " << PrintEndFlag(tj, 1);
        myprt << " AlgMods:";
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
          if (tj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
      }
      PrintTPHeader(someText);
      for (unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt)
        PrintTP(someText, slc, ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
    }
    else {
      // just print one traj point
      if (tPoint > tj.Pts.size() - 1) {
        mf::LogVerbatim("TC") << "Can't print non-existent traj point " << tPoint;
        return;
      }
      PrintTP(someText, slc, tPoint, tj.StepDir, tj.Pass, tj.Pts[tPoint]);
    }
  } // PrintTrajectory

  //////////////////////////////////////////
  void
  PrintTPHeader(std::string someText)
  {
    mf::LogVerbatim("TC") << someText
                          << " TRP     CTP  Ind  Stp Delta  RMS    Ang C   Err  Dir0  Dir1      Q  "
                             "  AveQ  Pull FitChi  NTPF KinkSig  Hits ";
  } // PrintTPHeader

  ////////////////////////////////////////////////
  void
  PrintTP(std::string someText,
          const TCSlice& slc,
          unsigned short ipt,
          short dir,
          unsigned short pass,
          const TrajPoint& tp)
  {
    mf::LogVerbatim myprt("TC");
    myprt << someText << " TRP" << std::fixed;
    myprt << pass;
    if (dir > 0) { myprt << "+"; }
    else {
      myprt << "-";
    }
    myprt << std::setw(6) << tp.CTP;
    myprt << std::setw(5) << ipt;
    myprt << std::setw(5) << tp.Step;
    myprt << std::setw(6) << std::setprecision(2) << tp.Delta;
    myprt << std::setw(6) << std::setprecision(2) << tp.DeltaRMS;
    myprt << std::setw(6) << std::setprecision(2) << tp.Ang;
    myprt << std::setw(2) << tp.AngleCode;
    myprt << std::setw(6) << std::setprecision(2) << tp.AngErr;
    myprt << std::setw(6) << std::setprecision(2) << tp.Dir[0];
    myprt << std::setw(6) << std::setprecision(2) << tp.Dir[1];
    myprt << std::setw(7) << (int)tp.Chg;
    myprt << std::setw(8) << (int)tp.AveChg;
    myprt << std::setw(6) << std::setprecision(1) << tp.ChgPull;
    myprt << std::setw(7) << tp.FitChi;
    myprt << std::setw(6) << tp.NTPsFit;
    myprt << std::setw(7) << std::setprecision(3) << tp.KinkSig;
    // print the hits associated with this traj point
    if (tp.Hits.size() > 16) {
      // don't print too many hits (e.g. from a shower Tj)
      myprt << " " << tp.Hits.size() << " shower hits";
    }
    else {
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
        myprt << " " << hit.WireID().Wire << ":" << (int)hit.PeakTime();
        if (tp.UseHit[ii]) {
          // Distinguish used hits from nearby hits
          myprt << "_";
        }
        else {
          myprt << "x";
        }
        myprt << "T" << slc.slHits[iht].InTraj;
      } // iht
      if (tp.InPFP > 0) myprt << " inP" << tp.InPFP;
    }
    // print Environment
    if (tp.Environment.any()) myprt << " Env: " << TPEnvString(tp);
  } // PrintTP

  /////////////////////////////////////////
  std::string
  TPEnvString(const TrajPoint& tp)
  {
    // Print environment bits in human-readable format
    std::string str = "";
    for (unsigned short ib = 0; ib < 8; ++ib) {
      // There aren't any bit names for Environment_t
      if (!tp.Environment[ib]) continue;
      if (ib == kEnvNotGoodWire) str += " NoGdwire";
      if (ib == kEnvNearMuon) str += " NearMuon";
      if (ib == kEnvNearShower) str += " NearShower";
      if (ib == kEnvOverlap) str += " Overlap";
      if (ib == kEnvUnusedHits) str += " UnusedHits";
      if (ib == kEnvNearSrcHit) str += " NearSrcHit";
      if (ib == kEnvFlag) str += " Flag";
    } // ib
    return str;
  } // TPEnvironment

  /////////////////////////////////////////
  void
  PrintPFP(std::string someText, TCSlice& slc, const PFPStruct& pfp, bool printHeader)
  {
    mf::LogVerbatim myprt("TC");
    if (printHeader) {
      myprt << someText;
      myprt << "  PFP sVx  ________sPos_______ EF _______sDir______ ____sdEdx_____ eVx  "
               "________ePos_______ EF _______eDir______ ____edEdx____   Len nTp3 MCSMom ShLike? "
               "PDG Par Prim\n";
    }
    myprt << someText;
    std::string pid = "P" + std::to_string(pfp.ID);
    myprt << std::setw(5) << pid;
    // start and end stuff
    for (unsigned short end = 0; end < 2; ++end) {
      myprt << std::setw(4) << pfp.Vx3ID[end];
      myprt << std::fixed << std::right << std::setprecision(1);
      auto tp3d = EndTP3D(pfp, end);
      auto pos = tp3d.Pos;
      myprt << std::setw(7) << pos[0];
      myprt << std::setw(7) << pos[1];
      myprt << std::setw(7) << pos[2];
      // print characters that encode the EndFlag
      std::string ef;
      if (pfp.EndFlag[end][kEndOutFV]) { ef = "O"; }
      else {
        ef = "I";
      }
      if (pfp.EndFlag[end][kEndBragg]) ef += "B";
      myprt << std::setw(6) << ef;
      myprt << std::fixed << std::right << std::setprecision(2);
      auto dir = tp3d.Dir;
      myprt << std::setw(6) << dir[0];
      myprt << std::setw(6) << dir[1];
      myprt << std::setw(6) << dir[2];
      for (auto& dedx : pfp.dEdx[end]) {
        if (dedx < 50) { myprt << std::setw(5) << std::setprecision(1) << dedx; }
        else {
          myprt << std::setw(5) << std::setprecision(0) << dedx;
        }
      } // dedx
      if (pfp.dEdx[end].size() < 3) {
        for (size_t i = 0; i < 3 - pfp.dEdx[end].size(); ++i) {
          myprt << std::setw(6) << ' ';
        }
      }
    } // startend
    // global stuff
    float length = Length(pfp);
    if (length < 100) { myprt << std::setw(5) << std::setprecision(1) << length; }
    else {
      myprt << std::setw(5) << std::setprecision(0) << length;
    }
    myprt << std::setw(5) << std::setprecision(2) << pfp.TP3Ds.size();
    myprt << std::setw(7) << MCSMom(slc, pfp.TjIDs);
    myprt << std::setw(5) << IsShowerLike(slc, pfp.TjIDs);
    myprt << std::setw(5) << pfp.PDGCode;
    myprt << "      NA";
    myprt << std::setw(4) << pfp.ParentUID;
    myprt << std::setw(5) << PrimaryUID(slc, pfp);
    if (!pfp.TjIDs.empty()) {
      for (auto& tjID : pfp.TjIDs)
        myprt << " T" << tjID;
    }
    if (!pfp.DtrUIDs.empty()) {
      myprt << " dtrs";
      for (auto& dtrUID : pfp.DtrUIDs)
        myprt << " P" << dtrUID;
    }
  } // PrintPFP

  /////////////////////////////////////////
  void
  PrintPFPs(std::string someText, TCSlice& slc)
  {
    if (slc.pfps.empty()) return;

    mf::LogVerbatim myprt("TC");
    myprt << someText;
    myprt
      << "  PFP sVx  ________sPos_______  ______sDir______  ______sdEdx_____ eVx  "
         "________ePos_______  ______eDir______  ______edEdx_____ BstPln PDG TruPDG Par Prim E*P\n";
    bool printHeader = true;
    for (auto& pfp : slc.pfps) {
      PrintPFP(someText, slc, pfp, printHeader);
      printHeader = false;
    } // im

  } // PrintPFPs

  /////////////////////////////////////////
  std::string
  PrintEndFlag(const PFPStruct& pfp, unsigned short end)
  {
    if (end > 1) return "Invalid end";
    std::string tmp;
    bool first = true;
    for (unsigned short ib = 0; ib < EndFlagNames.size(); ++ib) {
      if (pfp.EndFlag[end][ib]) {
        if (first) {
          tmp = std::to_string(end) + ":" + EndFlagNames[ib];
          first = false;
        }
        else {
          tmp += "," + EndFlagNames[ib];
        }
      }
    } // ib
    if (first) tmp = " none";
    return tmp;
  } // PrintEndFlag

  /////////////////////////////////////////
  std::string
  PrintEndFlag(const Trajectory& tj, unsigned short end)
  {
    if (end > 1) return "Invalid end";
    std::string tmp;
    bool first = true;
    for (unsigned short ib = 0; ib < EndFlagNames.size(); ++ib) {
      if (tj.EndFlag[end][ib]) {
        if (first) {
          tmp = std::to_string(end) + ":" + EndFlagNames[ib];
          first = false;
        }
        else {
          tmp += "," + EndFlagNames[ib];
        }
      }
    } // ib
    return tmp;
  } // PrintEndFlag

  /////////////////////////////////////////
  std::string
  PrintHitShort(const TCHit& tch)
  {
    if (tch.allHitsIndex > (*evt.allHits).size() - 1) return "NA";
    auto& hit = (*evt.allHits)[tch.allHitsIndex];
    return std::to_string(hit.WireID().Plane) + ":" + std::to_string(hit.WireID().Wire) + ":" +
           std::to_string((int)hit.PeakTime());
  } // PrintHit

  /////////////////////////////////////////
  std::string
  PrintHit(const TCHit& tch)
  {
    if (tch.allHitsIndex > (*evt.allHits).size() - 1) return "NA";
    auto& hit = (*evt.allHits)[tch.allHitsIndex];
    return std::to_string(hit.WireID().Plane) + ":" + std::to_string(hit.WireID().Wire) + ":" +
           std::to_string((int)hit.PeakTime()) + "_" + std::to_string(tch.InTraj);
  } // PrintHit

  /////////////////////////////////////////
  std::string
  PrintPos(const TCSlice& slc, const TrajPoint& tp)
  {
    return std::to_string(DecodeCTP(tp.CTP).Plane) + ":" + PrintPos(slc, tp.Pos);
  } // PrintPos

  /////////////////////////////////////////
  std::string
  PrintPos(const TCSlice& slc, const Point2_t& pos)
  {
    unsigned int wire = 0;
    if (pos[0] > -0.4) wire = std::nearbyint(pos[0]);
    int time = std::nearbyint(pos[1] / tcc.unitsPerTick);
    return std::to_string(wire) + ":" + std::to_string(time);
  } // PrintPos


  ////////////////////////////////////////////////

} // namespace tca

