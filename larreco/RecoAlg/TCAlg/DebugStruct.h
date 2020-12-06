////////////////////////////////////////////////////////////////////////
//
//
// TCAlg debug struct
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGDEBUGSTRUCT_H
#define TRAJCLUSTERALGDEBUGSTRUCT_H

#include <limits.h>

// LArSoft libraries
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace tca {

  struct DebugStuff {
    int Cryostat {0};    ///< Select Cryostat
    int TPC {0};    ///< Select TPC
    int Plane {-1}; ///< Select plane
    CTP_t CTP {UINT_MAX}; ///< set to an invalid CTP
    int Wire {-1};  ///< Select hit Wire for debugging
    int Tick {-1};   ///< Select hit PeakTime for debugging (< 0 for vertex finding)
    unsigned int Hit {UINT_MAX};    ///< set to the hit index in evt.allHits if a Plane:Wire:Tick match is found
    int WorkID {0}; ///< Select the StartWorkID for debugging
    int UID {0};     ///< or alternatively the trajectory unique ID
    unsigned int MVI {UINT_MAX}; ///< MatchVec Index for detailed 3D matching
    unsigned short MVI_Iter {USHRT_MAX}; ///< MVI iteration - see FindPFParticles
    int Slice {-1};
  };
  extern DebugStuff debug;

  bool DecodeDebugString(std::string ctpwt);
  void DumpTj();
  bool InTrajOK(TCSlice& slc, std::string someText);
  void PrintDebugMode();
  void PrintAll(detinfo::DetectorPropertiesData const& detProp, std::string someText);
  void PrintP(std::string someText, mf::LogVerbatim& myprt, PFPStruct& pfp, bool& printHeader);
  void Print3V(detinfo::DetectorPropertiesData const& detProp,
               std::string someText,
               mf::LogVerbatim& myprt,
               Vtx3Store& vx3,
               bool& printHeader);
  void Print2V(std::string someText, mf::LogVerbatim& myprt, VtxStore& vx2, bool& printHeader);
  void PrintT(std::string someText, mf::LogVerbatim& myprt, Trajectory& tj, bool& printHeader);
  std::string PackEndFlags(const Trajectory& tj, unsigned short end);
  void PrintTrajectory(std::string someText,
                       const TCSlice& slc,
                       const Trajectory& tj,
                       unsigned short tPoint);
  void PrintAllTraj(detinfo::DetectorPropertiesData const& detProp,
                    std::string someText,
                    TCSlice& slc,
                    unsigned short itj,
                    unsigned short ipt,
                    bool printVtx = true);
  void PrintTPHeader(std::string someText);
  void PrintTP(std::string someText,
               const TCSlice& slc,
               unsigned short ipt,
               short dir,
               unsigned short pass,
               const TrajPoint& tp);
  std::string TPEnvString(const TrajPoint& tp);
  void PrintPFP(std::string someText, TCSlice& slc, const PFPStruct& pfp, bool printHeader);
  void PrintPFPs(std::string someText, TCSlice& slc);
  // Print clusters after calling MakeAllTrajClusters
  void PrintClusters();
  // Print a single hit in the standard format
  std::string PrintHit(const TCHit& hit);
  std::string PrintHitShort(const TCHit& hit);
  // Print Trajectory position in the standard format
  std::string PrintPos(const TCSlice& slc, const TrajPoint& tp);
  std::string PrintPos(const TCSlice& slc, const Point2_t& pos);
  std::string PrintEndFlag(const Trajectory& tj, unsigned short end);
  std::string PrintEndFlag(const PFPStruct& pfp, unsigned short end);

} // namespace tca

#endif // ifndef TRAJCLUSTERALGDEBUGSTRUCT_H
