////////////////////////////////////////////////////////////////////////
//
//
// TCAlg data structs
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGDATASTRUCT_H
#define TRAJCLUSTERALGDATASTRUCT_H


// C/C++ standard libraries
#include <array>
#include <vector>
#include <bitset>

// LArSoft libraries
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "art/Persistency/Common/Ptr.h"
#include "lardata/RecoBase/Hit.h"

namespace tca {
  
  // some functions to handle the CTP_t type
  typedef unsigned int CTP_t;
  
  constexpr unsigned int CTPpad = 1000; // alignment for CTP sub-items
  
  inline CTP_t EncodeCTP
  (unsigned int cryo, unsigned int tpc, unsigned int plane)
  { return cryo * CTPpad*CTPpad + tpc * CTPpad + plane; }
  inline CTP_t EncodeCTP(const geo::PlaneID& planeID)
  { return EncodeCTP(planeID.Cryostat, planeID.TPC, planeID.Plane); }
  inline CTP_t EncodeCTP(const geo::WireID& wireID)
  { return EncodeCTP(wireID.Cryostat, wireID.TPC, wireID.Plane); }
  
  inline geo::PlaneID DecodeCTP(CTP_t CTP)
  { return { CTP / (CTPpad*CTPpad), CTP / CTPpad % CTPpad, CTP % CTPpad }; }
  
  /// @{
  /// @name Data structures for the reconstruction results
  
  /// struct of temporary clusters
  struct ClusterStore {
    short ID;         // Cluster ID. ID < 0 = abandoned cluster
    CTP_t CTP;        // Cryostat/TPC/Plane code
    unsigned short PDG; // PDG-like code shower-like or line-like
    unsigned short ParentCluster;
    float BeginWir;   // begin wire
    float BeginTim;   // begin tick
    float BeginAng;   // begin angle
    float BeginChg;   // beginning average charge
    short BeginVtx; 	// ID of the Begin vertex
    float EndWir;     // end wire
    float EndTim;     // end tick
    float EndAng;     // end angle
    float EndChg;     // ending average charge
    short EndVtx;     // ID of the end vertex
    std::vector<unsigned int> tclhits; // hits on the cluster
  }; // ClusterStore
  
  /// struct of temporary 2D vertices (end points)
  struct VtxStore {
    float Wire {0};
    float WireErr {2};
    float Time {0};
    float TimeErr {0.1};
    unsigned short NTraj {0};  // = 0 for abandoned vertices
    float ChiDOF {0};
    short Topo {0}; 			// 1 = US-US, 2 = US-DS, 3 = DS-US, 4 = DS-DS, 5 = Star, 6 = hammer, 7 = photon conversion
    CTP_t CTP {0};
    bool Fixed {false};                 // Vertex position fixed (should not be re-fit)
  };
  
  /// struct of temporary 3D vertices
  struct Vtx3Store {
    std::array<short, 3> Ptr2D; // pointers to 2D vertices in each plane
    float X;                    // x position
    float XErr;                 // x position error
    float Y;                    // y position
    float YErr;                 // y position error
    float Z;                    // z position
    float ZErr;                 // z position error
    short Wire;                 // wire number for an incomplete 3D vertex
    unsigned short CStat;
    unsigned short TPC;
    unsigned short ProcCode;
  };
  
  struct TrajPoint {
    CTP_t CTP {0};                   ///< Cryostat, TPC, Plane code
    std::array<float, 2> HitPos {{0,0}}; // Charge weighted position of hits in wire equivalent units
    std::array<float, 2> Pos {{0,0}}; // Trajectory position in wire equivalent units
    std::array<float, 2> Dir {{0,0}}; // Direction
    float HitPosErr2 {0};         // Uncertainty^2 of the hit position perpendiclar to the direction
    // HitPosErr2 < 0 = HitPos not defined because no hits used
    float Ang {0};                // Trajectory angle (-pi, +pi)
    float AngErr {0.1};             // Trajectory angle error
    float KinkAng {-1};            // Just what it says
    float Chg {0};                // Charge
    float AveChg {-1};             // Average charge of last ~20 TPs
    float ChgPull {0.1};          //  = (Chg - AveChg) / ChgRMS
    float Delta {0};              // Deviation between trajectory and hits (WSE)
    float DeltaRMS {0};           // RMS of Deviation between trajectory and hits (WSE)
    unsigned short NTPsFit {2}; // Number of trajectory points fitted to make this point
    unsigned short Step {0};      // Step number at which this TP was created
    float FitChi {0};             // Chi/DOF of the fit
    std::vector<unsigned int> Hits; // vector of fHits indices
    std::vector<bool> UseHit; // set true if the hit is used in the fit
  };
  
  // Global information for the trajectory
  struct Trajectory {
    short ID;
    CTP_t CTP {0};                      ///< Cryostat, TPC, Plane code
    unsigned short Pass {0};            ///< the pass on which it was created
    short StepDir {0};                 /// -1 = going US (CC proper order), 1 = going DS
    unsigned short ClusterIndex {USHRT_MAX};   ///< Index not the ID...
    std::bitset<32> AlgMod;        ///< Bit set if algorithm AlgBit_t modifed the trajectory
    unsigned short PDG {0};            ///< shower-like or track-like {default is track-like}
    unsigned short ParentTraj {USHRT_MAX};     ///< index of the parent (if PDG = 12)
    float AveChg {0};                   ///< Calculated using ALL hits
    float ChgRMS {1};                 /// Normalized RMS using ALL hits. Assume it is 100% to start
    int TruPDG {0};                    ///< MC truth
    int TruKE {0};                     ///< MeV
    float EffPur {0};                     ///< Efficiency * Purity
    std::array<short, 2> Vtx {{-1,-1}};      ///< Index of 2D vertex
    std::array<unsigned short, 2> EndPt {{0,0}}; ///< First and last point in the trajectory that has a hit
    std::vector<TrajPoint> Pts;    ///< Trajectory points
    std::vector<TrajPoint> EndTP {(2)} ;    ///< Trajectory point at each end for merging
  };
  
  // Trajectory "intersections" used to search for superclusters (aka showers)
  struct TrjInt {
    unsigned short itj1;
    unsigned short ipt1;
    unsigned short itj2;
    unsigned short ipt2;
    float sep2;   // separation^2 at closest point
    float dang;   // opening angle at closest point
    float vw;     // intersection wire
    float vt;     // intersection time
  };
  
  struct TjPairHitShare {
    // Trajectories in two different trials that share hits
    unsigned short iTrial;
    unsigned short iTj;
    unsigned short jTrial;
    unsigned short jTj;
    unsigned short nSameHits;
  };
  
  // Algorithm modification bits
  typedef enum {
    kMaskedWorkHits,
    kGottaKink,     ///< GottaKink found a kink
    kCWKink,        ///< kink found in CheckWork
    kCWStepChk,
    kMaybeDeltaRay,
    kModifyShortTraj,
    kTryWithNextPass,
    kRevProp,
    kRecovery1,
    kRecovery2,
    kManyHitsAdded,
    kSplitTraj,
    kComp3DVx,
    kHiEndDelta,
    kHammer2DVx,
    kJunkTj,
    kKilled,
    kStopAtVtx,
    kMerged,
    kTrimHits,
    kUseHiMultEndHits,
    kAlgBitSize     ///< don't mess with this line
  } AlgBit_t;
  
  extern const std::vector<std::string> AlgBitNames;
  
  struct TjStuff {
    std::vector<Trajectory> allTraj; ///< vector of all trajectories in each plane
    std::vector<short> inTraj;       ///< Hit -> trajectory ID (0 = unused)
    std::vector<art::Ptr<recob::Hit>> fHits;
    std::vector<short> inClus;    ///< Hit -> cluster ID (0 = unused)
    std::vector< ClusterStore > tcl; ///< the clusters we are creating
    std::vector< VtxStore > vtx; ///< 2D vertices
    std::vector< Vtx3Store > vtx3; ///< 3D vertices
    float UnitsPerTick;     ///< scale factor from Tick to WSE equivalent units
  };

} // namespace tca

#endif // ifndef TRAJCLUSTERALGDATASTRUCT_H
