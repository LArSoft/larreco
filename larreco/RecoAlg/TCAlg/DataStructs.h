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
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

namespace tca {
  
  // some functions to handle the CTP_t type
  typedef unsigned int CTP_t;
  constexpr unsigned int CTPpad = 1000; // alignment for CTP sub-items
  inline CTP_t EncodeCTP(unsigned int cryo, unsigned int tpc, unsigned int plane) { return cryo * CTPpad*CTPpad + tpc * CTPpad + plane; }
  inline CTP_t EncodeCTP(const geo::PlaneID& planeID) { return EncodeCTP(planeID.Cryostat, planeID.TPC, planeID.Plane); }
  inline CTP_t EncodeCTP(const geo::WireID& wireID) { return EncodeCTP(wireID.Cryostat, wireID.TPC, wireID.Plane); }
  inline geo::PlaneID DecodeCTP(CTP_t CTP) { return { CTP / (CTPpad*CTPpad), CTP / CTPpad % CTPpad, CTP % CTPpad }; }

  /// @{
  /// @name Data structures for the reconstruction results
  
  /// struct of temporary clusters
  struct ClusterStore {
    short ID {0};         // Cluster ID. ID < 0 = abandoned cluster
    CTP_t CTP {0};        // Cryostat/TPC/Plane code
    unsigned short PDGCode {0}; // PDG-like code shower-like or line-like
    unsigned short ParentCluster {0};
    float BeginWir {0};   // begin wire
    float BeginTim {0};   // begin tick
    float BeginAng {0};   // begin angle
    float BeginChg {0};   // beginning average charge
    short BeginVtx {-1}; 	// ID of the Begin vertex
    float EndWir {0};     // end wire
    float EndTim {0};     // end tick
    float EndAng {0};     // end angle
    float EndChg {0};     // ending average charge
    short EndVtx {-1};     // ID of the end vertex
    std::vector<unsigned int> tclhits; // hits on the cluster
  }; // ClusterStore
  
  /// struct of temporary 2D vertices (end points)
  struct VtxStore {
    std::array<float, 2> Pos {{0,0}};
    std::array<float, 2> PosErr {{2,1}};
    unsigned short NTraj {0};  // = 0 for abandoned vertices
    unsigned short Pass {0};   // Pass in which this vertex was created
    float ChiDOF {0};
    short Topo {0}; 			// 0 = end0-end0, 1 = end0(1)-end1(0), 2 = end1-end1, 5 = Star, 6 = hammer, 7 = photon conversion, 8 = dead region
    CTP_t CTP {0};
    unsigned short ID {0};
    short Ptr3D {SHRT_MAX};
    std::bitset<16> Stat {0};        ///< Vertex status bits using kVtxBit_t
  };
  
  typedef enum {
    kFixed,           ///< vertex position fixed manually - no fitting done
    kVtxTrjTried,     ///< FindVtxTraj algorithm tried
    kOnDeadWire,
    kVtxRefined,
    kVtxBitSize     ///< don't mess with this line
  } VtxBit_t;
  
  /// struct of temporary 3D vertices
  struct Vtx3Store {
    float X {0};                    // x position
    float XErr {0};                 // x position error
    float Y {0};                    // y position
    float YErr {0};                 // y position error
    float Z {0};                    // z position
    float ZErr {0};                 // z position error
    short Wire {-1};                 // wire number for an incomplete 3D vertex
    unsigned short CStat {0};
    unsigned short TPC {0};
    std::array<short, 3> Ptr2D {{-1, -1, -1}}; // pointers to 2D vertices in each plane
  };
  
  struct TrajPoint {
    CTP_t CTP {0};                   ///< Cryostat, TPC, Plane code
    std::array<float, 2> HitPos {{0,0}}; // Charge weighted position of hits in wire equivalent units
    std::array<float, 2> Pos {{0,0}}; // Trajectory position in wire equivalent units
    std::array<double, 2> Dir {{0,0}}; // Direction cosines in the StepDir direction
    double HitPosErr2 {0};         // Uncertainty^2 of the hit position perpendiclar to the direction
    // HitPosErr2 < 0 = HitPos not defined because no hits used
    double Ang {0};                // Trajectory angle (-pi, +pi)
    double AngErr {0.1};             // Trajectory angle error
    float Chg {0};                // Charge
    float AveChg {-1};             // Average charge of last ~20 TPs
    float ChgPull {0.1};          //  = (Chg - AveChg) / ChgRMS
    float Delta {0};              // Deviation between trajectory and hits (WSE)
    float DeltaRMS {0.02};           // RMS of Deviation between trajectory and hits (WSE)
    float FitChi {0};             // Chi/DOF of the fit
    unsigned short NTPsFit {2}; // Number of trajectory points fitted to make this point
    unsigned short Step {0};      // Step number at which this TP was created
    unsigned short AngleCode {0};          // 0 = small angle, 1 = large angle, 2 = very large angle
    std::vector<unsigned int> Hits; // vector of fHits indices
    std::bitset<16> UseHit {0};   // set true if the hit is used in the fit
  };
  
  // Global information for the trajectory
  struct Trajectory {
    std::vector<TrajPoint> Pts;    ///< Trajectory points
    CTP_t CTP {0};                      ///< Cryostat, TPC, Plane code
    std::bitset<64> AlgMod;        ///< Bit set if algorithm AlgBit_t modifed the trajectory
    unsigned short PDGCode {0};            ///< shower-like or track-like {default is track-like}
    unsigned short ParentTrajID {0};     ///< ID of the parent (if PDG = 12)
    float AveChg {0};                   ///< Calculated using ALL hits
    float ChgRMS {0.5};                 /// Normalized RMS using ALL hits. Assume it is 50% to start
    short MCSMom {-1};         //< Crude 2D estimate to use for shower-like vs track-like discrimination
    int TruPDG {0};                    ///< MC truth
    int TruKE {0};                     ///< MeV
    float EffPur {0};                     ///< Efficiency * Purity
    std::array<unsigned short, 2> VtxID {{0,0}};      ///< ID of 2D vertex
    std::array<unsigned short, 2> EndPt {{0,0}}; ///< First and last point in the trajectory that has charge
    short ID;
    unsigned short ClusterIndex {USHRT_MAX};   ///< Index not the ID...
    unsigned short Pass {0};            ///< the pass on which it was created
    short StepDir {0};                 ///< -1 = going US (CC proper order), 1 = going DS
    short Dir {0};                     ///< direction determined by dQ/ds, delta ray direction, etc
                                        ///< 1 (-1) = in (opposite to)the  StepDir direction, 0 = don't know
    int WorkID {0};
    std::array<std::bitset<8>, 2> StopFlag {};  // Bitset that encodes the reason for stopping
  };
  
  // Local version of recob::Hit
  struct TCHit {
    raw::TDCtick_t StartTick {0};
    raw::TDCtick_t EndTick {0};
    float PeakTime {0};     ///< Note that this the time in WSE units - NOT ticks
    float SigmaPeakTime {1};
    float PeakAmplitude {1};
    float SigmaPeakAmp {1};
    float Integral {1};
    float SigmaIntegral {1};
    float RMS {1};
    float GoodnessOfFit {0};
    unsigned short NDOF {0};
    unsigned short Multiplicity {1};
    unsigned short LocalIndex {0};
    geo::WireID WireID;
    short InTraj {0};
  };

  // Struct for 3D trajectory matching
  struct MatchStruct {
    // IDs of Trajectories that match in all planes
    std::vector<unsigned short> TjIDs;
    std::vector<unsigned int> SeedHit;
    // Count of the number of time-matched hits
    int Count;
  };

  // Algorithm modification bits
  typedef enum {
    kMaskHits,
    kCTKink,        ///< kink found in CheckWork
    kCTStepChk,
    kTryWithNextPass,
    kRevProp,
    kChkHiMultHits,
    kSplitTraj,
    kComp3DVx,
    kComp3DVxIG,
    kHiEndDelta,
    kHamVx,
    kHamVx2,
    kJunkTj,
    kKilled,
    kEndMerge,
    kTrimEndPts,
    kChkHiMultEndHits,
    kFillGap,
    kUseGhostHits,
    kChkInTraj,
    kFixBegin,
    kFixEnd,
    kUseUnusedHits,
    kVtxTj,
    kRefineVtx,
    kMaskBadTPs,
    kNoKinkChk,
    kSoftKink,
    kChkStop,
    kChkAllStop,
    kFTBRevProp,
    kStopAtTj,
    kMatch3D,
    kVtxHitsSwap,
    kAlgBitSize     ///< don't mess with this line
  } AlgBit_t;
  
  // Stop flag bits
  typedef enum {
    kSignal,
    kAtKink,
    kAtVtx,
    kBragg,
    kRvPrp,
    kAtTj,
    kFlagBigSize     ///< don't mess with this line
  } StopFlag_t; 
  
  extern const std::vector<std::string> AlgBitNames;
  extern const std::vector<std::string> StopFlagNames;
  
  struct TjStuff {
    // These variables don't change in size from event to event
    float UnitsPerTick;     ///< scale factor from Tick to WSE equivalent units
    std::vector<unsigned int> NumWires;
    std::vector<float> MaxPos0;
    std::vector<float> MaxPos1;
    std::vector<unsigned int> FirstWire;    ///< the first wire with a hit
    std::vector<unsigned int> LastWire;      ///< the last wire with a hit
    // The variables below do change in size from event to event
    std::vector<Trajectory> allTraj; ///< vector of all trajectories in each plane
    std::vector<TCHit> fHits;
    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector<std::vector< std::pair<int, int>>> WireHitRange;
    unsigned short WireHitRangeCstat;
    unsigned short WireHitRangeTPC;
    std::vector<short> inClus;    ///< Hit -> cluster ID (0 = unused)
    std::vector< ClusterStore > tcl; ///< the clusters we are creating
    std::vector< VtxStore > vtx; ///< 2D vertices
    std::vector< Vtx3Store > vtx3; ///< 3D vertices
    std::vector<MatchStruct> matchVec; ///< 3D matching vector
    std::vector<std::vector<unsigned short>> MatchedTjIDs;
    std::vector<std::vector<unsigned short>> MatchedClusters;
    unsigned short NumPlanes;
    float YLo; // fiducial volume of the current tpc
    float YHi;
    float ZLo;
    float ZHi;
   };

} // namespace tca

#endif // ifndef TRAJCLUSTERALGDATASTRUCT_H
