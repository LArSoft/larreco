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
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "TVector3.h"

namespace tca {
  
  // some functions to handle the CTP_t type
  typedef unsigned int CTP_t;
  constexpr unsigned int Tpad = 10; // alignment for CTP sub-items - TPC
  constexpr unsigned int Cpad = 10000; // alignment for CTP sub-items - Cryostat
  
  inline CTP_t EncodeCTP(unsigned int cryo, unsigned int tpc, unsigned int plane) { return cryo * Cpad + tpc * Tpad + plane; }
  inline CTP_t EncodeCTP(const geo::PlaneID& planeID) { return EncodeCTP(planeID.Cryostat, planeID.TPC, planeID.Plane); }
  inline CTP_t EncodeCTP(const geo::WireID& wireID) { return EncodeCTP(wireID.Cryostat, wireID.TPC, wireID.Plane); }
  geo::PlaneID DecodeCTP(CTP_t CTP);

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
    short Topo {0}; 			// 0 = end0-end0, 1 = end0(1)-end1(0), 2 = end1-end1, 3 = CI3DV, 4 = C3DIVIG, 5 = FHV, 6 = FHV2, 7 = SHCH
    CTP_t CTP {0};
    unsigned short ID {0};          ///< set to 0 if killed
    unsigned short Vtx3ID {0};
    float Score {0};
    float TjChgFrac {0};            ///< Fraction of charge near the vertex that is from hits on the vertex Tjs
    std::bitset<16> Stat {0};        ///< Vertex status bits using kVtxBit_t
  };
  
  typedef enum {
    kVtxTrjTried,     ///< FindVtxTraj algorithm tried
    kFixed,           ///< vertex position fixed manually - no fitting done
    kOnDeadWire,
    kVtxRefined,
    kHiVx3Score,      ///< matched to a high-score 3D vertex
    kVtxTruMatch,      ///< tagged as a vertex between Tjs that are matched to MC truth neutrino interaction particles
    kVtxMerged,
    kVtxBitSize     ///< don't mess with this line
  } VtxBit_t;
  
  /// struct of temporary 3D vertices
  struct Vtx3Store {
    float X {0};                    // x position
    float XErr {0.5};                 // x position error
    float Y {0};                    // y position
    float YErr {0.5};                 // y position error
    float Z {0};                    // z position
    float ZErr {0.5};                 // z position error
    float Score {0};
    short Wire {-1};                 // wire number for an incomplete 3D vertex
    geo::TPCID TPCID;
    std::array<unsigned short, 3> Vx2ID {{0}}; // List of 2D vertex IDs in each plane
    unsigned short ID {0};          // 0 = obsolete vertex
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
    int WorkID {0};
    int ParentTrajID {0};     ///< ID of the parent
    float AveChg {0};                   ///< Calculated using ALL hits
    float ChgRMS {0.5};                 /// Normalized RMS using ALL hits. Assume it is 50% to start
    short MCSMom {-1};         //< Crude 2D estimate to use for shower-like vs track-like discrimination
    int TruPDG {0};                    ///< MC truth
    int TruKE {0};                     ///< MeV
    float EffPur {0};                     ///< Efficiency * Purity
    std::array<float, 2> dEdx {{0,0}};      ///< dE/dx for 3D matched trajectories
    std::array<unsigned short, 2> VtxID {{0,0}};      ///< ID of 2D vertex
    std::array<unsigned short, 2> EndPt {{0,0}}; ///< First and last point in the trajectory that has charge
    int ID;
    unsigned short PDGCode {0};            ///< shower-like or track-like {default is track-like}
    unsigned int ClusterIndex {USHRT_MAX};   ///< Index not the ID...
    unsigned short Pass {0};            ///< the pass on which it was created
    short StepDir {0};                 ///< -1 = going US (CC proper order), 1 = going DS
    short TjDir {0};                     ///< direction determined by dQ/ds, delta ray direction, etc
                                        ///< 1 = in the StepDir direction, -1 in the opposite direction, 0 = don't know
    unsigned short MCPartListIndex {USHRT_MAX};
    unsigned short NNeighbors {0};    /// number of neighbors within window defined by ShowerTag
    std::array<std::bitset<8>, 2> StopFlag {};  // Bitset that encodes the reason for stopping
  };
  
  // Local version of recob::Hit
  struct TCHit {
    raw::TDCtick_t StartTick {0};
    raw::TDCtick_t EndTick {0};
    float PeakTime {0};
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
    int InTraj {0};
    unsigned short MCPartListIndex {USHRT_MAX};
  };

  // Struct for 3D trajectory matching
  struct MatchStruct {
    // IDs of Trajectories that match in all planes
    std::vector<int> TjIDs;
    // Count of the number of X-matched hits
    int Count {0};                    // Set to 0 if matching failed
    // Start is 0, End is 1
    std::array<std::array<float, 3>, 2> XYZ;        // XYZ position at both ends (cm)
    std::array<TVector3, 2> Dir;
    std::array<TVector3, 2> DirErr;
    std::array<std::vector<float>, 2> dEdx;
    std::array<std::vector<float>, 2> dEdxErr;
    std::array<unsigned short, 2> Vx3ID {0, 0};
    int BestPlane {INT_MAX};
    // stuff for constructing the PFParticle
    int PDGCode {0};
    std::vector<size_t> DtrIndices;
    size_t ParentMSIndex {0};       // Parent MatchStruct index (or index of self if no parent exists)
    geo::TPCID TPCID;
    float EffPur {0};                     ///< Efficiency * Purity
    unsigned short MCPartListIndex {USHRT_MAX};
  };

  struct ShowerPoint {
    std::array<float, 2> Pos;       // Hit Position in the normal coordinate system
    std::array<float, 2> RotPos;    // Position rotated into the shower coordinate system (0 = along, 1 = transverse)
    float Chg {0};                      // Charge of this point
    unsigned int HitIndex;                       // the hit index
    unsigned short TID;             // The ID of the tj the point (hit) is in. TODO eliminate this redundant variable
  };

  // A temporary structure that defines a 2D shower-like cluster of trajectories
  struct ShowerStruct {
    CTP_t CTP;
    int ShowerTjID {0};      // ID of the shower Trajectory composed of many InShower Tjs
    std::vector<int> TjIDs;  // list of InShower Tjs
    std::vector<int> NearTjIDs;   // list of Tjs that are not InShower but satisfy the maxSep cut
    std::vector<int> MatchedTjIDs;  /// list of Tjs in the other planes that are 3D matched to Tjs in this shower
    std::vector<ShowerPoint> ShPts;    // Trajectory points inside the shower
    float Angle {0};                   // Angle of the shower axis
    float AngleErr {3};                 // Error
    float AspectRatio {1};              // The ratio of charge weighted transverse/longitudinal positions
    float DirectionFOM {1};
    std::vector<std::array<float, 2>> Envelope; // Vertices of a polygon that encompasses the shower
    float EnvelopeArea {0};
    float ChgDensity {0};                   // Charge density inside the Envelope
    float Energy {0};
    float ParentFOM {10};
    int ParentID {0};  // The ID of an external parent Tj that was added to the shower
    bool NeedsUpdate {false};       // This is set true whenever the shower needs to be updated
    unsigned short TruParentID {0};
  };
  
  // Shower variables filled in MakeShowers. These are in cm and radians
  struct ShowerStruct3D {
    TVector3 Dir;
    TVector3 DirErr;
    TVector3 Pos;
    TVector3 PosErr;
    double Len {1};
    double OpenAngle {0.2};
    std::vector<double> Energy;
    std::vector<double> EnergyErr;
    std::vector<double> MIPEnergy;
    std::vector<double> MIPEnergyErr;
    std::vector<double> dEdx;
    std::vector<double> dEdxErr;
    int BestPlane;
    int ID;
    std::vector<unsigned short> TjIDs;
    std::vector<unsigned int> Hits;
  };

  struct ShowerTreeVars {

    std::vector<float> BeginWir;   // begin wire
    std::vector<float> BeginTim;   // begin tick
    std::vector<float> BeginAng;   // begin angle
    std::vector<float> BeginChg;   // beginning average charge
    std::vector<short> BeginVtx;   // ID of begin vertex
    std::vector<float> EndWir;   // end wire
    std::vector<float> EndTim;   // end tick
    std::vector<float> EndAng;   // end angle
    std::vector<float> EndChg;   // ending average charge
    std::vector<short> EndVtx;   //ID of end vertex
    
    std::vector<short> MCSMom;

    std::vector<short> PlaneNum; 
    
    std::vector<int> ShowerID; // shower ID associated w/ trajectory. -1 = no shower
    std::vector<int> IsShowerParent;
    std::vector<int> StageNum; // stage of reconstruction

    // envelope information
    std::vector<float> Envelope;
    std::vector<int>EnvPlane;
    std::vector<int>EnvStage;
    std::vector<int>EnvShowerID;

    int nStages;
    unsigned short nPlanes;

  };

  // Algorithm modification bits
  typedef enum {
    kMaskHits,
    kMaskBadTPs,
    kCTKink,        ///< kink found in CheckTraj
    kCTStepChk,
    kTryWithNextPass,
    kRvPrp,
    kChkHiMultHits,
    kSplit,
    kComp3DVx,
    kComp3DVxIG,
    kHED, // High End Delta
    kHamVx,
    kHamVx2,
    kJunkTj,
    kKilled,
    kMerge,
    kTEP,
    kCHMEH,
    kFillGap,
    kUseGhostHits,
    kChkInTraj,
    kStopBadFits,
    kFixBegin,
    kFixEnd,
    kUUH,
    kVtxTj,
    kRefineVtx,
    kNoKinkChk,
    kSoftKink,
    kChkStop,
    kFTBRvProp,
    kStopAtTj,
    kMat3D,
    kMat3DMerge,
    kTjHiVx3Score,
    kVtxHitsSwap,
    kSplitHiChgHits,
    kInShower,
    kShowerTj,
    kShwrParent,
    kMergeOverlap,
    kMergeSubShowers,
    kMergeNrShowers,
    kAlgBitSize     ///< don't mess with this line
  } AlgBit_t;
  
  // Stop flag bits
  typedef enum {
    kSignal,
    kAtKink,
    kAtVtx,
    kBragg,
    kAtTj,
    kFlagBitSize     ///< don't mess with this line
  } StopFlag_t; 
  
  extern const std::vector<std::string> AlgBitNames;
  extern const std::vector<std::string> StopFlagNames;
  extern const std::vector<std::string> VtxBitNames;
  
  struct TjStuff {
    // These variables don't change in size from event to event
    float UnitsPerTick;     ///< scale factor from Tick to WSE equivalent units
    geo::TPCID TPCID;
    std::vector<unsigned int> NumWires;
    std::vector<float> MaxPos0;
    std::vector<float> MaxPos1;
    std::vector<unsigned int> FirstWire;    ///< the first wire with a hit
    std::vector<unsigned int> LastWire;      ///< the last wire with a hit
    unsigned short NumPlanes;
    float XLo; // fiducial volume of the current tpc
    float XHi;
    float YLo;
    float YHi;
    float ZLo;
    float ZHi;
    std::vector<float> AveHitRMS;      ///< average RMS of an isolated hit
    // The variables below do change in size from event to event
    ShowerTreeVars stv; // 
    bool SaveShowerTree;
    
    std::vector<Trajectory> allTraj; ///< vector of all trajectories in each plane
    std::vector<TCHit> fHits;
    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector<std::vector< std::pair<int, int>>> WireHitRange;
    std::vector<float> AngleRanges; ///< list of max angles for each angle range
    std::vector<short> inClus;    ///< Hit -> cluster ID (0 = unused)
    std::vector< ClusterStore > tcl; ///< the clusters we are creating
    std::vector< VtxStore > vtx; ///< 2D vertices
    std::vector< Vtx3Store > vtx3; ///< 3D vertices
    std::vector<MatchStruct> matchVec; ///< 3D matching vector
    std::vector<unsigned short> matchVecPFPList;  /// list of matchVec entries that will become PFPs
    std::vector<ShowerStruct> cots;       // Clusters of Trajectories that define 2D showers
    std::vector<ShowerStruct3D> showers;  // 3D showers
    std::vector<float> Vertex2DCuts; ///< Max position pull, max Position error rms
    float Vertex3DChiCut;   ///< 2D vtx -> 3D vtx matching cut (chisq/dof)
    std::vector<float> VertexScoreWeights;
    std::vector<short> DeltaRayTag; ///< min length, min MCSMom and min separation (WSE) for a delta ray tag
    std::vector<short> MuonTag; ///< min length and min MCSMom for a muon tag
    std::vector<float> ShowerTag; ///< [min MCSMom, max separation, min # Tj < separation] for a shower tag
    std::vector<float> Match3DCuts;  ///< 3D matching cuts
    std::vector<float> MatchTruth;     ///< Match to MC truth
    std::vector<const simb::MCParticle*> MCPartList;
    std::bitset<64> UseAlg;  ///< Allow user to mask off specific algorithms
    const geo::GeometryCore* geom;
    const detinfo::DetectorProperties* detprop;
    calo::CalorimetryAlg* caloAlg;
    bool IgnoreNegChiHits;
   };

} // namespace tca

#endif // ifndef TRAJCLUSTERALGDATASTRUCT_H
