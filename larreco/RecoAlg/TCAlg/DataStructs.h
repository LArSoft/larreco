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

namespace tca {
  
  using Point3_t = std::array<double, 3>;
  using Vector3_t = std::array<double, 3>;
  using Point2_t = std::array<float, 2>;
  using Vector2_t = std::array<double, 2>;

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
  
/*
    Associations
    Hit.InTraj <-> tj.Pts.Hits
    Tj.ParentID -> ID of parent tj
    PFParticle.TjIDs -> IDs of tjs in each plane that define the PFParticle
    PFParticle.ParentID -> PFParticle ID of parent
    Shower.TjIDs -> IDs of InShower tjs, (tj.AlgMod[kInShower] set true)
    Shower.ShowerTjID -> ID of the shower tj (3 pts for start, chg center, end) (tj.AlgMod[kShowerTj] set true)
    Shower.ParentID -> ID of the tj identified as the shower parent (tj.AlgMod[kShwrParent] set true)
*/ 
  
  /// struct of temporary clusters
  struct ClusterStore {
    int ID {0};         // Cluster ID. ID < 0 = abandoned cluster
    CTP_t CTP {0};        // Cryostat/TPC/Plane code
    unsigned short PDGCode {0}; // PDG-like code shower-like or line-like
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
    Point2_t Pos {{0,0}};
    Point2_t PosErr {{2,1}};
    unsigned short NTraj {0};  
    unsigned short Pass {0};   // Pass in which this vertex was created
    float ChiDOF {0};
    // Topo: 0 = end0-end0, 1 = end0(1)-end1(0), 2 = end1-end1, 3 = CI3DV, 
    //       4 = C3DIVIG, 5 = FHV, 6 = FHV2, 7 = SHCH, 8 = CTBC, 9 = Junk, 10 = 3D split
    short Topo {0}; 			
    CTP_t CTP {0};
    unsigned short ID {0};          ///< set to 0 if killed
    unsigned short Vx3ID {0};
    float Score {0};
    float TjChgFrac {0};            ///< Fraction of charge near the vertex that is from hits on the vertex Tjs
    std::bitset<16> Stat {0};        ///< Vertex status bits using kVtxBit_t
  };
  
  typedef enum {
    kVtxTrjTried,     ///< FindVtxTraj algorithm tried
    kFixed,           ///< vertex position fixed manually - no fitting done
    kOnDeadWire,
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
    bool Primary {false};
    bool Neutrino {false};
  };
  
  // A temporary struct for matching trajectory points; 1 struct for each TP for
  // each trajectory. These are put into mallTraj which is then sorted by increasing xlo
  struct Tj2Pt{
    std::array<double, 2> dir;
    unsigned int wire;
    // x range spanned by hits on the TP
    float xlo;
    float xhi;
    CTP_t ctp;
    // the Trajectory ID
    unsigned short id;
    unsigned short ipt; // The trajectory point
    // the number of points in the Tj so that the minimum Tj length cut (MatchCuts[2]) can be made
    unsigned short npts;
    short score; // 0 = Tj with nice vertex, 1 = high quality Tj, 2 = normal, -1 = already matched
    bool inShower;
  };

  struct TrajPoint {
    CTP_t CTP {0};                   ///< Cryostat, TPC, Plane code
    Point2_t HitPos {{0,0}}; // Charge weighted position of hits in wire equivalent units
    Point2_t Pos {{0,0}}; // Trajectory position in wire equivalent units
    std::array<double, 2> Dir {{0,0}}; // Direction cosines in the StepDir direction
    double HitPosErr2 {0};         // Uncertainty^2 of the hit position perpendiclar to the direction
    // HitPosErr2 < 0 = HitPos not defined because no hits used
    double Ang {0};                // Trajectory angle (-pi, +pi)
    double AngErr {0.1};             // Trajectory angle error
    float Chg {0};                // Chargetj2pt
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
    std::bitset<4> Environment {0};    // TPEnvironment_t bitset that describes the environment, e.g. nearby showers or other Tjs
  };
  
  // Global information for the trajectory
  struct Trajectory {
    std::vector<TrajPoint> Pts;    ///< Trajectory points
    CTP_t CTP {0};                      ///< Cryostat, TPC, Plane code
    std::bitset<64> AlgMod;        ///< Bit set if algorithm AlgBit_t modifed the trajectory
    int WorkID {0};
    int ParentID {-1};     ///< ID of the parent, or the ID of the Tj this one was merged with if it is killed
    float AveChg {0};                   ///< Calculated using ALL hits
    float TotChg {0};                   ///< Total including an estimate for dead wires
    float ChgRMS {0.5};                 /// Normalized RMS using ALL hits. Assume it is 50% to start
    float DirFOM {0.5};         ///< confidence level that the Tj points are ordered correctly using  charge pattern
    short MCSMom {-1};         //< Crude 2D estimate to use for shower-like vs track-like discrimination
    float EffPur {0};                     ///< Efficiency * Purity
    Point2_t dEdx {{0,0}};      ///< dE/dx for 3D matched trajectories
    std::array<unsigned short, 2> VtxID {{0,0}};      ///< ID of 2D vertex
    std::array<unsigned short, 2> EndPt {{0,0}}; ///< First and last point in the trajectory that has charge
    int ID;
    unsigned short PDGCode {0};            ///< shower-like or track-like {default is track-like}
    unsigned int ClusterIndex {USHRT_MAX};   ///< Index not the ID...
    unsigned short Pass {0};            ///< the pass on which it was created
    short StepDir {0};                 ///< -1 = going US (-> small wire#), 1 = going DS (-> large wire#)
    unsigned int MCPartListIndex {UINT_MAX};
    std::array<std::bitset<8>, 2> StopFlag {};  // Bitset that encodes the reason for stopping
    bool NeedsUpdate {false};          ///< Set true when the Tj needs to be updated (only for the TP Environment right now)
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
    art::Ptr<recob::Hit> ArtPtr;
    unsigned short NDOF {0};
    unsigned short Multiplicity {1};
    unsigned short LocalIndex {0};
    int InTraj {0};
    unsigned int MCPartListIndex {UINT_MAX};
  };
  
  // struct used for TrajCluster 3D trajectory points
  struct TrajPoint3 {
    Point3_t Pos {0};
    Vector3_t Dir {0};
    std::vector<Tj2Pt> Tj2Pts;  // list of trajectory points
    float dEdx {0};             // The charge is stored here before dE/dx is calculated
    float dEdxErr {1};
    float ChiDOF {10};             // Chi/DOF of the fit < 0 = not valid
    unsigned short nPtsFit {2}; 
    bool IsValid {true};     // Is consistent with the position/angle of nearby space points
  };

  // Struct for 3D trajectory matching
  struct MatchStruct {
    // IDs of Trajectories that match in all planes
    std::vector<int> TjIDs;
    std::vector<float> TjCompleteness;  // fraction of TP points that are 3D-matched
    // Count of the number of X-matched hits and de-weight by angle
    float Count {0};                    // Set to 0 if matching failed
    Point3_t Pos;               // Position center using 3D-matched points on the Tjs - 3D fit
    Vector3_t Dir;              // Direction using 3D-matched points on the Tjs - 3D fit
    float AspectRatio {-1};           // Aspect ratio calculated when doing a 3D fit
  };
  
  struct PFPStruct {
    std::vector<int> TjIDs;
    std::vector<float> TjCompleteness;  // fraction of TP points that are 3D-matched
    std::vector<TrajPoint3> Tp3s;    // TrajCluster 3D trajectory points
    // Start is 0, End is 1
    std::array<Point3_t, 2> XYZ;        // XYZ position at both ends (cm)
    std::array<Vector3_t, 2> Dir;
    std::array<Vector3_t, 2> DirErr;
    std::array<std::vector<float>, 2> dEdx;
    std::array<std::vector<float>, 2> dEdxErr;
    std::array<unsigned short, 2> Vx3ID {0, 0};
    int BestPlane {-1};
    // stuff for constructing the PFParticle
    int PDGCode {-1};
    std::vector<int> DtrIDs;
    size_t ParentID {0};       // Parent PFP ID (or ID of self if no parent exists)
    geo::TPCID TPCID;
    float EffPur {0};                     ///< Efficiency * Purity
    unsigned int MCPartListIndex {UINT_MAX};
    unsigned short MatchVecIndex {USHRT_MAX};
    float CosmicScore{0};
    float AspectRatio {0};
    unsigned short ID {0};
    std::array<std::bitset<8>, 2> StopFlag {};  // Bitset that encodes the reason for stopping
    bool Primary;             // PFParticle is attached to a primary vertex
    bool NeedsUpdate {true};    // Set true if the PFParticle needs to be (re-)defined
    bool DirectionFixed {false};  // Fix the direction of the pfp if it is small angle
  };

  struct ShowerPoint {
    Point2_t Pos;       // Hit Position in the normal coordinate system
    Point2_t RotPos;    // Position rotated into the shower coordinate system (0 = along, 1 = transverse)
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
    std::vector<ShowerPoint> ShPts;    // Trajectory points inside the shower
    float Angle {0};                   // Angle of the shower axis
    float AngleErr {3};                 // Error
    float AspectRatio {1};              // The ratio of charge weighted transverse/longitudinal positions
    float DirectionFOM {1};
    std::vector<Point2_t> Envelope; // Vertices of a polygon that encompasses the shower
    float EnvelopeArea {0};
    float ChgDensity {0};                   // Charge density inside the Envelope
    float Energy {0};
    float ParentFOM {10};
    int ID {0}; 
    int ParentID {0};  // The ID of an external parent Tj that was added to the shower
    unsigned short TruParentID {0};
    unsigned short SS3ID {0};     // ID of a ShowerStruct3D to which this 2D shower is matched
    bool NeedsUpdate {true};       // This is set true whenever the shower needs to be updated
    bool Constraint3D {false};    // Some properties are defined by the 3D shower to which it is associated
  };
  
  // Shower variables filled in MakeShowers. These are in cm and radians
  struct ShowerStruct3D {
    Vector3_t Dir;              //
    Vector3_t DirErr;           // DirErr is hijacked to store the shower rms at the start, center and end sections
    Point3_t Pos;               //
    Point3_t PosErr;            // PosErr is hijacked to temporarily store the charge in the three sections
    Point3_t ChgPos;            // position of the center of charge
    Point3_t EndPos;            // end position
    double Len {1};
    double OpenAngle {0.12};
    std::vector<double> Energy;
    std::vector<double> EnergyErr;
    std::vector<double> MIPEnergy;
    std::vector<double> MIPEnergyErr;
    std::vector<double> dEdx;
    std::vector<double> dEdxErr;
    geo::TPCID TPCID;
    std::vector<unsigned short> PFPIDs;  // list of indices of InShower PFParticles 
    std::vector<unsigned short> CotIndices;  // list of indices of 2D showers in tjs.cots
    std::vector<unsigned int> Hits;
    int BestPlane;
    int ID;
    float FOM;
    unsigned short PFPIndex {USHRT_MAX};    // The index of the PFParticle for this shower
    unsigned short Vx3ID {0};
  };

  struct ShowerTreeVars {
    // run, subrun, and event are also saved to this tree

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
   
    std::vector<int> TjID;
    std::vector<int> IsShowerTj; // indicates tj is an shower trajectory
    std::vector<int> ShowerID; // shower ID associated w/ trajectory. -1 = no shower
    std::vector<int> IsShowerParent; // this tj was chosen as a parent tj
    std::vector<int> StageNum; // stage of reconstruction
    std::vector<std::string> StageName; // stage name

    // envelope information
    std::vector<float> Envelope;
    std::vector<int>EnvPlane;
    std::vector<int>EnvStage;
    std::vector<int>EnvShowerID;

    int nStages {0};
    unsigned short nPlanes {0};

  };

  struct CRTreeVars {
    std::vector<int>   cr_origin;
    std::vector<float> cr_pfpxmin;
    std::vector<float> cr_pfpxmax;
    std::vector<float> cr_pfpyzmindis;
  };

  // Algorithm modification bits
  typedef enum {
    kMaskHits,
    kMaskBadTPs,
    kMichel,
    kDeltaRay,
    kCTKink,        ///< kink found in CheckTraj
    kCTStepChk,
    kTryWithNextPass,
    kRvPrp,
    kCHMUH,
    kSplit,
    kComp3DVx,
    kComp3DVxIG,
    kHED, // High End Delta
    kHamVx,
    kHamVx2,
    kJunkVx,
    kJunkTj,
    kKilled,
    kMerge,
    kTEP,
    kCHMEH,
    kFillGap,
    kUseGhostHits,
    kMrgGhost,
    kChkInTraj,
    kStopBadFits,
    kFixBegin,
    kFTBChg,
    kBeginChg,
    kFixEnd,
    kUUH,
    kVtxTj,
    kChkVxTj,
    kMisdVxTj,
    kPhoton,
    kNoFitToVx,
    kVxMerge,
    kNoKinkChk,
    kSoftKink,
    kChkStop,
    kChkStopEP,
    kChkChgAsym,
    kFTBRvProp,
    kStopAtTj,
    kMat3D,
    kMat3DMerge,
    kSplit3DKink,
    kTjHiVx3Score,
    kVtxHitsSwap,
    kSplitHiChgHits,
    kInShower,
    kKillInShowerVx,
    kShowerTj,
    kShwrParent,
    kChkShwrParEnd,
    kKillShwrNuPFP,
    kMergeOverlap,
    kMergeSubShowers,
    kMergeNrShowers,
    kMergeShChain,
    kSplitTjCVx,
    kSetDir,
    kAlgBitSize     ///< don't mess with this line
  } AlgBit_t;
  
  // Stop flag bits
  typedef enum {
    kSignal,
    kAtKink,
    kAtVtx,
    kBragg,
    kAtTj,
    kOutFV,
    kFlagBitSize     ///< don't mess with this line
  } StopFlag_t; 
  
  typedef enum {
    kEnvNearTj,
    kEnvNearShower,
    kEnvOverlap,
    kEnvUnusedHits
  } TPEnvironment_t;
  
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
    float WirePitch;
    std::vector<float> AveHitRMS;      ///< average RMS of an isolated hit
    // The variables below do change in size from event to event
    ShowerTreeVars stv; // 
    bool SaveShowerTree;

    // Save histograms to develop cosmic removal tools
    CRTreeVars crt;
    bool SaveCRTree;
    bool TagCosmics;

    std::vector<Trajectory> allTraj; ///< vector of all trajectories in each plane
    std::vector<Tj2Pt> mallTraj;      ///< vector of trajectory points ordered by increasing X
    std::vector<TCHit> fHits;
    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector<std::vector< std::pair<int, int>>> WireHitRange;
    std::vector<std::pair<float, float>> srcHits;   ///< (lower,upper) time of each hit in the source collection
    std::vector<std::vector< std::pair<int, int>>> srcWireHitRange;
    std::vector<float> AngleRanges; ///< list of max angles for each angle range
    std::vector< ClusterStore > tcl; ///< the clusters we are creating
    std::vector< VtxStore > vtx; ///< 2D vertices
    std::vector< Vtx3Store > vtx3; ///< 3D vertices
    std::vector<MatchStruct> matchVec; ///< 3D matching vector
    std::vector<PFPStruct> pfps;
    std::vector<ShowerStruct> cots;       // Clusters of Trajectories that define 2D showers
    std::vector<ShowerStruct3D> showers;  // 3D showers
    std::vector<float> Vertex2DCuts; ///< Max position pull, max Position error rms
    std::vector<float> Vertex3DCuts;   ///< 2D vtx -> 3D vtx matching cuts 
    std::vector<float> VertexScoreWeights;
    std::vector<short> DeltaRayTag; ///< min length, min MCSMom and min separation (WSE) for a delta ray tag
    std::vector<short> MuonTag; ///< min length and min MCSMom for a muon tag
    std::vector<float> ShowerTag; ///< [min MCSMom, max separation, min # Tj < separation] for a shower tag
    std::vector<float> KinkCuts; ///< kink angle, nPts fit, (alternate) kink angle significance
    std::vector<float> Match3DCuts;  ///< 3D matching cuts
    std::vector<float> MatchTruth;     ///< Match to MC truth
    std::vector<float> ChargeCuts;
    std::vector<simb::MCParticle*> MCPartList;
    unsigned int EventsProcessed;
    unsigned int Run;
    unsigned int SubRun;
    unsigned int Event;
    std::bitset<64> UseAlg;  ///< Allow user to mask off specific algorithms
    const geo::GeometryCore* geom;
    const detinfo::DetectorProperties* detprop;
    calo::CalorimetryAlg* caloAlg;
    short StepDir;        ///< the normal user-defined stepping direction = 1 (US -> DS) or -1 (DS -> US)
    short NPtsAve;         /// number of points to find AveChg
    bool SelectEvent;     ///< select this event for use in the performance metric, writing out, etc
    bool TestBeam;      ///< Expect tracks entering from the front face. Don't create neutrino PFParticles
   };

} // namespace tca

#endif // ifndef TRAJCLUSTERALGDATASTRUCT_H
