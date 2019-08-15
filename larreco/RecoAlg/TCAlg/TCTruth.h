////////////////////////////////////////////////////////////////////////
//
//
// TCAlg MC Truth
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGTRUTH_H
#define TRAJCLUSTERALGTRUTH_H

#include <array>
#include <vector>

#include "larreco/RecoAlg/TCAlg/DataStructs.h"
namespace geo { struct TPCID; }
namespace tca { struct HistStuff; }

namespace tca {

  class TruthMatcher
  {

    public:

    TruthMatcher() {
      EPCnts.fill(0);
      TSums.fill(0.0);
      EPTSums.fill(0.0);
      TruVxCounts.fill(0);
      MCP_TSum = 0;
      MCP_EPTSum = 0;
      MCP_Cnt = 0;
      MCP_PFP_Cnt = 0;
      Prim_TSum = 0;
      Prim_EPTSum = 0;
      PFP_Cnt = 0;
      nBadEP = 0;
      nLongInPln = 0;
      nLongMCP = 0;
      nGoodLongMCP = 0;
   }

    void Initialize();
    void MatchTruth(std::vector<unsigned int>& tpcList);
    void MatchAndSum(std::vector<unsigned int>& tpcList);
    void PrintResults(int eventNum) const;
    bool CanReconstruct(std::vector<unsigned int> mcpHits, unsigned short nDimensions, const geo::TPCID& inTPCID);
    void StudyShowerParents(TCSlice& slc, HistStuff& hist);
    void StudyElectrons(TCSlice& slc, const HistStuff& hist);
    void StudyPiZeros(TCSlice& slc, const HistStuff& hist);

    // Variables for summing Eff*Pur for electrons, muons, pions, kaons and protons for Trajectories
    std::array<short, 5> EPCnts {{0}};
    std::array<float, 5> TSums;    // sum of kinetic energy
    std::array<float, 5> EPTSums;   // E*P sum weighted by kinetic energy for 5 particle types

    float MCP_TSum;                // T sum of MCParticles that should be reconstructed in 3D
    float MCP_EPTSum;              // E*P weighted T sum of MCParticles that ARE reconstructed in 3D
    float MCP_Cnt;                 // Count of MCParticles that should be reconstructed in 3D
    float MCP_PFP_Cnt;             // Count of MCParticles that are matched to a PFParticle
    float Prim_TSum;               // T sum of Primary MCParticles that should be reconstructed in 3D
    float Prim_EPTSum;             // E*P weighted T sum of primary MCParticles that ARE reconstructed in 3D
    float PFP_Cnt;                 // Count of ALL PFParticles

    unsigned short nBadEP;      // Number of MCParticles that have >= MatchTruth[3] hits in a plane  && EP < MatchTruth[2]
    unsigned short nLongInPln;   // Number of MCParticles that have >= MatchTruth[3] hits in a plane
    unsigned short nLongMCP;   // Number of MCParticles that have >= MatchTruth[3] hits
    unsigned short nGoodLongMCP;   // Number of MCParticles that have >= 2 * MatchTruth[3] hits with EP > 0.8

    // Counts of:
    // [0] = the number of true vertices that are reconstructable (> 0 primary MCParticles)
    // [1] = [0] + a vertex was reconstructed within 1 cm of the true vertex position
    // [2] = [1] + the vertex is attached to a neutrino PFParticle
    std::array<unsigned short, 3> TruVxCounts;
  }; // TruthMatcher class

  class MCParticleListUtils
  {
  public:

    MCParticleListUtils(TCSlice& my_slc);
    void MakeTruTrajPoint(TCSlice& slc, unsigned int MCParticleListIndex, TrajPoint& tp);
    ShowerStruct3D MakeCheatShower(TCSlice& slc, unsigned int mcpIndex, Point3_t primVx, int& truParentPFP);
    bool PrimaryElectronStart(Point3_t& start, Vector3_t& dir, float& energy);
    int PrimaryElectronPFPID(TCSlice& slc);
    int PrimaryElectronTjID(TCSlice& slc, CTP_t inCTP);
    int MCParticleStartTjID(TCSlice& slc, unsigned int MCParticleListIndex, CTP_t inCTP);
    unsigned int GetMCPListIndex(TCSlice& slc, const TrajPoint& tp);
    unsigned int GetMCPListIndex(TCSlice& slc, const Trajectory& tj, unsigned short& nTruHits);
    unsigned int GetMCPListIndex(TCSlice& slc, const ShowerStruct& ss, unsigned short& nTruHits);

  }; // MCParticleListUtils class


} // namespace tca

#endif // ifndef TRAJCLUSTERALGTRUTH_H
