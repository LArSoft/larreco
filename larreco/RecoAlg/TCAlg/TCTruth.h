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
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/TCHist.h"

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
      MCP_TSum = 0;
      MCP_EPTSum = 0;
      MCP_Cnt = 0;
      nBadT = 0;
      nBadP = 0;
   }

    void Initialize();
    void MatchTruth();
    void MatchTAndSum();
    void PrintResults(int eventNum) const;

    // Variables for summing Eff*Pur for electrons, muons, pions, kaons and protons for Trajectories
    std::array<short, 5> EPCnts {{0}};
    std::array<float, 5> TSums;    // sum of kinetic energy
    std::array<float, 5> EPTSums;   // E*P sum weighted by kinetic energy for 5 particle types

    float MCP_TSum;                // T sum of MCParticles that should be reconstructed in 3D
    float MCP_EPTSum;              // E*P weighted T sum of MCParticles that ARE reconstructed in 3D
    float MCP_Cnt;                 // Count of MCParticles that should be reconstructed in 3D

    float nBadT;    // Number of reconstructed Tjs with a bad EP
    float nBadP;    // Number of reconstructed pfps with a bad EP

    HistStuff hist;

  }; // TruthMatcher class
/*
  class MCParticleListUtils
  {
  public:

    MCParticleListUtils(TCSlice& my_slc);
    ShowerStruct3D MakeCheatShower(TCSlice& slc, unsigned int mcpIndex, Point3_t primVx, int& truParentPFP);
    bool PrimaryElectronStart(Point3_t& start, Vector3_t& dir, float& energy);
    int PrimaryElectronPFPID(TCSlice& slc);
    int PrimaryElectronTjID(TCSlice& slc, CTP_t inCTP);
    int MCParticleStartTjID(TCSlice& slc, unsigned int MCParticleListIndex, CTP_t inCTP);
    unsigned int GetMCPListIndex(TCSlice& slc, const TrajPoint& tp);
    unsigned int GetMCPListIndex(TCSlice& slc, const Trajectory& tj, unsigned short& nTruHits);
    unsigned int GetMCPListIndex(TCSlice& slc, const ShowerStruct& ss, unsigned short& nTruHits);

  }; // MCParticleListUtils class
*/

} // namespace tca

#endif // ifndef TRAJCLUSTERALGTRUTH_H
