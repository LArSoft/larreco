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


#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/TCHist.h"


namespace tca {
  
  class TruthMatcher
  {
    
    public:
    
    TruthMatcher(TjStuff& my_tjs) : tjs(my_tjs) {
      EPCnts.fill(0); 
      TSums.fill(0.0);
      EPTSums.fill(0.0);
      TruVxCounts.fill(0);
      MCP_TSum = 0;
      MCP_EPTSum = 0;
      MCP_Cnt = 0;
      PFP_CntGoodMat = 0;
      PFP_Cnt = 0;
      nBadEP = 0;
   }
    
    void Initialize();
    void MatchTrueHits(const HistStuff& hist);
    void MatchTruth(const HistStuff& hist, bool fStudyMode);
    void PrintResults(int eventNum) const;
    
    TjStuff& tjs;
    // Variables for summing Eff*Pur for electrons, muons, pions, kaons and protons
    std::array<short, 5> EPCnts {{0}};
    std::array<float, 5> TSums;    // sum of kinetic energy
    std::array<float, 5> EPTSums;   // E*P sum weighted by kinetic energy for 5 particle types
    float MCP_TSum;                   // T sum
    float MCP_EPTSum;                  // E*P weighted sum
    float MCP_Cnt;
    float PFP_CntGoodMat;              // Count of well matched MCParticles (> 1 plane)
    float PFP_Cnt;
    
    float fNeutrinoEnergy;
    float fSourceParticleEnergy; //< in MeV
    unsigned short nBadEP;
    // Counts of:
    // [0] = the number of true neutrino interaction vertices in the fiducial volume
    // [1] = the number of those [0] that are reconstructable by TrajCluster
    // [2] = the number of those [1] in which a reconstructed vertex is within 1 cm of the true vertex
    // [3] = the number of those [2] in which the reconstructed vertex is identified as the true vertex
    std::array<unsigned short, 4> TruVxCounts;
/*
    // number of reconstructable primary particles in the event
    unsigned short nTruPrimaryOK;
    // number of reconstructable neutrino vertices in ALL events
    unsigned short nTruPrimaryVtxOK;
    // number of reconstructable neutrino vertices in ALL events that were reconstructed
    unsigned short nTruPrimaryVtxReco;
*/
  };
  
  class MCParticleListUtils
  {
  public:

    MCParticleListUtils(TjStuff& my_tjs) : tjs(my_tjs) {}
    TjStuff& tjs;
    void MakeTruTrajPoint(unsigned short MCParticleListIndex, TrajPoint& tp);
    unsigned short MCParticleStartTjID(unsigned short MCParticleListIndex, CTP_t inCTP);
    unsigned short GetMCPartListIndex(const Trajectory& tj, unsigned short& nTruHits);
    unsigned short GetMCPartListIndex(const ShowerStruct& ss, unsigned short& nTruHits);

  };


} // namespace tca

#endif // ifndef TRAJCLUSTERALGTRUTH_H
