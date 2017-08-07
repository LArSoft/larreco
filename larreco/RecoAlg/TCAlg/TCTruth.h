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
      MCP_TSum = 0;
      MCP_EPTSum = 0;
      MCP_Cnt = 0;
      PFP_CntGoodMat = 0;
      PFP_Cnt = 0;
   }
    
    void Initialize();
    void MatchTrueHits();
    void MatchTruth(const HistStuff& hist, unsigned int fEventsProcessed);
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
    
    // number of primary particles in the event
    unsigned short nTruPrimary;
    float fNeutrinoEnergy;
    float fSourceParticleEnergy; //< in MeV
    // number of reconstructable primary particles in the event
    unsigned short nTruPrimaryOK;
    // number of reconstructable neutrino vertices in ALL events
    unsigned short nTruPrimaryVtxOK;
    // number of reconstructable neutrino vertices in ALL events that were reconstructed
    unsigned short nTruPrimaryVtxReco;
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
