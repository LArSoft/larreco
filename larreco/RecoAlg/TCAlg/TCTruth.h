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
      EPSums.fill(0.0);
      EPTSums.fill(0.0);
    }
    
    void Initialize();
    void MatchTrueHits();
    void MatchTruth(const HistStuff& hist, unsigned int& fEventsProcessed);
    void PrintResults(int eventNum) const;
    
    TjStuff& tjs;
    // Variables for summing Eff*Pur for electrons, muons, pions, kaons and protons
    std::array<short, 5> EPCnts {{0}};
    std::array<float, 5> EPSums;
    std::array<float, 5> EPTSums;
    
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
