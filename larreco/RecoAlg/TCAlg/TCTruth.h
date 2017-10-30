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
      RecoVx2Count = 0;
      MCP_TSum = 0;
      MCP_EPTSum = 0;
      MCP_Cnt = 0;
      PFP_CntMat = 0;
      PFP_Cnt = 0;
      nBadEP = 0;
   }
    
    void Initialize();
    void MatchTrueHits(const HistStuff& hist);
    void MatchTruth(const HistStuff& hist, bool fStudyMode);
    void AccumulatePFPSums();
    void PrintResults(int eventNum) const;
    bool CanReconstruct(unsigned short mcpIndex, unsigned short nDimensions);
    // Put hits matched to a MCParticle in CTP into a vector
    std::vector<unsigned int> PutMCPHitsInVector(unsigned short mcpIndex, CTP_t inCTP);
    
    TjStuff& tjs;
    // Variables for summing Eff*Pur for electrons, muons, pions, kaons and protons for Trajectories
    std::array<short, 5> EPCnts {{0}};
    std::array<float, 5> TSums;    // sum of kinetic energy
    std::array<float, 5> EPTSums;   // E*P sum weighted by kinetic energy for 5 particle types
    
    float MCP_TSum;                // T sum of MCParticles that should be reconstructed in 3D
    float MCP_EPTSum;              // E*P weighted T sum of MCParticles that ARE reconstructed in 3D
    float MCP_Cnt;                 // Count of MCParticles that should be reconstructed in 3D
    float PFP_CntMat;          // Count of PFParticles that are matched to an MCParticle
    float PFP_Cnt;                 // Count of ALL PFParticles
    
    float fNeutrinoEnergy;
    float fSourceParticleEnergy; //< in MeV
    unsigned short nBadEP;
    // Counts of:
    // [0] = the number of true neutrino interaction vertices in the fiducial volume
    // [1] = the number of those [0] that are reconstructable by TrajCluster
    // [2] = the number of those [1] in which a reconstructed vertex is within 1 cm of the true vertex
    // [3] = the number of those [2] in which the reconstructed vertex is identified as the true vertex
    std::array<unsigned short, 4> TruVxCounts;
    unsigned int RecoVx2Count;
  }; // TruthMatcher class
  
  class MCParticleListUtils
  {
  public:

    MCParticleListUtils(TjStuff& my_tjs) : tjs(my_tjs) {}
    TjStuff& tjs;
    void MakeTruTrajPoint(unsigned short MCParticleListIndex, TrajPoint& tp);
    unsigned short MCParticleStartTjID(unsigned short MCParticleListIndex, CTP_t inCTP);
    unsigned short GetMCPartListIndex(const Trajectory& tj, unsigned short& nTruHits);
    unsigned short GetMCPartListIndex(const ShowerStruct& ss, unsigned short& nTruHits);

  }; // MCParticleListUtils class


} // namespace tca

#endif // ifndef TRAJCLUSTERALGTRUTH_H
