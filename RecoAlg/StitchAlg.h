/////////////////////////////////////////////////////////////////
//  \StitchAlg.h
//  echurch@fnal.gov
////////////////////////////////////////////////////////////////////

#ifndef STITCHALG_H
#define STITCHALG_H

// C/C++ standard libraries
#include <vector>

// art libraries
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Principal/Handle.h"

// LArSoft libraries
#include "RecoBase/Track.h"

namespace fhicl { class ParameterSet; }
namespace art { class Event; }


namespace trkf{


class StitchAlg 
{


 public:
  StitchAlg (fhicl::ParameterSet const& pset)  ;
  virtual ~StitchAlg ()  {};

  void reconfigure(fhicl::ParameterSet const& pset) ;

  void FindHeadsAndTails( const art::Event& e, const std::string& t);
  void FirstStitch(const std::vector<art::PtrVector <recob::Track>>::iterator itvvArg, const std::vector <recob::Track>::iterator itvArg);
  void WalkStitch();
  bool CommonComponentStitch();

  void GetTrackComposites(std::vector <art::PtrVector <recob::Track> > & c) { c = fTrackComposite;};
  void GetTracks(std::vector <recob::Track>& t) { t = fTrackVec ;};
  art::Handle< std::vector< recob::Track > > ftListHandle;

 private:

  std::vector <std::tuple <std::string, int, int, double, double> > fh;
  std::vector <std::tuple <std::string, int, int, double, double> > ft;
  int ftNo;
  double fCosAngTol;
  double fSepTol;


  std::vector <art::PtrVector <recob::Track> > fTrackComposite;
  std::vector <recob::Track> fTrackVec;
  std::vector <std::vector <std::string>>  fHT;

};



} // namespace

#endif // ifndef STITCHALG_H
