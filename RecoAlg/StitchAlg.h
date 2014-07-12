/////////////////////////////////////////////////////////////////
//  \StitchAlg.h
//  echurch@fnal.gov
////////////////////////////////////////////////////////////////////

#ifndef STITCHALG_H
#define STITCHALG_H
#include <vector>
#include <cmath>
#include <iostream>
#include <stdint.h>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

#include "Geometry/Geometry.h"
#include "RecoBase/Track.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"

#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVectorF.h"
#include "TVector.h"
#include "TH1.h"



//namespace recob { class Hit; }
//namespace recob { class Track; }

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
  void CommonComponentStitch();

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
