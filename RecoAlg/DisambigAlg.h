/////////////////////////////////////////////////////////////////
//  \fileDisambigAlg.h
//  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////
#ifndef DisambigAlg_H
#define DisambigAlg_H
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <stdint.h>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Principal/Event.h"

#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoAlg/APAGeometryAlg.h"

#include "MCCheater/BackTracker.h"

#include "TMatrixD.h"
#include "TVectorD.h"



namespace apa{



  //--------------------------------------------------------------- 
  class DisambigAlg {
  public:
    
    
    DisambigAlg(fhicl::ParameterSet const& pset);
    virtual ~DisambigAlg();
    
    void               reconfigure(fhicl::ParameterSet const& p);

    void               RunDisambig( art::Handle< std::vector<recob::Hit> > GausHits );
                                                                  ///< Run disambiguation as currently configured

    void               TrivialDisambig    ( unsigned int apa );   ///< Make the easiest and safest disambiguations in apa
    void               Crawl              ( unsigned int apa );   ///< Extend what we disambiguation we do have in apa
    unsigned int       FindChanTimeEndPts ( unsigned int apa );   ///< Basic endpoint-hit finder per apa
    void               UseEndPts          ( unsigned int apa );   ///< Try to associate endpoint hits and crawl from there
    unsigned int       CompareViews       ( unsigned int apa );   ///< Compare U and V to see if one says something about the other
    void               AssessDisambigSoFar( unsigned int apa );   ///< See how much disambiguation has been done in this apa so far


    std::map<unsigned int, double>         fUeffSoFar;
    std::map<unsigned int, double>         fVeffSoFar;
    std::map<unsigned int, unsigned int>   fnUSoFar;
    std::map<unsigned int, unsigned int>   fnVSoFar;
    std::map<unsigned int, unsigned int>   fnDUSoFar;
    std::map<unsigned int, unsigned int>   fnDVSoFar;

    std::vector< std::pair<art::Ptr<recob::Hit>, geo::WireID> > fDisambigHits;
                                                                   ///< The final list of hits to pass back to be made

    private:

    // other classes we will use
    apa::APAGeometryAlg                           fAPAGeo;
    art::ServiceHandle<geo::Geometry>             geom;
    art::ServiceHandle<util::DetectorProperties>  detprop;
    art::ServiceHandle<cheat::BackTracker> bt;                     ///< For *TEMPORARY* monitering of potential problems

    // Hits organization
    std::map< uint32_t, std::vector< art::Ptr< recob::Hit > > >        fChannelToHits; 
    std::map< unsigned int, std::vector< art::Ptr< recob::Hit > > >    fAPAToUVHits, fAPAToZHits; 
    std::map< unsigned int, std::vector< art::Ptr< recob::Hit > > >    fAPAToHits; 
                                                                   ///\ todo: Channel/APA to hits can be done in a unified way
    std::map< unsigned int, std::vector< art::Ptr< recob::Hit > > >    fAPAToEndPHits;
    std::map< unsigned int, std::vector< std::pair<art::Ptr<recob::Hit>, geo::WireID> > >  fAPAToDHits;
                                                                   ///< Hold the disambiguations per APA



    // data/function to keep track of disambiguation along the way
    std::map<std::pair<double,double>, geo::WireID>                      fChanTimeToWid;
                                    ///< If a hit is disambiguated, map its chan and peak time to the chosen wireID
    std::map< unsigned int, std::map<std::pair<double,double>, bool> >   fHasBeenDisambiged;
                                    ///< Convenient way to keep track of disambiguation so far
    void          MakeDisambigHit( art::Ptr<recob::Hit> hit, 
				   geo::WireID, 
				    unsigned int apa);
                                    ///< Makes a disambiguated hit while keeping track of what has already been disambiguated



    // Functions that support disambiguation methods
    unsigned int  MakeCloseHits        (int ext, geo::WireID wid, double Dmin, double Dmax);
                                    ///< Having disambiguated a time range on a wireID, extend to neighboring channels
    bool          HitsOverlapInTime    ( art::Ptr<recob::Hit> hitA, art::Ptr<recob::Hit> hitB );
    bool          HitsReasonablyMatch  ( art::Ptr<recob::Hit> hitA, art::Ptr<recob::Hit> hitB );
                                    ///\ todo: Write function that compares hits more detailedly



    // Configure the disambiguation
    bool         fCrawl;
    bool         fUseEndP;
    bool         fCompareViews;
    unsigned int fNChanJumps;       ///< Number of channels the crawl can jump over
    double       fCloseHitsRadius;  ///< Distance (cm) away from a hit to look when checking if it's an endpoint 
    double       fMaxEndPDegRange;  ///< Within the close hits radius, how spread can the majority 
                                    ///< of the activity be around a possible endpoint

  }; // class DisambigAlg

} // namespace apa

#endif // ifndef DisambigAlg_H
