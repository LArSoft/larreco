////////////////////////////////////////////////////////////////////////
//
// DisambigAlg.cxx
//
// tylerdalion@gmail.com
//
// description
//
//
// Since we are letting the module make the final column of hits, a disambiguated 
// hit in this algorithm is a pair< art::Ptr<recob::Hit>, geo::WireID >
//
//
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/DisambigAlg.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

// **tomporarily** here to look at performance without noise hits
#include "MCCheater/BackTracker.h"


#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>


namespace apa{

DisambigAlg::DisambigAlg(fhicl::ParameterSet const& pset) 
  : fAPAGeo(pset.get< fhicl::ParameterSet >("APAGeometryAlg"))
{
  this->reconfigure(pset); 
}

//----------------------------------------------------------
DisambigAlg::~DisambigAlg()
{
}

//----------------------------------------------------------
void DisambigAlg::reconfigure(fhicl::ParameterSet const& p)
{

  fCrawl            =  p.get< bool >("Crawl");
  fUseEndP          =  p.get< bool >("UseEndP");
  fCompareViews     =  p.get< bool >("CompareViews");
  fCloseHitsRadius  =  p.get< double >("CloseHitsRadius");
  fMaxEndPDegRange  =  p.get< double >("MaxEndPDegRange");
  fNChanJumps       =  p.get< unsigned int >("NChanJumps");

}


//----------------------------------------------------------
//----------------------------------------------------------
void DisambigAlg::RunDisambig(   art::Handle< std::vector<recob::Hit> > ChannelHits   )
{

  // **tomporarily** here to look at performance without noise hits
  art::ServiceHandle<cheat::BackTracker> bt;

  fUeffSoFar.clear();
  fVeffSoFar.clear();
  fnUSoFar.clear();
  fnVSoFar.clear();
  fnDUSoFar.clear();
  fnDVSoFar.clear();
  fChannelToHits.clear();
  fAPAToUVHits.clear();
  fAPAToZHits.clear();
  fAPAToHits.clear(); 
  fAPAToEndPHits.clear();
  fAPAToDHits.clear();
  fDisambigHits.clear();
  fChanTimeToWid.clear();
  fHasBeenDisambiged.clear();


  std::vector< art::Ptr<recob::Hit> >  ChHits;
  art::fill_ptr_vector(ChHits, ChannelHits);


  fHasBeenDisambiged.clear();
  unsigned int skipNoise(0);
  // Map hits by channel/APA, initialize the disambiguation status map
  for( size_t h = 0; h < ChHits.size(); h++ ){
    art::Ptr<recob::Hit> hit = ChHits[h];

    // **temporary** option to skip noise hits
    try{ bt->HitToXYZ(hit); }
    catch(...){ 
      skipNoise++;
      continue;
    }

    geo::View_t view = hit->View();
    unsigned int apa(0), cryo(0);
    fAPAGeo.ChannelToAPA(hit->Channel(), apa, cryo);
    fAPAToHits[apa].push_back(hit);
    if( view == geo::kZ ){ 
      fAPAToZHits[apa].push_back(hit);
      continue;
    } else if ( view==geo::kU || view==geo::kV ){
      std::pair<double,double> ChanTime( hit->Channel()*1., hit->PeakTime()*1. );
      this->fHasBeenDisambiged[apa][ChanTime] = false;
      fChannelToHits[ hit->Channel() ].push_back( hit );
      fAPAToUVHits[apa].push_back(hit);
    }
  }

  if(skipNoise>0)
    mf::LogWarning("DisambigAlg")<<"\nSkipped "<< skipNoise <<" induction noise hits using the BackTracker.\n"
				 <<"This is only to temporarily deal with the excessive amount of noise due to the bad deconvolution.\n"; 

  mf::LogVerbatim("RunDisambig")<<"\n~~~~~~~~~~~ Running Disambiguation ~~~~~~~~~~~\n";

  std::map<unsigned int, std::vector< art::Ptr< recob::Hit> > >::iterator APA_it;
  for( APA_it = fAPAToUVHits.begin(); APA_it != fAPAToUVHits.end(); APA_it++ ){
    unsigned int apa = APA_it->first;

    mf::LogVerbatim("RunDisambig")<<"APA "<<apa<<":";

    fUeffSoFar[apa] = 0.;   fVeffSoFar[apa] = 0.;
    fnUSoFar[apa]   = 0;    fnVSoFar[apa]   = 0;
    fnDUSoFar[apa]  = 0;    fnDVSoFar[apa]  = 0;

    // Always run this...
    this->TrivialDisambig(apa);     
    this->AssessDisambigSoFar(apa);
    mf::LogVerbatim("RunDisambig") << "  Trivial Disambig -->  " 
				   << fnDUSoFar[apa] << " / " << fnUSoFar[apa] << " U,  "
				   << fnDVSoFar[apa] << " / " << fnVSoFar[apa] << " V";
    
    
    // ... and pick the rest with the configurations.
    if( fCrawl ){
      this->Crawl(apa);     
      this->AssessDisambigSoFar(apa);
      mf::LogVerbatim("RunDisambig") << "  Crawl            -->  " 
				     << fnDUSoFar[apa] << " / " << fnUSoFar[apa] << " U,  "
				     << fnDVSoFar[apa] << " / " << fnVSoFar[apa] << " V";
    }
    

    if(fUseEndP){
      this->FindChanTimeEndPts(apa);
      this->UseEndPts(apa); // does the crawl from inside
      this->AssessDisambigSoFar(apa);
      mf::LogVerbatim("RunDisambig") << "  Endpoint Crawl   -->  " 
				     << fnDUSoFar[apa] << " / " << fnUSoFar[apa] << " U,  "
				     << fnDVSoFar[apa] << " / " << fnVSoFar[apa] << " V";
    }


    if(fCompareViews){
      unsigned int nDisambig(1);
      while(nDisambig > 0){
	nDisambig = 0;
	nDisambig = this->CompareViews(apa);
	this->Crawl(apa);
      }
      this->AssessDisambigSoFar(apa);
      mf::LogVerbatim("RunDisambig") << "  Compare Views    -->  " 
				     << fnDUSoFar[apa] << " / " << fnUSoFar[apa] << " U,  "
				     << fnDVSoFar[apa] << " / " << fnVSoFar[apa] << " V";
    }      


    //this->GatherLeftoverHits()

    // For now just buld a simple list to get from the module
    for(size_t i=0; i<fAPAToDHits[apa].size(); i++)
      fDisambigHits.push_back(fAPAToDHits[apa][i]); 

  } // end loop through APA


}




//-------------------------------------------------
//-------------------------------------------------
void DisambigAlg::MakeDisambigHit( art::Ptr<recob::Hit> hit, 
					  geo::WireID wid,
					  unsigned int apa )
{ 

  std::pair<double,double> ChanTime( hit->Channel()*1., hit->PeakTime()*1. );
  if( fHasBeenDisambiged[apa][ChanTime] ) return; 

  if( !wid.isValid ){ 
    mf::LogWarning("InvalidWireID") << "wid is invalid, hit not being made\n";
    return; }
    
  std::pair<art::Ptr<recob::Hit>,geo::WireID> Dhit(hit, wid);
  fAPAToDHits[apa].push_back(Dhit);
  fHasBeenDisambiged[apa][ChanTime] = true;
  fChanTimeToWid[ChanTime] = wid;
  return;
  
}



//----------------------------------------------------------
//----------------------------------------------------------
  bool DisambigAlg::HitsOverlapInTime( art::Ptr<recob::Hit> hitA, 
					      art::Ptr<recob::Hit> hitB )
{
    double AsT = hitA->StartTime();
    double AeT = hitA->EndTime();
    double BsT = hitB->StartTime();
    double BeT = hitB->EndTime();

    if( hitA->View() == geo::kU ){ AsT -= detprop->TimeOffsetU(); AeT -= detprop->TimeOffsetU(); }
    else if( hitA->View() == geo::kV ){ AsT -= detprop->TimeOffsetV(); AeT -= detprop->TimeOffsetV(); }
    else if( hitA->View() == geo::kZ ){ AsT -= detprop->TimeOffsetZ(); AeT -= detprop->TimeOffsetZ(); }
    if( hitB->View() == geo::kU ){ BsT += detprop->TimeOffsetU(); BeT -= detprop->TimeOffsetU(); }
    else if( hitB->View() == geo::kV ){ BsT -= detprop->TimeOffsetV(); BeT -= detprop->TimeOffsetV(); }
    else if( hitA->View() == geo::kZ ){ AsT -= detprop->TimeOffsetZ(); AeT -= detprop->TimeOffsetZ(); }

    return (AsT<=BsT && BsT<=AeT) || (AsT<=BeT && BeT<=AeT) || 
           (BsT<=AsT && AsT<=BeT) || (BsT<=AeT && AeT<=BeT);
}



//----------------------------------------------------------
//----------------------------------------------------------
void DisambigAlg::TrivialDisambig( unsigned int apa )
{

  // Loop through ambiguous hits (U/V) in this APA
  for( size_t h=0; h<fAPAToUVHits[apa].size(); h++ ){
    const art::Ptr<recob::Hit> hit = fAPAToUVHits[apa][h];
    uint32_t chan = hit->Channel();
    unsigned int peakT = hit->PeakTime();

    std::vector<geo::WireID> hitwids = geom->ChannelToWire(chan);
    std::vector<bool> IsReasonableWid(hitwids.size(),false);
    unsigned short nPossibleWids(0);
    for(size_t w=0; w<hitwids.size(); w++){
      geo::WireID wid = hitwids[w];

      double xyzStart[3] = {0.};  double xyzEnd[3] = {0.};
      geom->WireEndPoints(wid.Cryostat, wid.TPC, wid.Plane, wid.Wire, xyzStart, xyzEnd);
      unsigned int side(wid.TPC%2), cryo(wid.Cryostat);
      double zminPos(xyzStart[2]), zmaxPos(xyzEnd[2]);

      // get appropriate x and y with tpc center
      TVector3 tpcCenter(0,0,0);
      unsigned int tpc = 2*apa + side - cryo*geom->NTPC(); // apa number does not reset per cryo
      tpcCenter = geom->Cryostat(cryo).TPC(tpc).LocalToWorld(tpcCenter);
      
      // get channel range
      TVector3 Min(tpcCenter); Min[2] = zminPos;
      TVector3 Max(tpcCenter); Max[2] = zmaxPos;
      uint32_t ZminChan = geom->NearestChannel( Min, 2, tpc, cryo );
      uint32_t ZmaxChan = geom->NearestChannel( Max, 2, tpc, cryo );
      
      for( size_t z=0; z < fAPAToZHits[apa].size(); z++ ){
	uint32_t chan = fAPAToZHits[apa][z]->Channel();
	if( chan <= ZminChan || ZmaxChan <= chan ) continue;
	art::Ptr<recob::Hit> zhit = fAPAToZHits[apa][z];

// 	try{ bt->HitToXYZ(zhit); }
// 	catch(...){ 
// 	  mf::LogWarning("DisambigAlg")<<"skipping a noise collection hit"; 
// 	  continue;
// 	}

	if( this->HitsOverlapInTime(hit,zhit) ){
	  IsReasonableWid[w] = true; nPossibleWids++; break; }
      } 

    } // end hit chan-wid loop

      

    if(nPossibleWids==0){
      std::vector<double> xyz;
      try{ xyz = bt->HitToXYZ(hit); } // TEMPORARY
      catch(...){ continue; }
      ///\ todo: Figure out why sometimes non-noise hits dont match any Z hits at all.
      mf::LogWarning ("UniqueTimeSeg") << "U/V hit inconsistent with Z info; peak time is " 
				       << peakT << " in APA " << apa << " on channel " << hit->Channel(); 
    }
    else if(nPossibleWids==1){
      for(size_t d=0; d<hitwids.size(); d++) 
	if(IsReasonableWid[d]) this->MakeDisambigHit( hit, hitwids[d], apa );
    }
    else if(nPossibleWids==2){
      ///\ todo: Add mechanism to at least eliminate the wids that aren't even possible, for the benefit of future methods
    }
    

  } // end ambig hits loop  


}



//----------------------------------------------------------
//----------------------------------------------------------
unsigned int DisambigAlg::MakeCloseHits( int ext, geo::WireID Dwid, double Dmin, double Dmax)
{

  // Function to look, on a channel *ext* channels away from
  // a disambiguated hit channel, for hits with time windows 
  // touching range *Dmin to Dmax*. If found, make such a
  // hit to have a wireID adjacent to supplied *wid*.
  // Returns number of NEW hits made.



  uint32_t Dchan = geom->PlaneWireToChannel(Dwid.Plane, Dwid.Wire, Dwid.TPC, Dwid.Cryostat);
  geo::View_t view = geom->View(Dchan);
  if(view==geo::kZ) 
    throw cet::exception("MakeCloseHits") << "Function not meant for non-wrapped channels.\n";

  // Account for wrapping
  uint32_t firstChan = fAPAGeo.FirstChannelInView(view, Dchan);
  unsigned int ChanPerView = fAPAGeo.ChannelsInView(view);
  int tempchan = Dchan + ext; // need sign for the set of channels starting with channel 0
  if( tempchan < (int)firstChan ) tempchan += ChanPerView;
  if( tempchan > (int)(firstChan+ChanPerView-1) ) tempchan -= ChanPerView;
  uint32_t chan = (uint32_t)(tempchan);

  // There may just be no hits
  if( fChannelToHits.count(chan) == 0 ) return 0;

  // There are close channel hits, so for each
  unsigned int apa(0), cryo(0);
  fAPAGeo.ChannelToAPA(chan, apa, cryo);
  unsigned int MakeCount(0);
  for(size_t i=0; i<fChannelToHits[chan].size(); i++){
    art::Ptr< recob::Hit > closeHit = fChannelToHits[chan][i];
    double st = closeHit->StartTime();
    double et = closeHit->EndTime();
    std::vector<geo::WireID> wids = geom->ChannelToWire(chan);

    if( !(Dmin <= st && st <= Dmax) && !(Dmin <= et && et <= Dmax) ) continue;


    // Found hit with window overlapping given range, 
    // now find the only reasonable wireID.
    for(size_t w=0; w<wids.size(); w++){ 
      if( wids[w].TPC != Dwid.TPC ) continue;
      if( (int)(wids[w].Wire)-(int)(Dwid.Wire) != ext ) continue;    

      // In this case, we have a unique wireID.
      // Check to see if it has already been made - if so, do not incriment count
      std::pair<double,double> ChanTime( closeHit->Channel()*1., closeHit->PeakTime()*1. );
      if( !fHasBeenDisambiged[apa][ChanTime] ){ 
	this->MakeDisambigHit(closeHit, wids[w], apa); 
	MakeCount++;
	//std::cout << "     Close hit found on channel " << chan << ", time " << st<<"-"<<et << "... \n";
	//std::cout << " ... giving it wireID ("<< Dwid.Cryostat <<"," << Dwid.TPC 
	//	  << "," << Dwid.Plane << "," << Dwid.Wire << ")\n";
      }
      break;
    } // end find right wireID
    

  } // end loop through all hits on chan
  
  return MakeCount;

}






//----------------------------------------------------------
//----------------------------------------------------------
void DisambigAlg::Crawl( unsigned int apa )
{

  std::vector<art::Ptr<recob::Hit> > hits = fAPAToUVHits[apa];

  // repeat this method until stable  
  unsigned int nExtended(1);
  while( nExtended > 0 ){
    nExtended = 0;
    
    // Look for any disambiguated hit ... 
    for(size_t h=0; h < hits.size(); h++){
      std::pair<double,double> ChanTime( hits[h]->Channel()*1., hits[h]->PeakTime()*1. );
      if( !fHasBeenDisambiged[apa][ChanTime] ) continue;
      double stD = hits[h]->StartTime();
      double etD = hits[h]->EndTime();
      double hitWindow = hits[h]->EndTime() - hits[h]->StartTime();
      geo::WireID Dwid = fChanTimeToWid[ChanTime];
            
      // ... and if any neighboring-channel hits are close enough in time,
      // extend the disambiguation to the neighboring wire.      
      unsigned int extensions = 0;
      for( unsigned int ext=1; ext < fNChanJumps+1; ext++ ){
	///\ todo: Evaluate how aggressive we can be here. How far should we jump? In what cases should we quit out?
	unsigned int N(0);
	double timeExt = hitWindow*ext; 
	N += this->MakeCloseHits( (int)(-ext), Dwid, stD-5-timeExt, etD+5+timeExt);
	N += this->MakeCloseHits(  (int)(ext), Dwid, stD-5-timeExt, etD+5+timeExt);
	extensions += N;
      }  
      nExtended += extensions;
      
    } // end UV hits loop
    
  } // end while still disambiguating***
  
  // *** nested while loops allow disambigauation to hop the channel-wrap-boundary   


}



//----------------------------------------------------------
//----------------------------------------------------------
unsigned int DisambigAlg::FindChanTimeEndPts( unsigned int apa )
{

  ///\ todo: Clean up and break down into two functions.
  ///\ todo: Make the conditions more robust to some spotty hits around a potential endpoint.

  double pi = 3.14159265;
  double fMaxEndPRadRange = fMaxEndPDegRange/180. * (2*pi);
  
  for(size_t h=0; h<fAPAToHits[apa].size(); h++){
    art::Ptr<recob::Hit> centhit = fAPAToHits[apa][h];
    geo::View_t view = centhit->View();
    unsigned int plane = 0; if(view==geo::kV){ plane = 1; } else if(view==geo::kZ) plane = 2;
    std::vector<double> ChanTimeCenter(2, 0.);
    unsigned int relchan = centhit->Channel() - fAPAGeo.FirstChannelInView(centhit->Channel());
    ChanTimeCenter[0] = relchan*geom->WirePitch(view);
    ChanTimeCenter[1] = detprop->ConvertTicksToX( centhit->PeakTime(), 
						  plane, 
						  apa*2,  // tpc doesnt matter
						  centhit->WireID().Cryostat ); 
    //std::vector< art::Ptr<recob::Hit> > CloseHits;
    std::vector<std::vector<double> > CloseHitsChanTime;
    std::vector<double>  FurthestCloseChanTime(2,0.); //double maxDist = 0.;
    std::vector<double>  ClosestChanTime(2,0.);       double minDist = fCloseHitsRadius+1.;
    double ChanDistRange = fAPAGeo.ChannelsInView(view)*geom->WirePitch(view);

    for(size_t c=0; c<fAPAToHits[apa].size(); c++){
      art::Ptr<recob::Hit> closehit = fAPAToHits[apa][c];
      if(view!=closehit->View()) continue;
      if(view==geo::kZ && centhit->WireID().TPC != closehit->WireID().TPC ) continue;
      unsigned int plane = 0;  if(view==geo::kV){ plane = 1; } else if(view==geo::kZ) plane = 2;
      std::vector<double> ChanTimeClose(2, 0.);
      unsigned int relchanclose = closehit->Channel() - fAPAGeo.FirstChannelInView(closehit->Channel());
      ChanTimeClose[0] = relchanclose*geom->WirePitch(view);
      ChanTimeClose[1] = detprop->ConvertTicksToX( closehit->PeakTime(), 
						   plane, 
						   apa*2,  // tpc doesnt matter
						   closehit->WireID().Cryostat ); 
      if(ChanTimeClose == ChanTimeCenter) continue; // move on if the same one

      double ChanDist = ChanTimeClose[0]-ChanTimeCenter[0];
      if( ChanDist > ChanDistRange/2 ) ChanDist = ChanDistRange - ChanDist;
      
      double distance = std::sqrt(   std::pow(ChanDist,2) 
				   + std::pow(ChanTimeClose[1]-ChanTimeCenter[1],2)  );

      if( distance <= fCloseHitsRadius ) CloseHitsChanTime.push_back(ChanTimeClose);

      if( distance < minDist ){ 
	ClosestChanTime = ChanTimeClose;
	minDist = distance;
      }

    } // end close-by hit loop

    if(CloseHitsChanTime.size()<5) continue; // quick fix, to-be improved

    double minRad(2*pi+1.), maxRad(0.);
    bool CloseToNegPi(false), CloseToPosPi(false);
    for(size_t i=0; i<CloseHitsChanTime.size(); i++){
      std::vector<double> ThisChanTime(CloseHitsChanTime[i]);
      //if(ThisChanTime==ClosestChanTime) continue;
      double ChanDist = ThisChanTime[0]-ChanTimeCenter[0];
      if( ChanDist > ChanDistRange/2 ) ChanDist = ChanDistRange - ChanDist;
      double hitrad = std::atan2( ThisChanTime[1]-ChanTimeCenter[1] , 
				  ChanDist );
       if( hitrad > maxRad ) 	maxRad = hitrad;
       if( hitrad < minRad ) 	minRad = hitrad;
       if( hitrad + fMaxEndPRadRange > pi  )       CloseToPosPi = true;
       else if( hitrad - fMaxEndPRadRange < -pi )  CloseToNegPi = true;
    }

    // activity at this boundary automatically kills the test, move boundary and redo
    if(CloseToPosPi&&CloseToNegPi){
      for(size_t i=0; i<CloseHitsChanTime.size(); i++){
	std::vector<double> ThisChanTime(CloseHitsChanTime[i]);
	//if(ThisChanTime==ClosestChanTime) continue;
	double ChanDist = ThisChanTime[0]-ChanTimeCenter[0];
	if( ChanDist > ChanDistRange/2 ) ChanDist = ChanDistRange - ChanDist;
	double hitrad = std::atan2( ThisChanTime[1]-ChanTimeCenter[1] , 
				    ChanDist );
	if( hitrad > 0 ) hitrad = pi - hitrad; // reflec across pi/2 line
	if( hitrad < 0 ) hitrad = -pi - hitrad;
	if( hitrad > maxRad ) 	maxRad = hitrad;
	if( hitrad < minRad ) 	minRad = hitrad;
      }
    }

    if( maxRad - minRad < fMaxEndPRadRange ) fAPAToEndPHits[apa].push_back( centhit );

  } // end UV hit loop

  if(fAPAToEndPHits[apa].size()==0) return 0;
  mf::LogVerbatim("FindChanTimeEndPts") << "          Found " << fAPAToEndPHits[apa].size() 
					<< " endpoint hits in apa " << apa << std::endl;
  for(size_t ep=0; ep<fAPAToEndPHits[apa].size(); ep++){
    art::Ptr<recob::Hit> epHit = fAPAToEndPHits[apa][ep];
    mf::LogVerbatim("FindChanTimeEndPts") << "           endP on channel " << epHit->Channel()
					  << " at time " << epHit->PeakTime() << std::endl;
  }

  return fAPAToEndPHits[apa].size();

}





//----------------------------------------------------------
//----------------------------------------------------------
void DisambigAlg::UseEndPts( unsigned int apa )
{

  ///\ todo: This function could be made much cleaner and more compact

  if(fAPAToEndPHits[apa].size()==0){
    mf::LogVerbatim("UseEndPts") << "          APA " << apa << " has no endpoints.";
    return; }
  std::vector< art::Ptr<recob::Hit> > endPts = fAPAToEndPHits[apa];


  std::vector<std::vector< art::Ptr<recob::Hit> > > EndPMatch;
  unsigned short nZendPts(0);
  for(size_t ep=0; ep<endPts.size(); ep++){
    if(endPts[ep]->View()!=geo::kZ) continue;
    art::Ptr<recob::Hit> Zhit = endPts[ep];
    art::Ptr<recob::Hit> Uhit(Zhit), Vhit(Zhit);
    unsigned short Umatch(0), Vmatch(0);
    nZendPts++;
    
    // look for U and V hits overlapping in time
    for(size_t p=0; p<endPts.size(); p++){
      if(endPts[p]->View()==geo::kU){
	if( this->HitsOverlapInTime(Zhit, endPts[p]) ){
	  Uhit = endPts[p];
	  //endPts.erase(endPts.begin()+p);
	  Umatch++; }
      } else if(endPts[p]->View()==geo::kV){
	if( this->HitsOverlapInTime(Zhit, endPts[p]) ){
	  Vhit = endPts[p];
	  //endPts.erase(endPts.begin()+p);
	  Vmatch++; }
      }
    }


    TVector3 tpcCenter(0,0,0);
    unsigned int tpc(endPts[ep]->WireID().TPC), cryo(endPts[ep]->WireID().Cryostat);  
    tpcCenter = geom->Cryostat(cryo).TPC(tpc).LocalToWorld(tpcCenter);


    if(Umatch == 1 && Vmatch == 1){

      std::vector<double> yzEndPt 
	= fAPAGeo.ThreeChanPos(Uhit->Channel(), Vhit->Channel(), endPts[ep]->Channel());
      double intersect[3] = {tpcCenter[0], yzEndPt[0], yzEndPt[1]};

      geo::WireID Uwid = fAPAGeo.NearestWireIDOnChan(intersect, Uhit->Channel(), 0, tpc, cryo);
      geo::WireID Vwid = fAPAGeo.NearestWireIDOnChan(intersect, Vhit->Channel(), 1, tpc, cryo);
      this->MakeDisambigHit(Uhit, Uwid, apa);
      this->MakeDisambigHit(Vhit, Vwid, apa);
	
    } else if(Umatch == 1 && Vmatch != 1){

      std::vector< geo::WireIDIntersection > widIntersects;
      fAPAGeo.APAChannelsIntersect(Uhit->Channel(), Zhit->Channel(), widIntersects);
      if( widIntersects.size() == 0 ) continue;
      else if( widIntersects.size() == 1 ){
	double intersect[3] = {tpcCenter[0], widIntersects[0].y, widIntersects[0].z};
	geo::WireID Uwid = fAPAGeo.NearestWireIDOnChan(intersect, Uhit->Channel(), 0, tpc, cryo);
	this->MakeDisambigHit(Uhit, Uwid, apa);
      } else {
	for(size_t i=0; i < widIntersects.size(); i++){
	  // compare to V hit times, see if only one makes sense
	}
      }

    } else if(Umatch == 1 && Vmatch != 1){

      std::vector< geo::WireIDIntersection > widIntersects;
      fAPAGeo.APAChannelsIntersect(Vhit->Channel(), Zhit->Channel(), widIntersects);
      if( widIntersects.size() == 0 ) continue;
      else if( widIntersects.size() == 1 ){
	double intersect[3] = {tpcCenter[0], widIntersects[0].y, widIntersects[0].z};
	geo::WireID Vwid = fAPAGeo.NearestWireIDOnChan(intersect, Vhit->Channel(), 0, tpc, cryo);
	this->MakeDisambigHit(Vhit, Vwid, apa);
      } else {
	for(size_t i=0; i < widIntersects.size(); i++){
	  // compare to V hit times, see if only one makes sense
	}
      }

    }
    
    //endPts.erase(endPts.begin()+ep);
  }


  if( nZendPts==0 && endPts.size()==2 && this->HitsOverlapInTime(endPts[0],endPts[1]) ){
    std::vector< geo::WireIDIntersection > widIntersects;
    fAPAGeo.APAChannelsIntersect(endPts[0]->Channel(), endPts[1]->Channel(), widIntersects);
    if( widIntersects.size() == 1 ){
      TVector3 tpcCenter(0,0,0);
      unsigned int cryo = endPts[0]->WireID().Cryostat;
      unsigned int tpc = widIntersects[0].TPC;
      tpcCenter = geom->Cryostat(cryo).TPC(tpc).LocalToWorld(tpcCenter);
      double intersect[3] = {tpcCenter[0], widIntersects[0].y, widIntersects[0].z};
      unsigned int plane0(0), plane1(0);
      if( endPts[0]->View()==geo::kV ) plane0=1;
      if( endPts[1]->View()==geo::kV ) plane1=1;
      geo::WireID wid0 = fAPAGeo.NearestWireIDOnChan(intersect, endPts[0]->Channel(), plane0, tpc, cryo);
      this->MakeDisambigHit(endPts[0], wid0, apa);
      geo::WireID wid1 = fAPAGeo.NearestWireIDOnChan(intersect, endPts[1]->Channel(), plane1, tpc, cryo);
      this->MakeDisambigHit(endPts[1], wid1, apa);
    } else if( widIntersects.size() > 1 ){
	for(size_t i=0; i < widIntersects.size(); i++){
	  // compare to V hit times, see if only one makes sense
	}
      }
  }

 
  this->Crawl(apa);
 
}


//----------------------------------------------------------
//----------------------------------------------------------
void DisambigAlg::AssessDisambigSoFar( unsigned int apa )
{

  unsigned int nU(0), nV(0);
  for(size_t h=0; h < fAPAToUVHits[apa].size(); h++){
    art::Ptr<recob::Hit> hit = fAPAToUVHits[apa][h];
    if(hit->View()==geo::kU) nU++;
    else if(hit->View()==geo::kV) nV++;
  }

  unsigned int nDU(0), nDV(0);
  for(size_t h=0; h < fAPAToDHits[apa].size(); h++){
    art::Ptr<recob::Hit> hit = fAPAToDHits[apa][h].first;
    if(hit->View()==geo::kU) nDU++;
    else if(hit->View()==geo::kV) nDV++;
  }

  fUeffSoFar[apa] = (nDU*1.)/(nU*1.);
  fVeffSoFar[apa] = (nDV*1.)/(nV*1.);
  fnUSoFar[apa] = nU;
  fnVSoFar[apa] = nV;
  fnDUSoFar[apa] = nDU;
  fnDVSoFar[apa] = nDV ;


}


//----------------------------------------------------------
//----------------------------------------------------------
unsigned int DisambigAlg::CompareViews( unsigned int apa )
{

  unsigned int nDisambiguations(0);

  // loop through all hits that are still ambiguous
  for(size_t h=0; h < fAPAToUVHits[apa].size(); h++){
    art::Ptr<recob::Hit>      ambighit  = fAPAToUVHits[apa][h];
    uint32_t                  ambigchan = ambighit->Channel();
    std::pair<double,double>  ambigChanTime(ambigchan*1.,ambighit->PeakTime());
    if( fHasBeenDisambiged[apa][ambigChanTime] ) continue;
    geo::View_t               view      = ambighit->View();
    std::vector<geo::WireID>  ambigwids = geom->ChannelToWire(ambigchan);
    std::vector<unsigned int> widDcounts  (ambigwids.size(), 0);
    std::vector<unsigned int> widAcounts  (ambigwids.size(), 0);



    // loop through hits in the other view which are close in time
    for(size_t i=0; i < fAPAToUVHits[apa].size(); i++){
      art::Ptr<recob::Hit> hit = fAPAToUVHits[apa][i];
      if(hit->View()==view || !this->HitsOverlapInTime(ambighit, hit)) continue;

      // An other-view-hit overlaps in time, see what 
      // wids of the ambiguous hit's channels it overlaps
      uint32_t                  chan = hit->Channel();
      std::vector<geo::WireID>  wids = geom->ChannelToWire(chan);
      std::pair<double,double>  ChanTime(chan*1.,hit->PeakTime());
      geo::WireIDIntersection   widIntersect; // only so we can use the function
      if( fHasBeenDisambiged[apa][ChanTime] ){
	for(size_t a=0; a<ambigwids.size(); a++)
	  if( ambigwids[a].TPC == fChanTimeToWid[ChanTime].TPC  &&
	      geom->WireIDsIntersect(ambigwids[a], fChanTimeToWid[ChanTime], widIntersect) ) widDcounts[a]++;
      } else {
	// still might be able to glean disambiguation 
	// from the ambiguous hits at this time
      	for(size_t a=0; a<ambigwids.size(); a++)
	  for(size_t w=0; w<wids.size(); w++)
	    if( ambigwids[a].TPC == wids[w].TPC  &&
		geom->WireIDsIntersect(ambigwids[a], wids[w], widIntersect) ) widAcounts[a]++;
      }
    } // end loop through close-time hits


    // For now, just make a hit if either ambig or disambig hits 
    // unanimously intersect a single wireID
    unsigned int Dcount(0), Acount(0);
    for(size_t d=0; d<widDcounts.size(); d++) Dcount += widDcounts[d];
    for(size_t a=0; a<widAcounts.size(); a++) Acount += widAcounts[a];
    for(size_t d=0; d<widDcounts.size(); d++){
      if( Dcount == widDcounts[d] && Dcount>0 && Acount==0 ){ 
	this->MakeDisambigHit(ambighit, ambigwids[d], apa);
	nDisambiguations++;  } }
    for(size_t a=0; a<widAcounts.size(); a++){
      if( Acount == widAcounts[a] && Acount==1 ){
	//this->MakeDisambigHit(ambighit, ambigwids[a], apa);
	//nDisambiguations++;
      } }


  } // end loop through still ambiguous hits


  return nDisambiguations;
}


} //end namespace apa
