////////////////////////////////////////////////////////////////////////
//
// StitchAlg.cxx
//
// echurch@fnal.gov
//
// This alg is called by TrackStich_module.cc in order to stitch together
// tracks which point at each other and have gaps to within user-set
// tolerances.
////////////////////////////////////////////////////////////////////////


#include "larreco/RecoAlg/StitchAlg.h"

// C/C++ standard libraries
#include <cmath>
#include <cstdlib>
#include <vector>

//Framework includes:
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

trkf::StitchAlg::StitchAlg(fhicl::ParameterSet const& pset)
{
  ftNo = 0;
  ftListHandle.clear();
  this->reconfigure(pset);
}

//----------------------------------------------------------

//----------------------------------------------------------
void trkf::StitchAlg::reconfigure(fhicl::ParameterSet const& pset)
{

  fCosAngTol           = pset.get< double >("CosAngTolerance", 0.95);
  fSepTol              = pset.get< double >("SpptSepTolerance", 10.0); //cm
}


void trkf::StitchAlg::FindHeadsAndTails( const art::Event& EvtArg, const std::string& trackModuleLabelArg  )
{

 fTrackVec.clear();
 fTrackComposite.clear();
 fHT.clear();

  EvtArg.getByLabel(trackModuleLabelArg,ftListHandle);

    // An element of fh and ft for each outer track. Keep the cos and sep parameters of the match and a string that indicates whether it's the second track's head or tail that gives the match, along with ii, jj, the indices of the outer and inner tracks.
    ft.clear();
    fh.clear();

    int ntrack = ftListHandle->size();
    //    std::cout << "StitchAlg.FindHeadsAndTails: Number of tracks in " << ntrack << std::endl;
    for(int ii = 0; ii < ntrack; ++ii) {
      art::Ptr<recob::Track> ptrack1(ftListHandle, ii);
      const recob::Track& track1 = *ptrack1;
      const TVector3 start1(track1.Vertex<TVector3>());
      const TVector3 end1(track1.End<TVector3>());
      const TVector3 start1Dir(track1.VertexDirection<TVector3>());
      const TVector3 end1Dir(track1.EndDirection<TVector3>());
    // For each outer track, make a vector of 1 candidate track. Doesn't need to be a vector except for desire to have a 2-iteration history.
      std::vector< std::tuple< std::string, int, int, double, double> > headvv;
      std::vector< std::tuple< std::string, int, int, double, double> > tailvv;

      // For head/tail keep a vector of candidate (cos,sep)
      std::vector< std::vector<std::pair< double, double>> > matchhead;
      std::vector< std::vector<std::pair< double, double>> > matchtail;
      // Neither end of track1 is yet matched:
      bool head(false);
      bool tail(false);

      for(int jj = ii+1; jj < ntrack; ++jj) {
	art::Ptr<recob::Track> ptrack2(ftListHandle, jj);
	const recob::Track& track2 = *ptrack2;
	const TVector3& start2(track2.Vertex<TVector3>());
	const TVector3& end2(track2.End<TVector3>());
	const TVector3& start2Dir(track2.VertexDirection<TVector3>());
	const TVector3& end2Dir(track2.EndDirection<TVector3>());
	std::string sHT2("NA"); // track2 (receptor track) H or T is tagged as matched


	bool c12((std::abs(start1Dir.Dot(end2Dir))>fCosAngTol) && ((start1-end2).Mag()<fSepTol));
	bool c21((std::abs(end1Dir.Dot(start2Dir))>fCosAngTol) && ((start2-end1).Mag()<fSepTol));
	bool c11((std::abs(start1Dir.Dot(start2Dir))>fCosAngTol) && ((start1-start2).Mag()<fSepTol));
	bool c22((std::abs(end1Dir.Dot(end2Dir))>fCosAngTol) &&  ((end1-end2).Mag()<fSepTol));

	if ( c12 || c21 || c11 || c22 )
	  {

	    sHT2 = "NA";
	    if (c12||c11) { head=true; }
            if (c11) { sHT2 = "H"; } else if (c12) { sHT2 = "T"; }
	    if (c21||c22) { tail=true; }
            if (c21) { sHT2 = "H"; } else if (c22) { sHT2 = "T"; }

	    if (head && tail) // split the tie by distance
	      {
		head = false; tail = false;
		if ( ((start1-end2).Mag() < (start2-end1).Mag()) || ((start1-end2).Mag() < (start2-end2).Mag()) || ((start1-start2).Mag() < (start2-end1).Mag()) || ((start1-start2).Mag() < (start2-end2).Mag()) )
		  { head = true; tail = false; }
		else
		  { head = false; tail = true; }
	      }

	    if (head)
	      {
		// 2-deep vector, for head and tail of 2nd track
		std::vector< std::pair <double,double> > headv;
		headv.push_back(std::pair<double,double>(abs(start1Dir.Dot(start2Dir)), (start1-start2).Mag()) );
		headv.push_back(std::pair<double,double>(abs(start1Dir.Dot(end2Dir)), (start1-end2).Mag()) );

		matchhead.push_back( headv );
		// if inferior, drop the new elements; if superior, replace the old
		if ( ((matchhead.size() > 1) && ( (matchhead.back().at(0).second < matchhead.front().at(0).second) || (matchhead.back().at(1).second < matchhead.front().at(1).second ) ) ) || matchhead.size()==1 )
		  {
		    if (matchhead.size()>1) matchhead.erase(matchhead.begin());
		    if (headvv.size()>1) headvv.erase(headvv.begin());
		    if (!sHT2.compare("H"))
		      {
			auto tupTmp =std::make_tuple(sHT2,ii,jj,matchhead.back().at(0).first,matchhead.back().at(0).second);
			headvv.push_back (tupTmp);
		      }
		    else
		      {
			auto tupTmp = std::make_tuple(sHT2,ii,jj,matchhead.back().at(1).first,matchhead.back().at(1).second);
			headvv.push_back (tupTmp);
		      }
		  }
		else
		    matchhead.pop_back();



		//	        std::cout << "TrackStitcher: satisfied head. " << std::endl;
	      } // head

	    else if (tail) // else if to prevent same stitching candidate being
		               // allowed at both head and tail
	      {
		// 2-deep vector, for head and tail of 2nd track
		std::vector< std::pair <double,double> > tailv;
		tailv.push_back(std::pair<double,double>( abs(end1Dir.Dot(start2Dir)),(start2-end1).Mag() ) );
		tailv.push_back(std::pair<double,double>( abs(end1Dir.Dot(end2Dir)),(end1-end2).Mag() ) );
		matchtail.push_back( tailv );
		// if inferior, drop the new elements; if superior, replace the old
		if ( ((matchtail.size() > 1) && ( (matchtail.back().at(0).second < matchtail.front().at(0).second) || (matchtail.back().at(1).second < matchtail.front().at(1).second ) ) ) || matchtail.size()==1 )
		  {
		    if (matchtail.size()>1) matchtail.erase(matchtail.begin());
		    if (tailvv.size()>1) tailvv.erase(tailvv.begin());
		    if (!sHT2.compare("T"))
		      {
			auto tupTmp = std::make_tuple(sHT2,ii,jj,matchtail.back().at(0).first,matchtail.back().at(0).second);
			tailvv.push_back(tupTmp);
		      }
		    else
		      {
			auto tupTmp = std::make_tuple(sHT2,ii,jj,matchtail.back().at(1).first,matchtail.back().at(1).second);
			tailvv.push_back(tupTmp);
		      }
		  }
		else
		  matchtail.pop_back();

		//		std::cout << "TrackStitcher: satisfied tail. " << std::endl;
	      } //tail
	    /*
	    std::cout << "abs(start1Dir.Dot(end2Dir)) " << std::abs(start1Dir.Dot(end2Dir)) << ", start1-end2.Mag(): " << (start1-end2).Mag() << std::endl;
	    std::cout << "abs(end1Dir.Dot(start2Dir)) " << std::abs(end1Dir.Dot(start2Dir)) << ", start2-end1.Mag(): " << (start2-end1).Mag() << std::endl;
	    std::cout << "abs(start1Dir.Dot(start2Dir)) " << std::abs(start1Dir.Dot(start2Dir)) << ", start1-start2.Mag(): " << (start1-start2).Mag() << std::endl;
	    std::cout << "abs(end1Dir.Dot(end2Dir)) " << std::abs(end1Dir.Dot(end2Dir)) << ", end1-end2.Mag(): " << (end1-end2).Mag() << std::endl;
	    std::cout << "sHT2 " << sHT2 << std::endl;
	    */
	  } // end c11||c12||c21||c22

	// We've been careful to pick the best jj match for this iith track head and tail.
	// Now we need to be sure that for the jjth track head/tail we don't have two ii trks.
	if (headvv.size())
	  {
	    int otrk = std::get<2>(headvv.back()); // jj'th track for this iith trk
	    // H or T of this jj'th trk we're matched to.
	    std::string sotrkht(std::get<0>(headvv.back()));
	    for (int kk=0;kk<ii;++kk)
	      {
		if (std::get<2>(fh.at(kk)) == otrk && !sotrkht.compare(std::get<0>(fh.at(kk)) ) )
		  {
		    // check matching sep and pick the best one. Either erase this
		    // headvv (and it'll get null settings later below) or null out
		    // the parameters in fh.
		    if (std::get<4>(headvv.back()) < std::get<4>(fh.at(kk)) && std::get<4>(headvv.back())!=0.0)
		      {
			auto tupTmp2 = std::make_tuple(std::string("NA"),kk,-12,0.0,0.0);
			fh.at(kk) = tupTmp2;
		      }
		    else if (std::get<4>(headvv.back())!=0.0)
		      {
			headvv.pop_back();
			break;
		      }
		  }
	      }
	  }
	if (tailvv.size())
	  {
	    int otrk = std::get<2>(tailvv.back());  // jj'th track for this iith trk
	    // H or T of this jj'th trk we're matched to.
	    std::string sotrkht(std::get<0>(tailvv.back()));
	    for (int kk=0;kk<ii;++kk)
	      {
		if (std::get<2>(ft.at(kk)) == otrk && !sotrkht.compare(std::get<0>(ft.at(kk)) ) )
		  {
		    // check matching sep and pick the best one. erase either this
		    // tailvv or null out the parameters in ft.
		    if (std::get<4>(tailvv.back()) < std::get<4>(ft.at(kk)) && std::get<4>(tailvv.back())!=0.0)
		      {
			auto tupTmp2 = std::make_tuple(std::string("NA"),kk,-12,0.0,0.0);
			ft.at(kk) = tupTmp2;
		      }
		    else if (std::get<4>(tailvv.back())!=0.0)
		      {
			tailvv.pop_back();
			break;
		      }
		  }
	      }
	  }


      } // jj


      auto tupTmp2 = std::make_tuple(std::string("NA"),ii,-12,0.0,0.0);
      // We always have our best 1-element tailvv and headvv for trk o at this point
      if (!headvv.size()) headvv.push_back(tupTmp2);
      if (!tailvv.size()) tailvv.push_back(tupTmp2);
      //      std::cout << "StitchAlg::FindHeadsAndTails: headvv, tailvv .get<0> is " << std::get<0>(headvv.back()) << ", " << std::get<0>(tailvv.back()) << std::endl;
      fh.push_back(headvv.back());
      ft.push_back(tailvv.back());

    } // ii

    //    std::cout << "fh.size, ft.size are " << fh.size() << ", " << ft.size() << std::endl;

}

void trkf::StitchAlg::FirstStitch(const std::vector<art::PtrVector <recob::Track>>::iterator itvvArg, const std::vector <recob::Track>::iterator itvArg)
{
    // take the vector of tracks, walk through each track's vectors of xyz, dxdydz, etc
    // and concatenate them into longer vectors. Use those to instantiate one new
    // Stitched-together track.
    std::vector<recob::tracking::Point_t> xyz;
    std::vector<recob::tracking::Vector_t> dxdydz;
    std::vector<recob::tracking::SMatrixSym55> cov;
    std::vector<recob::TrackTrajectory::PointFlags_t> flgs;
    //art::PtrVector<recob::Track>::const_iterator

    bool hasMomentum = true; //true only if all tracks have momentum
    for (auto it = (*itvvArg).begin(); it!=(*itvvArg).end(); ++it)
      {
	if ((*it).get()->HasMomentum()==false) {
	  hasMomentum = false;
	  break;
	}
      }

    size_t cnt(0);
    for (auto it = (*itvvArg).begin(); it!=(*itvvArg).end(); ++it)
      {

	cnt++;
	//	std::cout << "Stitching track cnt is: " << cnt << std::endl;
	for (size_t pt = 0; pt!=(*it).get()->NumberTrajectoryPoints(); pt++)
	  {
	    size_t ptHere(pt);
	    // ask if 1st character of 2-character string is an H. If so, reverse order
	    // of the concatenation. Note however that I've had my notion of head & tail backwards
	    // throughout this whole class! start1, end1 are really T, H, not H, T as I've labeled 'em till now.
	    // hence flip the direction if this character is a T/H, when expecting H/T! -- EC, 29-June-2014.

	    auto itvfHT = fHT.begin() + size_t (itvArg - fTrackVec.begin());
	    if ( (*itvfHT).size() &&
		 (
		  // was fHT.back()
		  (cnt==1 && !(*itvfHT).at(cnt-1).compare(0,1,"H")) ||
		   (cnt>1  && !(*itvfHT).at(cnt-2).compare(1,1,"T"))
		  )
		 )
	      ptHere = (*it).get()->NumberTrajectoryPoints() - pt - 1;

	    try
	      {
		xyz.push_back((*it).get()->LocationAtPoint(ptHere));
		//		std::cout << "Stitching track number " << cnt << " with TrajPt at ptHere " << ptHere << " at x,y,z: " << xyz.back().X() << ", " << xyz.back().Y() << ", " << xyz.back().Z() << std::endl;
		dxdydz.push_back( (hasMomentum ? (*it).get()->MomentumVectorAtPoint(ptHere) : (*it).get()->DirectionAtPoint(ptHere)) );
		flgs.push_back((*it).get()->FlagsAtPoint(ptHere));
		cov.push_back( (ptHere==0 ? (*it).get()->VertexCovariance() : (*it).get()->EndCovariance()));
	      }
	    catch (cet::exception &e)
	      {
		mf::LogVerbatim("TrackStitcher bailing. ") << " One or more of xyz, dxdydz, cov, mom, dQdx elements from original Track is out of range..." << e.what() << __LINE__;
		break;
	      }
	  }
      }

    /// TODO: sort according to spacepoint distance separation.
    /// As is, we're not sure we're not forming a stitched track with a (some)
    /// jump(s) and a reversal(s) of direction in it.

    //const recob::Track t(xyz,dxdydz,cov,dQdx,mom,ftNo++);
    const recob::Track t(recob::TrackTrajectory(std::move(xyz), std::move(dxdydz),  std::move(flgs), hasMomentum),
			 0, -1., 0, cov.front(), cov.back(), ftNo++);
    //const art::Ptr<recob::Track> t(xyz,dxdydz,cov,dQdx,mom,ftNo++);
    fTrackVec.insert(itvArg,t);

}

void trkf::StitchAlg::WalkStitch()
{

    art::PtrVector<recob::Track> compTrack;
    std::vector <std::string> trackStatus (fh.size(),"NotDone"); // H or T
    std::vector <std::string> HT2 ; // H or T


    for (unsigned int ii=0; ii<fh.size(); ++ii) // same as t.size()

      {
	if (!trackStatus.at(ii).compare("Done")) continue;

	const art::Ptr<recob::Track> th(ftListHandle, ii);
	compTrack.push_back(th);
	// Should there not be an HT.push_back here??

	// start with track 1: see if head goes anywhere, walk with it. Add to compTrack.
	// Go until the other tuple's other vtx string says "NA." Then change status string
	// of vtxsJoined to "Done" for that track.
	bool chain(true);
	int walk(ii);
	std::string sh(std::get<0>(fh.at(walk)));
	std::string st("NA");
	while (chain)
	  {
	    int hInd = -12;
	    int tInd = -12;
	    if (walk!=(int)ii)
	      {
		sh = std::get<0>(fh.at(walk));
		st = std::get<0>(ft.at(walk));
	      }

	    //	    std::cout << "StichAlg::WalkStitch(): Inside head chain. walk(track), sh, st, " << walk << ", " << sh << ", "<<st <<", connected to tracks: " << std::get<2>(fh.at(walk))<<", "  << std::get<2>(ft.at(walk)) <<std::endl;
	    if (!sh.compare("H") || !sh.compare("T"))
	      {

		hInd = std::get<2>(fh.at(walk)); // index of track that walk is connected to.
		//		std::cout << "WalkStitch(): hInd is " << hInd << std::endl;

		const art::Ptr<recob::Track> th2( ftListHandle, hInd );
		compTrack.push_back(th2);
		HT2.push_back("H"+sh);
	      }

	    if ((!st.compare("H") || !st.compare("T")))
	      {
		tInd = std::get<2>(ft.at(walk));
		//		std::cout << "WalkStitch(): tInd is " << tInd << std::endl;
		const art::Ptr<recob::Track> th2( ftListHandle, tInd );
		compTrack.push_back(th2);
		HT2.push_back("T"+st); // since we will eventually read from 0th element forward
	      }
	    if (hInd!=-12) walk = hInd;
	    if (tInd!=-12) walk = tInd;
	    if (!sh.compare("NA") && !st.compare("NA"))
		chain = false;

	    trackStatus.at(walk) = "Done";
	  } // while

	// It is possible that our first (ii'th) track had a head _and_ a tail match. Thus, we must
	// see if tail goes anywhere. walk with it. Insert, don't push_back, to compTrack.
	chain = true;
        walk = ii;
	sh = "NA";
	st = std::get<0>(ft.at(walk));
	while (chain)
	  {
	    int hInd = -12;
	    int tInd = -12;
	    if (walk!=(int)ii)
	      {
		sh = std::get<0>(fh.at(walk));
		st = std::get<0>(ft.at(walk));
	      }

	    //	    std::cout << "StichAlg::WalkStitch(): Inside tail chain. walk(track), sh, st, " << walk << ", " << sh << ", "<<st <<", connected to tracks: " << std::get<2>(fh.at(walk))<<", "  << std::get<2>(ft.at(walk)) <<std::endl;
	    if (!sh.compare("H") || !sh.compare("T"))
	      {

		hInd = std::get<2>(fh.at(walk)); // index of track that walk is connected to.
		//		std::cout << "WalkStitch(): hInd is " << hInd << std::endl;

		const art::Ptr<recob::Track> th2( ftListHandle, hInd );
		compTrack.insert(compTrack.begin(),th2);
		HT2.insert(HT2.begin(),sh+"H");
	      }

	    if ((!st.compare("H") || !st.compare("T")))
	      {
		tInd = std::get<2>(ft.at(walk));
		//		std::cout << "WalkStitch(): tInd is " << tInd << std::endl;
		const art::Ptr<recob::Track> th2( ftListHandle, tInd );
		compTrack.insert(compTrack.begin(),th2);
		HT2.insert(HT2.begin(),st+"T"); // since we will eventually read from 0th element forward
	      }
	    if (hInd!=-12) walk = hInd;
	    if (tInd!=-12) walk = tInd;
	    if (!sh.compare("NA") && !st.compare("NA"))
		chain = false;

	    trackStatus.at(walk) = "Done";
	  } // while


	// inside FirstStitch() push_back onto the vec<vec> of components and the vec of stitched composite.
	if (compTrack.size())
	  {
	    // protect against stitching a component twice, as when somehow the same track has been matched to a
	    // head and a tail of another track. remove the repeat appearances of that track.
	    for (auto iit = compTrack.begin(); iit!=compTrack.end();++iit)
	      {
		int cjit(0);
		for (auto jit = iit+1; jit!=compTrack.end();++jit)
		  {
		    //		    std::cout<< " cjit is ." << cjit << std::endl;
		    if (*iit==*jit)
		      {
			//			std::cout<< " About to do cjit erase." << std::endl;
			compTrack.erase(jit); //std::cout<< "Survived the compTrack jit erase." << std::endl;
			HT2.erase(HT2.begin()+cjit); //std::cout<< "Survived the HT2 cjit erase." << std::endl;
			break;
		      }
		    cjit++;
		  }
	      }
	    fTrackComposite.push_back(compTrack);
	    fHT.push_back(HT2);
	    //	    std::cout << "WalkStitch:: calling FirstStitch(). fTrackComposite.size(), fHT.size() " << fTrackComposite.size()<< ", " << fHT.size()  << std::endl;
	    //std::cout << "WalkStitch:: calling FirstStitch(). fTrackComposite.back().size(), fHT.back().size() " << fTrackComposite.back().size()<< ", " << fHT.back().size()  << std::endl;
	    /*
	    for (unsigned int ll=0; ll<fHT.back().size() ; ++ll)
	      {
		std::cout << fHT.back().at(ll); std::cout << ", ";
	      }
		   std::cout << "." << std::endl;
	    */
	    // want the last vector of fTrackComposite, and will want to insert on fTrackVec, hence
	    // -1 and not -1, respectively.
	    FirstStitch(fTrackComposite.end()-1, fTrackVec.end());
	  }
	compTrack.clear();
	HT2.clear();


      } // end ii loop on head/tail vector of tuples.

}

    // This method will join separate composite tracks in which one common component is joined by its head in one
    // and by its tail in another. This can happen if the common component has a higher track index than either of
    // the the two it is separately stitched to. e.g., ____(1) -------(4) ___________(3).
    //
bool  trkf::StitchAlg::CommonComponentStitch()
{
  // "os" for outer scope.
  int osciit(-12), oscjit(-12);
  std::vector < art::PtrVector<recob::Track> >::iterator osiComposite, osjComposite;
  art::PtrVector < recob::Track >::iterator osiAgg, osjAgg;


  bool match(false);
  int ciit(0);
  for (auto iit = fTrackComposite.begin(); iit!=fTrackComposite.end()&&!match; ++iit)
    {
      ciit++;
      int cjit(ciit);
      for (auto jit = iit+1; jit!=fTrackComposite.end()&&!match; ++jit)
	{
	  cjit++;
	  for (auto iiit = iit->begin(); iiit!=iit->end()&&!match; ++iiit)
	    {
	      for (auto jjit = jit->begin(); jjit!=jit->end()&&!match; ++jjit)
		{
		  if (*iiit == *jjit) // 2 components from 2 different composites are the same
		    {
		      // head is attached to one trk and tail to another.
		      //		      std::cout << "StitchAlg::CommonComponentStitch: We have two aggregate tracks that have a common component and need to be further stitched. " << std::endl;

		      match = true;
		      osiComposite = iit;
		      osjComposite = jit;
		      osciit = ciit;
		      oscjit = cjit;
		      osiAgg = iiit;
		      osjAgg = jjit; // yes, unneeded, but we keep it for notational clarity

		    }
		}

	    }
	}
    }
  if (!match) return match;

  // Proceed to stitch 'em all together, dropping the one redundant component. Then
  // erase the first occurence of the composite and the matching aggregate trk.


  //  std::cout << "StitchAlg::CommonComponentStitch: pre erase: " << osiComposite->size() << std::endl;
  (*osiComposite).erase(osiAgg);                // erase redundant component track
                                                // do not erase this fHT element, however

  //  std::cout << "StitchAlg::CommonComponentStitch: post erase: " << osiComposite->size() << std::endl;
  // std::cout << "StitchAlg::CommonComponentStitch: fHT.size(): " << fHT.size() << std::endl;

  // Next is a loop over all remaining components in osiComposite.
  // insert the non-redundant osciit tracks onto front (back) of osjit
  // insert the non-redundant osciit vtx links onto front (back) of fHT.begin()+oscjit-1

  std::vector< std::vector<std::string> >::iterator siit(fHT.begin()+osciit-1);
  std::vector< std::vector<std::string> >::iterator sjit(fHT.begin()+oscjit-1);
  size_t itdiff(osiComposite->end()-osiComposite->begin());
  if (osjAgg == osjComposite->begin())
    {

      //      std::cout << "StitchAlg::begin insert starting: " << std::endl;
      // std::cout << "StitchAlg::CommonComponentStitch: itdiff: " << itdiff << std::endl;
      //std::cout << "StitchAlg::CommonComponentStitch: osiComposite.end-begin: " <<  osiComposite->end()-osiComposite->begin()<< std::endl;
      //std::cout << "StitchAlg::CommonComponentStitch: osjComposite.end-begin: " <<  osjComposite->end()-osjComposite->begin()<< std::endl;
      (*osjComposite).insert(osjComposite->begin(),osiComposite->begin(),osiComposite->begin()+itdiff);
      //std::cout << "StitchAlg::CommonComponentStitch: siit.end-begin: " <<  siit->end()-siit->begin()<< std::endl;
      //std::cout << "StitchAlg::CommonComponentStitch: sjit.end-begin: " <<  sjit->end()-sjit->begin()<< std::endl;
      (*sjit).insert(sjit->begin(),siit->begin(),siit->begin()+itdiff);
      //std::cout << "StitchAlg::begin insert done: " << std::endl;
    }
  else if (osjAgg == (osjComposite->end()-1))
    {
      //      std::cout << "StitchAlg::end insert starting: " << std::endl;
      (*osjComposite).insert(osjComposite->end(),osiComposite->begin(),osiComposite->begin()+itdiff);
      (*sjit).insert(sjit->end(),siit->begin(),siit->begin()+itdiff);
      //      std::cout << "StitchAlg::end insert done: " << std::endl;
    }

  //  std::cout << "StitchAlg:: 1: " << std::endl;
  fTrackVec.erase(fTrackVec.begin()+oscjit-1);   // erase old Stitched Track, which we'll recreate now...

  //  std::cout << "StitchAlg:: 2: " << std::endl;
  FirstStitch(osjComposite,fTrackVec.begin()+oscjit-1); // Create new Stitched Track
  //  fTrackComposite.insert(fTrackComposite.begin()+oscjit-1-1,*osjComposite);           // erase old composite Track
  //  std::cout << "StitchAlg:: 3: " << std::endl;
  fTrackVec.erase(fTrackVec.begin()+osciit-1);   // erase old Stitched Track
  //  std::cout << "StitchAlg:: 6: " << std::endl;
  fTrackComposite.erase(osiComposite);           // erase old composite Track
  //  std::cout << "StitchAlg:: 4: " << std::endl;
  fHT.erase(fHT.begin()+osciit-1);               // erase old vec of vtx links
  //  std::cout << "StitchAlg:: 5: " << std::endl;

  return match;

}   // end of bool CommonComponentStitch()
