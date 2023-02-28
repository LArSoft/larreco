////////////////////////////////////////////////////////////////////////
//
// \file SmallClusterFinderAlg.cxx
//
// \author corey.adams@yale.edu
//
// This algorithm is designed to find small clusters that could correspond to gammas
// or low energy electrons.
//
/*	There are two parameters that matter from the fcl file:
                fNHitsInClust is the number of hits that should be in these small clusters
                        ^-- Gamma events seem to rarely have more than 4 hits in the cluster
                        ^-- SN events are unclear.  Should this even be used for SN?
                fRadiusSizePar is the distance (in cm) between the small clusters and any other hits.

        This algorithm sorts the hits by plane, and then looks at each hit individually.  If
        there is a hit within RadiusSizePar, it gets added to a local list.  All other hits
        are ignored.  Then, if the number of hits that got added to the local list is greater
        then NHitsInClust, the original hit is ignored.  If it's less, the original hit is
        presumed to be part of a very small (or single hit) cluster.  So its added to the list
        of hits in the small cluster.

        All of the small clusters are then split apart into groups in the way you would expect.
        Each cluster is assigned an ID number to distinguish it, and the hits that aren't
        identified as small clusters all end up in the "leftover" cluster.  The numbering scheme
        is ID = 100*iPlane + Cluster on that plane, and the leftover hits are the first (0th)
        cluster written out.

        -Corey
*/
//
//
///////////////////////////////////////////////////////////////////////

#include <iostream>

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/SmallClusterFinderAlg.h"

// ***************** //

//constructor with parameters:
cluster::SmallClusterFinderAlg::SmallClusterFinderAlg(fhicl::ParameterSet const& pset)
{
  fNPlanes = geom->Nplanes();
  fRadiusSizePar = pset.get<double>("RadiusSizePar");
  fNHitsInClust = pset.get<double>("NHitsInClust");
  verbose = pset.get<bool>("Verbose");
  ClearandResizeVectors();
}

// ***************** //
// This method actually makes the clusters.
void cluster::SmallClusterFinderAlg::FindSmallClusters(
  util::GeometryUtilities const& gser,
  detinfo::DetectorClocksData const& clockData,
  detinfo::DetectorPropertiesData const& detProp,
  std::vector<art::Ptr<recob::Hit>> allHits)
{
  ///These lines determine the conversion factors to take wires and times to CMs
  fDriftVelocity = detProp.DriftVelocity(detProp.Efield(), detProp.Temperature());
  fWirePitch = geom->WirePitch();
  fTimeTick = sampling_rate(clockData) / 1000.;
  fWiretoCm = fWirePitch;
  fTimetoCm = fTimeTick * fDriftVelocity;
  fWireTimetoCmCm = (fTimeTick * fDriftVelocity) / fWirePitch;

  ClearandResizeVectors();

  // Catch the case were there are no hits in the event:
  if (allHits.size() == 0) {
    if (verbose) std::cout << " no hits received! exiting " << std::endl;
    return;
  }

  art::Ptr<recob::Hit> theHit;

  //sort the hits into hits by plane
  for (std::vector<art::Ptr<recob::Hit>>::iterator HitIter = allHits.begin();
       HitIter != allHits.end();
       HitIter++) {
    theHit = *HitIter;
    unsigned int p(0), w(0), t(0), cs(0);  //c=channel, p=plane, w=wire, but what is t?
    GetPlaneAndTPC(*HitIter, p, cs, t, w); //Find out what plane this hit is on.
    //add this hit to the list specific to this plane
    hitlistbyplane[p].push_back(theHit);
  } // End loop on hits.

  // Want to check that the wires are OK for each event, and by each wire.
  // This could be done in HitFinder, but I'm doing it here because I want to do
  // some things specifically for Small Clusters... Flags to look for:  Does
  // this
  // wire have more hits that the other wires that also have hits? Significantly
  // more?       Does the RMS of this hit lie well below the peak of the hit?
  // Significantly below?        More?   If a bad wire is found, remove all of
  // it's from the hitlists so they aren't clustered.

  //make the refined hit list for each plane.
  for (unsigned int ip = 0; ip < fNPlanes; ip++) {
    hitlistrefined[ip] = CreateHighHitlist(gser, hitlistbyplane[ip], hitlistleftover[ip]);

    //Check that the lists are populated correctly:
    if (verbose)
      std::cout << "At plane " << ip << ", found " << hitlistrefined[ip].size() << " hits with "
                << hitlistleftover[ip].size() << " leftover" << std::endl;
    //add the leftover hits to the correct object:
  }

  //Now we have all of the gammas.  This is good!
  //Want to split all of the gammas up into individual clusters.

  //I was going to use lists instead of vectors to do this, because lists are more efficient
  //at insertion and removal at arbitrary places.  But, the number of gammas is so small that
  //it just doesn't seem worth it
  //
  //Also, larsoft is so slow on its own that i think the difference is undetectable

  // The method to split the gammas into individual clusters that I want is:
  // Use the existing method below to find all the hits within a certain distance from a given hit
  // Take ALL of those hits and write them to a cluster.
  // Then, remove those hits from the list of gammas.

  // Repeat, until the list of gammas is empty.

  // There is an issue of knowing which hits to remove from the total gamma hit list.
  // To solve this, I'm going to overload SelectLocalHitList to take a reference to a vector
  // of ints as an argument.  When a hit is added to the local hit list,
  // it's index will be added to a vector that is returned by reference.
  // It's important to remove these hits in reverse order.  This is because removing a hit
  // changes the index of all of the hits after it

  //Now we need to take the lists of small clusters and sort it into the individual spots
  //going to end up with smallClustList[plane][iClust][Hit]
  //loop over planes of hits:

  for (unsigned int iplane = 0; iplane < fNPlanes; iplane++) {

    if (hitlistrefined[iplane].size() == 0) continue;

    //write the rest of the gammas one by one, in clusters:
    int i = 1;

    std::vector<art::Ptr<recob::Hit>> splittingVector = hitlistrefined[iplane];

    while (splittingVector.size() != 0) {

      // std::cout << "\nThe hits remaining to be spilt are:" << std::endl;
      //for (unsigned int j = 0; j < splittingVector.size();j++){
      //	std::cout << *splittingVector[j] << std::endl;
      //}

      //find the first small cluster of gammas:
      std::vector<int> index;
      std::vector<art::Ptr<recob::Hit>> thiscluster;
      thiscluster.clear();
      index.clear();

      //Just use the first hit in the list of gammas:

      art::Ptr<recob::Hit> theHit = splittingVector.front(); //grab a hit from the list

      double time = theHit->PeakTime();
      unsigned int plane(0), cstat(0), tpc(0), wire(0);
      GetPlaneAndTPC(theHit, plane, cstat, tpc, wire);

      SelectLocalHitlist(gser, splittingVector, thiscluster, wire, time, fRadiusSizePar, index);

      if (verbose)
        std::cout << "Done writing " << thiscluster.size() << " hits to cluster with ID "
                  << plane * 100 + i << std::endl;
      //make sure to add these hits to the object that stores them:
      smallClustList[plane].push_back(thiscluster);

      //Lastly, remove the gammas just clustered from the refinded hit list
      while (index.size() != 0) {
        splittingVector.erase(splittingVector.begin() + (index.back()));
        index.pop_back();
      }
      i++;
    }
  }
}

// ************************************* //
// Clear and resize - exactly what it sounds like
void cluster::SmallClusterFinderAlg::ClearandResizeVectors()
{
  smallClustList.clear();
  hitlistbyplane.clear();
  hitlistrefined.clear();
  hitlistleftover.clear();
  smallClustList.resize(fNPlanes);
  hitlistbyplane.resize(fNPlanes);
  hitlistrefined.resize(fNPlanes);
  hitlistleftover.resize(fNPlanes);
}

/*		This is the method takes a list of hits ("hitlist") and compares each hit to a 2D
        location that is specified by a wire ("wire_start") and time ("time_start").  If the
        hit being examined is farther away than a specified distance ("radlimit", in cm) then
        the hit is excluded.  If the hit is within that distance, it's added.
*/
void cluster::SmallClusterFinderAlg::SelectLocalHitlist(
  util::GeometryUtilities const& gser,
  std::vector<art::Ptr<recob::Hit>> hitlist,
  std::vector<art::Ptr<recob::Hit>>& hitlistlocal,
  double wire_start,
  double time_start,
  double radlimit) const
{
  //loop over the hits in "hitlist", which should contain the hits we're selecting from
  for (std::vector<art::Ptr<recob::Hit>>::const_iterator hitIter = hitlist.begin();
       hitIter != hitlist.end();
       hitIter++) {
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime();
    unsigned int plane, cstat, tpc, wire;
    GetPlaneAndTPC(theHit, plane, cstat, tpc, wire);
    //we now know which wire and what time this hit occurred on

    //calculate linear distance from start point and orthogonal distance from axis
    double linear_dist = gser.Get2DDistance(wire, time, wire_start, time_start);

    if (linear_dist < radlimit) hitlistlocal.push_back(theHit);
  }
  return;
}

//This method is identical to the other  method by the same name except that
//it keeps track of the location of the hits selected.  That is, index is filled with
//the indices of the selected hits in the hitlist vector (the input vector)
void cluster::SmallClusterFinderAlg::SelectLocalHitlist(
  util::GeometryUtilities const& gser,
  std::vector<art::Ptr<recob::Hit>> hitlist,
  std::vector<art::Ptr<recob::Hit>>& hitlistlocal,
  double wire_start,
  double time_start,
  double radlimit,
  std::vector<int>& index) const
{
  //loop over the hits in "hitlist", which should contain the hits we're selecting from
  int i = 0; //i keeps track of the index of the hit.
  for (std::vector<art::Ptr<recob::Hit>>::const_iterator hitIter = hitlist.begin();
       hitIter != hitlist.end();
       hitIter++) {
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime();
    unsigned int plane, cstat, tpc, wire;
    GetPlaneAndTPC(theHit, plane, cstat, tpc, wire);
    //we now know which wire and what time this hit occurred on

    //calculate linear distance from start point and orthogonal distance from axis
    double linear_dist = gser.Get2DDistance(wire, time, wire_start, time_start);

    if (linear_dist < radlimit) {
      hitlistlocal.push_back(theHit);
      index.push_back(i);
    }
    i++;
  }
  //std::cout << "in select local hit list, index is:" << std::endl;
  //for (int j = 0; j < index.size();j++) std::cout << index[j] << " ";

  //need to make sure the index array is in order.  Sort it!
  std::sort(index.begin(), index.end());

  return;
}

/*		This method sorts this hitlist (which should only be on a single plane)

*/
std::vector<art::Ptr<recob::Hit>> cluster::SmallClusterFinderAlg::CreateHighHitlist(
  util::GeometryUtilities const& gser,
  std::vector<art::Ptr<recob::Hit>> const& hitlist,
  std::vector<art::Ptr<recob::Hit>>& leftovers) const
{

  std::vector<art::Ptr<recob::Hit>>
    hitlist_total; //This is the final result, a list of hits that are small clusters

  std::vector<art::Ptr<recob::Hit>> hitlistlocal;

  for (unsigned int ix = 0; ix < hitlist.size(); ix++) {

    art::Ptr<recob::Hit> const& theHit = hitlist[ix]; //grab a hit from the list

    double time = theHit->PeakTime();
    unsigned int plane(0), cstat(0), tpc(0), wire(0);
    //std::cout << "The hit is " << (*theHit) << std::endl;
    GetPlaneAndTPC(theHit, plane, cstat, tpc, wire);
    //use the wire and time of this hit as a seed.
    // ^^^^^ This could probably be optimized?

    //get ALL of the hits from hitlist that are within the distance fRadiusSizePar of the seed hit.
    SelectLocalHitlist(gser, hitlist, hitlistlocal, (double)wire, time, fRadiusSizePar);

    if (hitlistlocal.size() < fNHitsInClust) {
      hitlist_total.push_back(theHit); //Add this hit if there are less than fNHitsInClust nearby.
      if (verbose)
        std::cout << " adding hit @ w,t " << wire << " " << time << " on plane " << plane
                  << std::endl;
    }
    else {
      //Add this hit to the leftover pile
      leftovers.push_back(theHit);
    }

    hitlistlocal.clear(); //clear the local hit list, and look at the next hit.
  }

  /*
        This method could definitely be optimized.  It creates a local hit list for each particle,
        while if there is a local hit list that is sufficiently separated from all others it should be OK to
        add them all at once.  This is a task for future coders!
*/

  return hitlist_total;
}

// ******************************* //
int cluster::SmallClusterFinderAlg::GetPlaneAndTPC(art::Ptr<recob::Hit> a, //the hit
                                                   unsigned int& p,        //plane
                                                   unsigned int& /*cs*/,   //cryostat
                                                   unsigned int& t,        //time
                                                   unsigned int& w) const  //wire
{
  art::ServiceHandle<geo::Geometry const> geom;
  unsigned int channel = a->Channel();
  geom->ChannelToWire(channel);
  p = a->WireID().Plane;
  t = a->PeakTime();
  w = a->WireID().Wire;

  return 0;
}

//Function to return the small clusters by plane
std::vector<std::vector<art::Ptr<recob::Hit>>>
cluster::SmallClusterFinderAlg::GetSmallClustersByPlane(unsigned int iPlane)
{
  if (iPlane < fNPlanes)
    return smallClustList[iPlane];
  else {
    std::vector<std::vector<art::Ptr<recob::Hit>>> vec;
    return vec;
  }
}

std::vector<art::Ptr<recob::Hit>> cluster::SmallClusterFinderAlg::GetLeftoversByPlane(
  unsigned int iPlane)
{
  if (iPlane < fNPlanes)
    return hitlistleftover[iPlane];
  else {
    std::vector<art::Ptr<recob::Hit>> vec;
    return vec;
  }
}
