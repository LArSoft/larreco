//////////////////////////////////////////////////////////////////////
///
/// ClusterCrawlerAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for crawling along a string of hits to make line clusters
/// Technical note in MicroBooNE docdb #2831
///
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm> // std::fill(), std::find(), std::sort()...

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/Exception.h"

// LArSoft libraries
#include "SimpleTypesAndConstants/RawTypes.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoAlg/ClusterCrawlerAlg.h"
#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"

namespace cluster {


////////////////////////////////////////////////
void ClusterCrawlerAlg::ClusterLoop2()
{
// Version 2 of ClusterLoop.
unsigned int ii, iwire, jwire, iht, jht;

unsigned int nwires = fLastWire - fFirstWire - 1;
bool clusterAdded, success;
std::vector<std::pair<unsigned short, unsigned short>> kinkIndex;

for(ii = 0; ii < nwires; ++ii) {
// decide which way to step given the sign of fStepCrawlDir
if(fStepCrawlDir > 0) {
// step DS
iwire = fFirstWire + ii;
jwire = iwire + 1;
} else {
// step US
iwire = fLastWire - ii;
jwire = iwire - 1;
}
if(iwire < fFirstWire || iwire > fLastWire - 1) {
std::cout<<"Bad iwire indexing "<<fFirstWire<<" "<<iwire<<" "<<fLastWire<<"\n";
exit(1);
}
if(jwire < fFirstWire || jwire > fLastWire - 1) {
std::cout<<"Bad jwire indexing "<<fFirstWire<<" "<<jwire<<" "<<fLastWire<<"\n";
exit(1);
}
// skip bad wires or no hits on the wire
if(WireHitRange[iwire].first < 0) continue;
if(WireHitRange[jwire].first < 0) continue;
unsigned int ifirsthit = (unsigned int)WireHitRange[iwire].first;
unsigned int ilasthit = (unsigned int)WireHitRange[iwire].second;
unsigned int jfirsthit = (unsigned int)WireHitRange[jwire].first;
unsigned int jlasthit = (unsigned int)WireHitRange[jwire].second;
for(iht = ifirsthit; iht < ilasthit; ++iht) {
prt = (fDebugPlane == (int)plane && (int)iwire == fDebugWire && std::abs((int)fHits[iht].PeakTime() - fDebugHit) < 20);
if(prt) mf::LogVerbatim("CC")<<"Found debug hit "<<plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" inClus "<<inClus[iht];
if(inClus[iht] != 0) continue;
clusterAdded = false;
for(jht = jfirsthit; jht < jlasthit; ++jht) {
if(inClus[iht] != 0) continue;
if(inClus[jht] != 0) continue;
// start a cluster with these two hits
StartTraj();
traj.clear();
fcl2hits.clear();
fcl2hits.push_back(iht);
fcl2hits.push_back(jht);
if(prt) mf::LogVerbatim("CC")<<" checking ClusterHitsOK with jht "<<plane<<":"<<fHits[jht].WireID().Wire<<":"<<(int)fHits[jht].PeakTime();
// Ensure that the hits StartTick and EndTick have the proper overlap
if(!ClusterHitsOK(-1)) continue;
StartTraj();
StepCrawl(success);
if(!success) {
mf::LogVerbatim("CC")<<"ClusterLoop2: Bad return from StepCrawl. Seed hit "<<plane<<":"<<fHits[fcl2hits[0]].WireID().Wire<<":"<<(int)fHits[fcl2hits[0]].PeakTime()<<" fcl2hits size "<<fcl2hits.size();
for(unsigned short kk = 0; kk < fcl2hits.size(); ++kk) {
//              mf::LogVerbatim("CC")<<"kk "<<kk<<" "<<fcl2hits[kk]<<" inClus "<<inClus[fcl2hits[kk]];
if(inClus[fcl2hits[kk]] != 0) mf::LogVerbatim("CC")<<"Bad inClus "<<kk<<" "<<fcl2hits[kk]<<" inClus "<<inClus[fcl2hits[kk]];
} // kk
continue;
}
//          mf::LogVerbatim("CC")<<"StepCrawl done: fcl2hits "<<fcl2hits.size()<<" traj size "<<traj.size()<<" seed hit "<<plane<<":"<<fHits[fcl2hits[0]].WireID().Wire<<":"<<(int)fHits[fcl2hits[0]].PeakTime();
if(prt) mf::LogVerbatim("CC")<<"StepCrawl done: fcl2hits "<<fcl2hits.size()<<" traj size "<<traj.size();
if(fcl2hits.size() < 3) continue;
if(traj.size() < 2) continue;
StepCrawlFindKinks(kinkIndex);
for(unsigned short iknk = 0; iknk < kinkIndex.size(); ++iknk)
StoreTraj(kinkIndex[iknk].first, kinkIndex[iknk].second, 990);
clusterAdded = true;
break;
} // jht
if(clusterAdded) break;
} // iht
} // iwire

prt = false;
traj.clear();
fcl2hits.clear();

} // ClusterLoop2

////////////////////////////////////////////////
void ClusterCrawlerAlg::StoreTraj(unsigned short fromIndex, unsigned short toIndex, unsigned short procCode)
{
// Creates a cluster using hits associated with trajectory points fromIndex to toIndex

bool properOrder;
unsigned short itj, end;
unsigned short lastPt = traj.size() - 1;
float slp;
std::array<float, 2> dir;

std::vector<TrajPoint> trajTmp;
std::vector<unsigned int> fcl2hitsTmp;
bool doRestore = false;

if(fromIndex > lastPt || toIndex > lastPt) {
mf::LogError("CC")<<"StoreTraj: Crazy indices "<<fromIndex<<" "<<toIndex<<" for traj size "<<traj.size();
return;
} // error

if(fromIndex > 0 && toIndex < lastPt) {
mf::LogError("CC")<<"StoreTraj: Invalid indices "<<fromIndex<<" "<<toIndex<<" for traj size "<<traj.size();
return;
} // error

// make a copy of traj and fcl2hits if we need to clobber parts of them
if(fromIndex > 0 || toIndex < lastPt) {
trajTmp = traj;
fcl2hitsTmp = fcl2hits;
doRestore = true;
if(fromIndex > 0) {
// erase the first part of traj
traj.erase(traj.begin(), traj.begin() + fromIndex - 1);
// erase the first part of fcl2hits
fcl2hits.erase(fcl2hits.begin(), fcl2hits.begin()+traj[0].NumHits);
} else {
// erase the last part of traj and fcl2hits
traj.resize(toIndex+1);
fcl2hits.resize(traj[toIndex].NumHits);
}
} // make a copy of traj and fcl2hits

itj = traj.size() - 1;
properOrder = (traj[itj].HitPos[0] < traj[0].HitPos[0]);
if(!properOrder) ReverseTraj();
// Define the clBegin (end = 0) and clEnd (end = 1) cluster parameters
for(end = 0; end < 2; ++end) {
if(end == 0) {
dir = traj[0].Dir;
} else {
dir = traj[itj].Dir;
}
if(std::abs(dir[0]) < 1e-3) {
// Slope way too large - fake it
if(dir[0] > 0) { slp = 999; } else { slp = 999; }
} else {
slp = dir[1] / dir[0];
} // std::abs(dir[0]) < 1e-3
// Define clBegin and clEnd parameters
// TODO: More careful treatment of slope errors
if(end == 0) {
clBeginSlp = slp; clBeginSlpErr = 0.1 * std::abs(slp); clBeginChgNear = 0;
} else {
clEndSlp = slp; clEndSlpErr = 0.1 * std::abs(slp); clEndChgNear = 0;
}
} // end
if(prt) mf::LogVerbatim("CC")<<"Store cluster with size "<<fcl2hits.size();
clProcCode = procCode;
// special stop code to indicate special crawling
clStopCode = 8;
TmpStore();

if(doRestore) {
traj = trajTmp;
fcl2hits = fcl2hitsTmp;
}
} // StoreTraj

////////////////////////////////////////////////
void ClusterCrawlerAlg::StartTraj()
{
// Start a simple (seed) trajectory using two hits in fcl2hits.
// The direction of the trajectory is from fcl2hits[0] to fcl2hits[1].
// fcl2hits[1] is then deleted, leaving one trajectory point.
// It will be re-found by StepCrawl along with any nearby hits that will
// constitute another trajectory point

traj.clear();
if(fcl2hits.size() != 2) return;

unsigned int iht, jht;
iht = fcl2hits[0];
jht = fcl2hits[1];

std::array<float, 2> pos, dir;

TrajPoint tp;
// position
pos[0] = (float)fHits[iht].WireID().Wire;
pos[1] = fHits[iht].PeakTime() * fScaleF;
tp.HitPos = pos; tp.Pos = pos;
// direction
dir[0] = fHits[jht].WireID().Wire - fHits[iht].WireID().Wire;
dir[1] = (fHits[jht].PeakTime() - fHits[iht].PeakTime()) * fScaleF;
float ur = sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
dir[0] /= ur; dir[1] /= ur;
tp.Dir = dir;
tp.Ang = atan2(dir[1], dir[0]);
tp.AngErr = 0.5;
tp.Chg = fHits[iht].Integral();
tp.NumHits = 1;
traj.push_back(tp);

// initialize other variables used by StepCrawl
fAveHitRMS = fHits[iht].RMS();
fAveChg = (fHits[iht].Integral() + fHits[jht].Integral()) / 2;
fAveChgRMS = 0.3 * fAveChg;

fcl2hits.pop_back();

} // StartTraj

} // namespace cluster