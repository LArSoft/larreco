/////////////////////////////////////////////////////////////////////
///
/// HoughBaseAlg class
///
/// Ben Carls, bcarls@fnal.gov
///
/// The Hough transform employed by fuzzy clustering is a heavily modified variant of the original 
/// Hough line code. It identifies lines in fuzzy clusters (e.g. muon tracks) and splits them off
/// into new clusters. 
///
/// The algorithm is based on the Progressive Probabilistic Hough Transform (PPHT).
/// See the following paper for details:
///
/// J. Matas et al., Robust Detection of Lines Using the Progressive Probabilistic Hough Transform,
/// Computer Vision and Image Understanding, Volume 78, Issue 1, April 2000, Pages 119–137
///
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include <boost/bind.hpp>

#include <TF1.h>
#include <TH1D.h>
#include <TStopwatch.h>

#include "CLHEP/Random/RandFlat.h"

#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
 
#include "Filters/ChannelFilter.h"
#include "RecoAlg/HoughBaseAlg.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#define PI TMath::Pi()
#define a0  0 /*-4.172325e-7f*/   /*(-(float)0x7)/((float)0x1000000); */
#define a1 1.000025f        /*((float)0x1922253)/((float)0x1000000)*2/Pi; */
#define a2 -2.652905e-4f    /*(-(float)0x2ae6)/((float)0x1000000)*4/(Pi*Pi); */
#define a3 -0.165624f       /*(-(float)0xa45511)/((float)0x1000000)*8/(Pi*Pi*Pi); */
#define a4 -1.964532e-3f    /*(-(float)0x30fd3)/((float)0x1000000)*16/(Pi*Pi*Pi*Pi); */
#define a5 1.02575e-2f      /*((float)0x191cac)/((float)0x1000000)*32/(Pi*Pi*Pi*Pi*Pi); */
#define a6 -9.580378e-4f    /*(-(float)0x3af27)/((float)0x1000000)*64/(Pi*Pi*Pi*Pi*Pi*Pi); */

#define _sin(x) ((((((a6*(x) + a5)*(x) + a4)*(x) + a3)*(x) + a2)*(x) + a1)*(x) + a0)
#define _cos(x) _sin(TMath::Pi()*0.5 - (x))

//------------------------------------------------------------------------------
cluster::HoughBaseAlg::HoughBaseAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

//------------------------------------------------------------------------------
cluster::HoughBaseAlg::~HoughBaseAlg()
{
}

//------------------------------------------------------------------------------
void cluster::HoughBaseAlg::reconfigure(fhicl::ParameterSet const& pset)
{ 
  fMaxLines                       = pset.get< int    >("MaxLines"                       );
  fMinHits                        = pset.get< int    >("MinHits"                        );
  fSaveAccumulator                = pset.get< int    >("SaveAccumulator"                );
  fNumAngleCells                  = pset.get< int    >("NumAngleCells"                  );
  fMaxDistance                    = pset.get< float  >("MaxDistance"                    );
  fMaxSlope                       = pset.get< float  >("MaxSlope"                       );
  fRhoZeroOutRange                = pset.get< int    >("RhoZeroOutRange"                );
  fThetaZeroOutRange              = pset.get< int    >("ThetaZeroOutRange"              );
  fRhoResolutionFactor            = pset.get< float  >("RhoResolutionFactor"            );
  fPerCluster                     = pset.get< int    >("HitsPerCluster"                 );
  fMissedHits                     = pset.get< int    >("MissedHits"                     );
  fMissedHitsDistance             = pset.get< float  >("MissedHitsDistance"             );
  fMissedHitsToLineSize           = pset.get< float  >("MissedHitsToLineSize"           );
  return;
}

//------------------------------------------------------------------------------
cluster::HoughTransform::HoughTransform()
{  
  //m_accum=NULL;
}


//------------------------------------------------------------------------------
size_t cluster::HoughBaseAlg::Transform(std::vector<art::Ptr<recob::Hit> > const& hits,
					std::vector<unsigned int>                 *fpointId_to_clusterId,
					unsigned int                               clusterId, // The id of the cluster we are examining
				        unsigned int                               *nClusters,
                                        std::vector<protoTrack>                    *linesFound
					)
{

  int nClustersTemp = *nClusters;
  
  HoughTransform c;

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  filter::ChannelFilter chanFilt;

  //  uint32_t     channel = hits[0]->Channel();
  unsigned int wire    = 0;
  unsigned int wireMax = 0;
  unsigned int cs=hits[0]->WireID().Cryostat;
  unsigned int t=hits[0]->WireID().TPC;
  unsigned int p=hits[0]->WireID().Plane;   


  //mf::LogInfo("HoughBaseAlg") << "nClusters is: " << *nClusters;


  //size_t cinctr = 0;
  //geo::View_t    view = geom->Cryostat(cs).TPC(t).Plane(p).View();
  //  geo::SigType_t sigt = geom->SignalType(channel);
  geo::SigType_t sigt = geom->Cryostat(cs).TPC(t).Plane(p).SignalType();   
  std::vector<int> skip;  
  
  double wire_pitch[3];
  wire_pitch[0]= geom->WirePitch(0,1,0);
  wire_pitch[1]= geom->WirePitch(0,1,1);
  wire_pitch[2]= geom->WirePitch(0,1,2);

  //factor to make x and y scale the same units
  double xyScale[3];
  xyScale[2]= .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
  xyScale[0]=xyScale[2] * detprop->SamplingRate()/wire_pitch[0];
  xyScale[1]=xyScale[2] * detprop->SamplingRate()/wire_pitch[1];
  xyScale[2]*= detprop->SamplingRate()/wire_pitch[2];

  float tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature())*1.e-3 * detprop->SamplingRate();
  //tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns


  //mf::LogInfo("HoughBaseAlg") << "xyScale: " << xyScale;
  //mf::LogInfo("HoughBaseAlg") << "tickToDist: " << tickToDist;
  
  int x, y;
  //unsigned int channel, plane, wire, tpc, cstat;
  //there must be a better way to find which plane a cluster comes from
  int dx = geom->Cryostat(hits[0]->WireID().Cryostat).TPC(hits[0]->WireID().TPC).Plane(hits[0]->WireID().Plane).Nwires();//number of wires 
  //int dy = hits[0]->Wire()->NSignal();//number of time samples. 
  int dy = detprop->NumberTimeSamples();//number of time samples. 
  skip.clear();
  skip.resize(hits.size());
  std::vector<int> listofxmax;
  std::vector<int> listofymax;  
  std::vector<int> hitsTemp;        //indecies ofcandidate hits
  std::vector<int> sequenceHolder; //wire numbers of hits in list
  std::vector<int> channelHolder; //channels of hits in list
  std::vector<int> currentHits;    //working vector of hits 
  std::vector<int> lastHits;       //best list of hits
  std::vector<art::Ptr<recob::Hit> > clusterHits;
  float indcolscaling = 0.;       //a parameter to account for the different 
        			   ////characteristic hit width of induction and collection plane
  /// \todo: the collection plane's characteristic hit width's are, 
  /// \todo: on average, about 5 time samples wider than the induction plane's. 
  /// \todo: this is hard-coded for now.
  if(sigt == geo::kInduction)
    indcolscaling = 0.;
  else
    indcolscaling = 1.;
  //indcolscaling = 0;
  float centerofmassx = 0;
  float centerofmassy = 0;
  float denom = 0; 
  float intercept=0.;
  float slope = 0.;
  //this array keeps track of the hits that have already been associated with a line. 
  int xMax = 0;
  int yMax = 0;
  int yClearStart;
  int yClearEnd;
  int xClearStart;
  int xClearEnd;
  int maxCell = 0;
  float rho;
  float theta; 
  int accDx(0), accDy(0);
  float pCornerMin[2];
  float pCornerMax[2];
  //float pCorner0[2];
  //float pCorner1[2];
  //bool newChannel = false;
  //unsigned int lastChannel;

  unsigned int fMaxWire = 0;
  int iMaxWire = 0;
  unsigned int fMinWire = 99999999;
  int iMinWire = -1;




  /// Outline of PPHT, J. Matas et. al. 
  /// ---------------------------------------
  /// 
  ///LOOP over hits, picking a random one
  ///  Enter the point into the accumulator
  ///  IF it is already in the accumulator or part of a line, skip it
  ///  Store it in a vector of points that have been chosen
  ///
  ///  Find max value in accumulator; IF above threshold, create a line
  ///    Subtract points in line from accumulator
  ///    
  ///
  ///END LOOP over hits, picking a random one
  ///

  ///Init specifies the size of the two-dimensional accumulator 
  ///(based on the arguments, number of wires and number of time samples). 
  c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);
  /// Adds all of the hits to the accumulator
  //mf::LogInfo("HoughBaseAlg") << "Beginning PPHT";

  c.GetAccumSize(accDy, accDx);

  // Define the prototrack object
  protoTrack protoTrackToLoad;


  // The number of lines we've found
  unsigned int nLinesFound = 0;
  std::vector<unsigned int> accumPoints(hits.size(),0);
  int nAccum = 0;

  /// count is how many points are left to randomly insert
  int count = 0;
  for(auto fpointId_to_clusterIdItr = fpointId_to_clusterId->begin(); fpointId_to_clusterIdItr != fpointId_to_clusterId->end();fpointId_to_clusterIdItr++)
    if(*fpointId_to_clusterIdItr == clusterId)
      count++;

  unsigned int randInd;

  /// Get the random number generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine & engine = rng -> getEngine();
  CLHEP::RandFlat flat(engine);
  TStopwatch w;
  //float timeTotal = 0;

  for( ; count > 0; ){

    /// The random hit we are examining
    ///unsigned int randInd = rand() % hits.size();
    randInd = (unsigned int)(flat.fire()*hits.size());

    //std::cout << count << " " << randInd << std::endl;
    
      /// If the point isn't in the current fuzzy cluster, skip it
    if(fpointId_to_clusterId->at(randInd) != clusterId)
      continue;

    --count;

    /// Skip if it's already in a line
    if(skip[randInd]==1){
      continue;
    }

    /// If we have already accumulated the point, skip it
    if(accumPoints[randInd])
      continue;
    accumPoints[randInd]=1;
   
    /// zeroes out the neighborhood of all previous lines  
    for(auto listofxmaxItr = listofxmax.begin(); listofxmaxItr != listofxmax.end(); ++listofxmaxItr) {
      yClearStart = listofymax[listofxmaxItr-listofxmax.begin()] - fRhoZeroOutRange;
      if (yClearStart < 0) yClearStart = 0;
      
      yClearEnd = listofymax[listofxmaxItr-listofxmax.begin()] + fRhoZeroOutRange;
      if (yClearEnd >= accDy) yClearEnd = accDy - 1;
      
      xClearStart = *listofxmaxItr - fThetaZeroOutRange;
      if (xClearStart < 0) xClearStart = 0;
      
      xClearEnd = *listofxmaxItr + fThetaZeroOutRange;
      if (xClearEnd >= accDx) xClearEnd = accDx - 1;
      
      for (y = yClearStart; y <= yClearEnd; ++y){
        for (x = xClearStart; x <= xClearEnd; ++x){
          c.SetCell(y,x,0);
        }
      }
    }/// end loop over size of listxmax
  
    /// Find the weightiest cell in the accumulator.
     uint32_t channel = hits[randInd]->Channel();
    wireMax = hits[randInd]->WireID().Wire;

    /// Add the randomly selected point to the accumulator
    //w.Start();
    std::vector<int> max = c.AddPointReturnMax(wireMax, (int)(hits[randInd]->PeakTime()));
    maxCell = max[0];
    xMax    = max[1];
    yMax    = max[2];
    //w.Stop();
    //std::cout << "Real Time: " << w.RealTime() << std::endl;
    //timeTotal += w.CpuTime();
    ++nAccum; 

    //mf::LogVerbatim("HoughBaseAlg") << "cout: " << count << " maxCell: " << maxCell << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "xMax: " << xMax << " yMax: " << yMax << std::endl;

    /// The threshold calculation, see http://www.via.cornell.edu/ece547/projects/Hough/webpage/statistics.html
    /// accDx is the number of rho bins,m_rowLength
    //TF1 *threshGaus = new TF1("threshGaus","(1/([0]*std::sqrt(2*TMath::Pi())))*exp(-0.5*pow(((x-[1])/[0]),2))");
    //double sigma = std::sqrt(((double)nAccum/(double)accDx)*(1-1/(double)accDx));
    //double mean = (double)nAccum/(double)accDx;
    //threshGaus->SetParameter(0,sigma);
    //threshGaus->SetParameter(1,mean);
    //mf::LogVerbatim("HoughBaseAlg") << "threshGaus mean: " << mean << " sigma: " << sigma << " accDx: " << accDx << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "nAccum: " << nAccum << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "threshGaus integral range: " << mean-2*sigma << " to " << maxCell << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "threshGaus integral: " << threshGaus->Integral(mean-2*sigma,maxCell) << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "threshGaus integral: " << threshGaus->Integral(0,maxCell) << std::endl;


    /// The threshold calculation using a Poisson distribution instead
    //double poisProbSum = 0;
    //for(int j = 0; j <= maxCell; j++){
    //double poisProb = TMath::Poisson(j,mean);
    //poisProbSum+=poisProb;
    //mf::LogVerbatim("HoughBaseAlg") << "Poisson: " << poisProb << std::endl;
    //}
    //mf::LogVerbatim("HoughBaseAlg") << "Poisson prob sum: " << poisProbSum << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "Probability it is higher: " << 1-poisProbSum << std::endl;

    // Continue if the probability of finding a point, (1-poisProbSum) is the probability of finding a 
    // value of maxCell higher than what it currently is
    //if( (1-poisProbSum) > 1e-13)
    //continue;


    // The threshold calculation using a Binomial distribution instead
    double binomProbSum = TMath::BinomialI(1/(double)accDx,nAccum,maxCell);
    //std::cout << "nAccum: " << nAccum << std::endl;
    //std::cout << "maxCell: " << maxCell << std::endl;
    //std::cout << "BinomialI: " << binomProbSum << std::endl;
    //std::cout << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "Probability it is higher: " << 1-binomProbSum << std::endl;
    //Continue if the probability of finding a point, (1-poisProbSum) is the probability of finding a
    //value of maxCell higher than what it currently is
    if( binomProbSum > 1e-21)
      continue;



    /// Continue if the biggest maximum for the randomly selected point is smaller than fMinHits
    //if (maxCell < fMinHits) 
      //continue;


    /// Find the center of mass of the 3x3 cell system (with maxCell at the center). 
    denom = centerofmassx = centerofmassy = 0;
  
    if(xMax > 0 && yMax > 0 && xMax + 1 < accDx && yMax + 1 < accDy){  
      for(int i = -1; i < 2; ++i){
        for(int j = -1; j < 2; ++j){
          denom += c.GetCell(yMax+i,xMax+j);
          centerofmassx += j*c.GetCell(yMax+i,xMax+j);
          centerofmassy += i*c.GetCell(yMax+i,xMax+j);
        }
      }
      centerofmassx /= denom;
      centerofmassy /= denom;      
    }
    else  centerofmassx = centerofmassy = 0;
  
    ///fill the list of cells that have already been found
    listofxmax.push_back(xMax);
    listofymax.push_back(yMax);

    // Find the lines equation
    c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
    slope = -1./tan(theta);    
    intercept = (rho/sin(theta));
    float distance;

    if(!std::isinf(slope) && !std::isnan(slope)){
      channelHolder.clear();
      sequenceHolder.clear();
      fMaxWire = 0;
      iMaxWire = 0;
      fMinWire = 99999999;
      iMinWire = -1;
      hitsTemp.clear();
      for(auto hitsItr = hits.cbegin(); hitsItr != hits.cend(); ++hitsItr){
        wire = (*hitsItr)->WireID().Wire;
        if(fpointId_to_clusterId->at(hitsItr - hits.begin()) != clusterId)
          continue;
        channel = (*hitsItr)->Channel();	
        distance = (std::abs((*hitsItr)->PeakTime()-slope*(float)((*hitsItr)->WireID().Wire)-intercept)/(std::sqrt(xyScale[(*hitsItr)->WireID().Plane]*xyScale[(*hitsItr)->WireID().Plane]*slope*slope+1)));
        if(distance < fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling && skip[hitsItr-hits.begin()]!=1){
          hitsTemp.push_back(hitsItr-hits.begin());
          sequenceHolder.push_back(wire);
          channelHolder.push_back(channel);
        }
      }/// end iterator over hits
   
      if(hitsTemp.size() < 5) continue;
      currentHits.clear();  
      lastHits.clear();
      int j; 
      currentHits.push_back(0);
      for(auto sequenceHolderItr = sequenceHolder.begin(); sequenceHolderItr+1 != sequenceHolder.end(); ++sequenceHolderItr) {
        j = 1;
        while((chanFilt.BadChannel(sequenceHolderItr-sequenceHolder.begin()+j)) == true) j++;
        if(sequenceHolder[sequenceHolderItr-sequenceHolder.begin()+1]-sequenceHolder[sequenceHolderItr-sequenceHolder.begin()] <= j + fMissedHits) currentHits.push_back(sequenceHolderItr-sequenceHolder.begin()+1);
        else if(currentHits.size() > lastHits.size()) {
          lastHits = currentHits;
          currentHits.clear();
        }
        else currentHits.clear();
      } 
      if(currentHits.size() > lastHits.size()) lastHits = currentHits;

      clusterHits.clear();    

      //if(lastHits.size() < 5) continue;
      if(lastHits.size() < (unsigned int)fMinHits) continue;







      //// Check if lastHits has hits with big gaps in it
      //// lastHits[i] is ordered in increasing channel and then increasing peak time,
      //// as a consequence, if the line has a negative slope and there are multiple hits in the line for a channel,
      //// we have to go back to the first hit (in terms of lastHits[i]) of that channel to find the distance
      //// between hits
      ////std::cout << "New line" << std::endl;
      ////std::cout << "slope: " << slope << std::endl;
      //int missedHits=0;
      //int lastHitsChannel = 0;
      //fMaxWire = 0;
      //iMaxWire = 0;
      //fMinWire = 99999999;
      //iMinWire = -1;
      //newChannel = false;
      //lastChannel = hits[hitsTemp[lastHits[0]]]->Wire()->RawDigit()->Channel();
      //for(auto lastHitsItr = lastHits.begin(); lastHitsItr != lastHits.end()-1; ++lastHitsItr) {

        //newChannel = false;
        //if(slope < 0){
          //if(hits[hitsTemp[*lastHitsItr+1]]->Wire()->RawDigit()->Channel() != lastChannel){
            //newChannel = true;
          //}           
        //}
        //if(slope > 0 || !newChannel){

          ////std::cout << hits[hitsTemp[lastHits[i]]]->Wire()->RawDigit()->Channel() << " " << ((hits[hitsTemp[lastHits[i]]]->StartTime()+hits[hitsTemp[lastHits[i]]]->EndTime())/2.) << std::endl;
          //pCorner0[0] = (hits[hitsTemp[*lastHitsItr]]->Wire()->RawDigit()->Channel())*wire_dist;
          //pCorner0[1] = ((hits[hitsTemp[*lastHitsItr]]->StartTime()+hits[hitsTemp[*lastHitsItr]]->EndTime())/2.)*tickToDist;
          //pCorner1[0] = (hits[hitsTemp[*lastHitsItr+1]]->Wire()->RawDigit()->Channel())*wire_dist;
          //pCorner1[1] = ((hits[hitsTemp[*lastHitsItr+1]]->StartTime()+hits[hitsTemp[*lastHitsItr+1]]->EndTime())/2.)*tickToDist;
          ////std::cout << std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) << std::endl;
          ////if(std::sqrt(pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) > fMissedHitsDistance)
          //if((pCorner0[0]-pCorner1[0])*(pCorner0[0]-pCorner1[0]) + (pCorner0[1]-pCorner1[1])*(pCorner0[1]-pCorner1[1]) > fMissedHitsDistance*fMissedHitsDistance)
            //missedHits++;
        ////} else if (slope < 0 && newChannel && nHitsInChannel > 1){
        //} else if (slope < 0 && newChannel){
          ////std::cout << lastHitsChannel << " " << lastHits[i+1] << " " << lastChannel << std::endl;
          ////std::cout << hits[hitsTemp[lastHits[lastHitsChannel]]]->Wire()->RawDigit()->Channel() << " " << ((hits[hitsTemp[lastHits[lastHitsChannel]]]->StartTime()+hits[hitsTemp[lastHits[lastHitsChannel]]]->EndTime())/2.) << std::endl;
          //pCorner0[0] = (hits[hitsTemp[lastHits[lastHitsChannel]]]->Wire()->RawDigit()->Channel())*wire_dist;
          //pCorner0[1] = ((hits[hitsTemp[lastHits[lastHitsChannel]]]->StartTime()+hits[hitsTemp[lastHits[lastHitsChannel]]]->EndTime())/2.)*tickToDist;
          //pCorner1[0] = (hits[hitsTemp[*lastHitsItr+1]]->Wire()->RawDigit()->Channel())*wire_dist;
          //pCorner1[1] = ((hits[hitsTemp[*lastHitsItr+1]]->StartTime()+hits[hitsTemp[*lastHitsItr+1]]->EndTime())/2.)*tickToDist;
          ////std::cout << std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) << std::endl;
          ////if(std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) > fMissedHitsDistance             )
          //if((pCorner0[0]-pCorner1[0])*(pCorner0[0]-pCorner1[0]) + (pCorner0[1]-pCorner1[1])*(pCorner0[1]-pCorner1[1]) > fMissedHitsDistance*fMissedHitsDistance)
            //missedHits++;
          //lastChannel=hits[hitsTemp[*lastHitsItr+1]]->Wire()->RawDigit()->Channel();
          //lastHitsChannel=lastHitsItr-lastHits.begin()+1;
        //}
      //}
      ////std::cout << "missedHits " << missedHits << std::endl;
      ////std::cout << "lastHits.size() " << lastHits.size() << std::endl;
      ////std::cout << "missedHits/lastHits.size() " << (double)missedHits/((double)lastHits.size()-1) << std::endl;
      //if((float)missedHits/((float)lastHits.size()-1) > fMissedHitsToLineSize)
        //continue;







      if(std::abs(slope)>fMaxSlope ){
        continue;
      }

      //std::cout << std::endl;
      //std::cout  << "Found line!" << std::endl
                                       //<< "Slope: " << slope << std::endl
                                       //<< "Intercept: " << intercept << std::endl
                                       //<< "Number of hits: " << lastHits.size() << std::endl
                                       //<< "Wire: " << fMinWire << " Peak time: " 
                                       //<< hits[iMinWire]->PeakTime() << std::endl
                                       //<< "Wire: " << fMaxWire << " Peak time: " 
                                       //<< hits[iMaxWire]->PeakTime() << std::endl;


      // Add new line to list of lines
      float totalQ = 0;
      std::vector< art::Ptr<recob::Hit> > lineHits;
      nClustersTemp++;
      ///std::cout << "nClusters: " << *nClusters << std::endl;
      for(auto lastHitsItr = lastHits.begin(); lastHitsItr != lastHits.end(); ++lastHitsItr) {
        fpointId_to_clusterId->at(hitsTemp[(*lastHitsItr)]) = nClustersTemp-1;
        //clusterHits.push_back(hits[hitsTemp[(*lastHitsItr)]]);
        //totalQ += clusterHits.back()->Charge();
        totalQ += hits[hitsTemp[(*lastHitsItr)]]->Charge();
        wire = hits[hitsTemp[(*lastHitsItr)]]->WireID().Wire;

        if(!accumPoints[hitsTemp[(*lastHitsItr)]])
          count--;
        
        skip[hitsTemp[(*lastHitsItr)]]=1;

        lineHits.push_back(hits[hitsTemp[(*lastHitsItr)]]);


        /// Subtract points from the accumulator that have already been used
        if(accumPoints[hitsTemp[(*lastHitsItr)]]) 
          c.SubtractPoint(wire, (int)(hits[hitsTemp[(*lastHitsItr)]]->PeakTime()));
        
        if(wire < fMinWire){
          fMinWire = wire;
          iMinWire = hitsTemp[(*lastHitsItr)];
        }
        if(wire > fMaxWire){
          fMaxWire = wire;
          iMaxWire = hitsTemp[(*lastHitsItr)];
        }
      }
      int pnum = hits[iMinWire]->WireID().Plane;    
      pCornerMin[0] = (hits[iMinWire]->WireID().Wire)*wire_pitch[pnum];
      pCornerMin[1] = ((hits[iMinWire]->StartTime()+hits[iMinWire]->EndTime())/2.)*tickToDist;
      pCornerMax[0] = (hits[iMaxWire]->WireID().Wire)*wire_pitch[pnum];
      pCornerMax[1] = ((hits[iMaxWire]->StartTime()+hits[iMaxWire]->EndTime())/2.)*tickToDist;

      ///std::cout << std::endl;
      ///std::cout << "pCornerMin[0]: " << pCornerMin[0] << " pCornerMin[1]: " << pCornerMin[1] << std::endl;
      ///std::cout << "pCornerMax[0]: " << pCornerMax[0] << " pCornerMax[1]: " << pCornerMax[1] << std::endl;
      protoTrackToLoad.Init(nClustersTemp-1,
	    pnum,
            slope,
            intercept,
            totalQ,
            pCornerMin[0],
            pCornerMin[1],
            pCornerMax[0],
            pCornerMax[1],
            iMinWire,
            iMaxWire,
            fMinWire,
            fMaxWire,
            lineHits);
      linesFound->push_back(protoTrackToLoad);
       
    }/// end if !std::isnan


    nLinesFound++;

    if(nLinesFound>(unsigned int)fMaxLines)
      break;

  }/// end loop over hits

  //std::cout << "Average cpu time: " << timeTotal/nAccum << std::endl;
  //std::cout << "Total cpu time: " << timeTotal << std::endl;





  lastHits.clear();

  listofxmax.clear();
  listofymax.clear();

  // saves a bitmap image of the accumulator (useful for debugging), 
  // with scaling based on the maximum cell value
  if(fSaveAccumulator){   
    unsigned char *outPix = new unsigned char [accDx*accDy];
    //finds the maximum cell in the accumulator for image scaling
    int cell, pix = 0, maxCell = 0;
    for (int y = 0; y < accDy; ++y){ 
      for (int x = 0; x < accDx; ++x){
        cell = c.GetCell(y,x);
        if (cell > maxCell) maxCell = cell;
      }
    }
    for (int y = 0; y < accDy; ++y){
      for (int x = 0; x < accDx; ++x){ 
        //scales the pixel weights based on the maximum cell value     
        if(maxCell > 0)
          pix = (int)((1500*c.GetCell(y,x))/maxCell);
        outPix[y*accDx + x] = pix;
      }
    }
                
    HLSSaveBMPFile("houghaccum.bmp", outPix, accDx, accDy);
    delete [] outPix;
  }// end if saving accumulator

  *nClusters=nClustersTemp;

  return 1; 
}


//------------------------------------------------------------------------------
cluster::HoughTransform::~HoughTransform()
{
}

//------------------------------------------------------------------------------
// returns a vector<int> where the first is the overall maximum,
// the second is the max x value, and the third is the max y value.
inline std::vector<int> cluster::HoughTransform::AddPointReturnMax(int x, 
								   int y)
{
  std::vector<int> max(3,0);
  if (x>m_dx || y>m_dy || x<0.0 || y<0.0)
    return max;
  return DoAddPointReturnMax(x, y);
}



//------------------------------------------------------------------------------
inline bool cluster::HoughTransform::SubtractPoint(int x, int y)
{
  if (x>m_dx || y>m_dy || x<0.0 || y<0.0)
    return false;
  return DoSubtractPoint(x, y);
}


//------------------------------------------------------------------------------
void cluster::HoughTransform::Init(int dx, 
				   int dy, 
				   int rhores,
				   int numACells)
{
  m_numAngleCells=numACells;
  m_rhoResolutionFactor = rhores;
  m_accum.clear();
  //m_accum.resize(m_numAngleCells);
  m_numAccumulated = 0;   
  //   m_cosTable.clear();
  //   m_sinTable.clear();
  //m_cosTable.resize(m_numAngleCells);
  //m_sinTable.resize(m_numAngleCells);
  //if (dx == m_dx && dy == m_dy)
  //return;
  m_dx = dx;
  m_dy = dy;
  m_rowLength = (int)(m_rhoResolutionFactor*2 * std::sqrt(dx*dx + dy*dy));
  m_accum.resize(m_numAngleCells);
  //for(int i = 0; i < m_numAngleCells; i++)
    //m_accum[i].resize((int)(m_rowLength));

  //int angleIndex;
  //float a, angleStep = TMath::Pi()/m_numAngleCells;
  //for (a=0.0, angleIndex = 0; angleIndex < m_numAngleCells; ++angleIndex){
    //m_cosTable[angleIndex] = (cos(a));
    //m_sinTable[angleIndex] = (sin(a));
    //a += angleStep;
  //}
}


//------------------------------------------------------------------------------
int cluster::HoughTransform::GetMax(int &xmax, 
				    int &ymax)
{
  int maxVal = -1;
  for(unsigned int i = 0; i < m_accum.size(); i++){
    for(auto rhoIter=m_accum[i].begin(); rhoIter!=m_accum[i].end(); ++rhoIter){
      if((*rhoIter).second > maxVal) {
      //if((*rhoIter) > maxVal) {
        maxVal = (*rhoIter).second;
	//maxVal = (*rhoIter);
	xmax = i;
        ymax = (*rhoIter).first;
	//ymax = rhoIter-m_accum[i].begin();
      }
    }
  }

  return maxVal;
}

//------------------------------------------------------------------------------
// returns a vector<int> where the first is the overall maximum,
// the second is the max x value, and the third is the max y value.
inline std::vector<int> cluster::HoughTransform::DoAddPointReturnMax(int x, 
								     int y)
{
  std::vector<int> max(3,-1);

  int distCenter = (int)(m_rowLength/2.);
  
  // prime the lastDist variable so our linear fill works below
  int lastDist = (int)(distCenter+(m_rhoResolutionFactor*x));

  //int max_val = minHits-1;
  int max_val = 0;
  int dist;
  int val;
  int stepDir;
  int cell;

  //std::pair<std::map<int,int>::iterator,bool> m_accumA;
  std::map<int,int>::iterator m_accumA;
  std::map<int,int>::iterator m_accumALast;

  // loop through all angles a from 0 to 180 degrees
  float angleStep = PI/m_numAngleCells;
  float angleStepInverse = m_numAngleCells/PI;
  float angleStepInt;
  //for (float a = angleStep; a < TMath::Pi(); a+=angleStep){
  for (float a = angleStep; a < PI; a+=angleStep){
   
    angleStepInt = (a*angleStepInverse);
    // Calculate the basic line equation dist = cos(a)*x + sin(a)*y.
    // Shift to center of row to cover negative values
    dist = (int)(distCenter+(m_rhoResolutionFactor*(cos(a)*x + sin(a)*y)));
    // sanity check to make sure we stay within our row
    //if (dist >= 0 && dist<m_rowLength){
      if(lastDist==dist){
        m_accumA = m_accum[angleStepInt].insert(m_accum[angleStepInt].end(),std::pair<int,int>(lastDist,0));
        //val = m_accumA.first->second++;
        val = ++(m_accumA->second);
        
        if( max_val < val){
	  max_val = val;
          max[0] = val;
          max[1] = lastDist;
          max[2] = (angleStepInt);
        }
        m_accumALast = m_accumA;
        //mf::LogVerbatim("HoughBaseAlg") << "First, a: " << a << " lastDist: " << lastDist << std::endl;
      }
      else{
        // fill in all values in row a, not just a single cell
        m_accumALast = m_accum[angleStepInt].begin();
        stepDir = dist>lastDist ? 1 : -1;
        for (cell=lastDist; cell!=dist; cell+=stepDir){  
          m_accumA = m_accum[angleStepInt].insert(m_accumALast,std::pair<int,int>(cell,0));
          //val = m_accumA.first->second++;
          val = ++(m_accumA->second);
          //// Note, m_accum is a vector of associative containers, "a" calls the vector element, "cell" is the container key, and the ++ iterates the value correspoding to the key
          if(max_val < val){
	    max_val = val;
            max[0] = val;
            max[1] = cell;
            max[2] = (angleStepInt);
          }
          m_accumALast = m_accumA;
        }      
      }
    //}
    lastDist = dist;
  }
  m_numAccumulated++;


  //mf::LogVerbatim("HoughBaseAlg") << "Add point says xmax: " << *xmax << " ymax: " << *ymax << std::endl;

  return max;
}



//------------------------------------------------------------------------------
inline bool cluster::HoughTransform::DoSubtractPoint(int x, int y)
{
  int distCenter = (int)(m_rowLength/2.);
 
  // prime the lastDist variable so our linear fill works below
  int lastDist = (int)(distCenter+(m_rhoResolutionFactor*x));
  int dist;
  int stepDir;
  // loop through all angles a from 0 to 180 degrees
  float angleStep = PI/m_numAngleCells;
  for (float a = angleStep; a < PI; a+=angleStep){
    // Calculate the basic line equation dist = cos(a)*x + sin(a)*y.
    // Shift to center of row to cover negative values
    dist = (int)(distCenter+(m_rhoResolutionFactor*((cos(a)*x + sin(a)*y))));
    // sanity check to make sure we stay within our row
      if(lastDist==dist)
	m_accum[(int)(a/angleStep)][lastDist]--;
      else{
	// fill in all values in row a, not just a single cell
	stepDir = dist>lastDist ? 1 : -1;
	for (int cell=lastDist; cell!=dist; cell+=stepDir){   
	  m_accum[(int)(a/angleStep)][cell]--;//maybe add weight of hit here?
          // Note, m_accum is a vector of associative containers, "a" calls the vector element, "cell" is the container key, and the -- iterates the value correspoding to the key
	}      
      }
    lastDist = dist;
  }
  m_numAccumulated--;

  return true;
}

//------------------------------------------------------------------------------
//this method saves a BMP image of the Hough Accumulator, which can be viewed with gimp
void cluster::HoughBaseAlg::HLSSaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy)
{
  ofstream bmpFile(fileName, std::ios::binary);
  bmpFile.write("B", 1);
  bmpFile.write("M", 1);
  int bitsOffset = 54 +256*4; 
  int size = bitsOffset + dx*dy; //header plus 256 entry LUT plus pixels
  bmpFile.write((const char *)&size, 4);
  int reserved = 0;
  bmpFile.write((const char *)&reserved, 4);
  bmpFile.write((const char *)&bitsOffset, 4);
  int bmiSize = 40;
  bmpFile.write((const char *)&bmiSize, 4);
  bmpFile.write((const char *)&dx, 4);
  bmpFile.write((const char *)&dy, 4);
  short planes = 1;
  bmpFile.write((const char *)&planes, 2);
  short bitCount = 8;
  bmpFile.write((const char *)&bitCount, 2);
  int i, temp = 0;
  for (i=0; i<6; i++)
    bmpFile.write((const char *)&temp, 4);  // zero out optional color info
  // write a linear LUT
  char lutEntry[4]; // blue,green,red
  lutEntry[3] = 0;  // reserved part
  for (i=0; i<256; i++)
    {
      lutEntry[0] = lutEntry[1] = lutEntry[2] = i;
      bmpFile.write(lutEntry, sizeof lutEntry);
    }
  // write the actual pixels
  bmpFile.write((const char *)pix, dx*dy);
}
 

//------------------------------------------------------------------------------
size_t cluster::HoughBaseAlg::FastTransform(const std::vector<art::Ptr<recob::Cluster> >         & clusIn,
					    std::vector<recob::Cluster>                    & ccol,  
					    std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
					    art::Event                                const& evt,
					    std::string                               const& label)
{
  std::vector<int> skip;  

  art::FindManyP<recob::Hit> fmh(clusIn, evt, label);

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  filter::ChannelFilter chanFilt;
  HoughTransform c;

  // Get the random number generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine & engine = rng -> getEngine();
  CLHEP::RandFlat flat(engine);

  std::vector< art::Ptr<recob::Hit> > hit;

  for(auto view : geom->Views() ){

    LOG_DEBUG("HoughBaseAlg") << "Analyzing view " << view;

    art::PtrVector<recob::Cluster>::const_iterator clusterIter = clusIn.begin();
    int clusterID = 0;//the unique ID of the cluster
    
    size_t cinctr = 0;
    while(clusterIter != clusIn.end()) {

      LOG_DEBUG("HoughBaseAlg") << "Analyzing cinctr=" << cinctr << " which is in view " << (*clusterIter)->View();

      hit.clear();
      if(fPerCluster){
	if((*clusterIter)->View() == view) hit = fmh.at(cinctr);
      }
      else{   
	while(clusterIter != clusIn.end()){
	  if( (*clusterIter)->View() == view ){
	    
	    hit = fmh.at(cinctr);
	  }// end if cluster is in correct view
	  clusterIter++;
	  ++cinctr;
	}//end loop over clusters
      }//end if not fPerCluster
      if(hit.size() == 0){ 
	if(fPerCluster){
	  clusterIter++;
	  ++cinctr;
	}
	continue;
      }
      
      LOG_DEBUG("HoughBaseAlg") << "We have all the hits..." << hit.size();
      
      /*
      //factor to make x and y scale the same units
      double xyScale  = .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
      xyScale        *= detprop->SamplingRate()/geom->WirePitch(0,1,p,t,cs);
      
      int x, y;
      int dx = geom->Cryostat(cs).TPC(t).Plane(p).Nwires();//number of wires 
      int dy = hit[0]->Wire()->NSignal();//number of time samples. 
      skip.clear();
      skip.resize(hit.size());
      std::vector<int> listofxmax;
      std::vector<int> listofymax;  
      std::vector<int> hitTemp;        //indecies ofcandidate hits
      std::vector<int> sequenceHolder; //channels of hits in list
      std::vector<int> currentHits;    //working vector of hits 
      std::vector<int> lastHits;       //best list of hits
      art::PtrVector<recob::Hit> clusterHits;
      double indcolscaling = 0.;       //a parameter to account for the different 
      //characteristic hit width of induction and collection plane
      double centerofmassx = 0;
      double centerofmassy = 0;
      double denom = 0; 
      double intercept=0.;
      double slope = 0.;
      //this array keeps track of the hits that have already been associated with a line. 
      int xMax = 0;
      int yMax = 0;
      double rho;
      double theta; 
      int accDx(0), accDy(0);
      
      
      //Init specifies the size of the two-dimensional accumulator 
      //(based on the arguments, number of wires and number of time samples). 
      //adds all of the hits (that have not yet been associated with a line) to the accumulator
      c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);
      
      
      // count is how many points are left to randomly insert
      unsigned int count = hit.size();
      std::vector<unsigned int> accumPoints;
      accumPoints.resize(hit.size());
      int nAccum = 0;
      unsigned int nLinesFound = 0;
      
      for( ; count > 0; count--){
      
      
      // The random hit we are examining
      unsigned int randInd = (unsigned int)(flat.fire()*hit.size());
      
      // Skip if it's already in a line
      if(skip[randInd]==1)
      continue;
      
      // If we have already accumulated the point, skip it
      if(accumPoints[randInd])
      continue;
      accumPoints[randInd]=1;
      
      // zeroes out the neighborhood of all previous lines  
      for(unsigned int i = 0; i < listofxmax.size(); ++i){
      int yClearStart = listofymax[i] - fRhoZeroOutRange;
      if (yClearStart < 0) yClearStart = 0;
      
      int yClearEnd = listofymax[i] + fRhoZeroOutRange;
      if (yClearEnd >= accDy) yClearEnd = accDy - 1;
      
      int xClearStart = listofxmax[i] - fThetaZeroOutRange;
      if (xClearStart < 0) xClearStart = 0;
      
      int xClearEnd = listofxmax[i] + fThetaZeroOutRange;
      if (xClearEnd >= accDx) xClearEnd = accDx - 1;
      
      for (y = yClearStart; y <= yClearEnd; ++y){
      for (x = xClearStart; x <= xClearEnd; ++x){
      c.SetCell(y,x,0);
      }
      }
      }// end loop over size of listxmax
      
      //find the weightiest cell in the accumulator.
      int maxCell = 0;
      unsigned int wireMax = hit[randInd]->WireID().Wire;
      xMax = 0;
      yMax = 0;
      
      // Add the randomly selected point to the accumulator
      maxCell = c.AddPointReturnMax(wireMax, (int)(hit[randInd]->PeakTime()), &yMax, &xMax, fMinHits);
      nAccum++; 
      
      // Continue if the biggest maximum for the randomly selected point is smaller than fMinHits
      if (maxCell < fMinHits) 
      continue;
      
      // Find the center of mass of the 3x3 cell system (with maxCell at the center). 
      denom = centerofmassx = centerofmassy = 0;
      
      if(xMax > 0 && yMax > 0 && xMax + 1 < accDx && yMax + 1 < accDy){  
      for(int i = -1; i < 2; ++i){
      for(int j = -1; j < 2; ++j){
      denom += c.GetCell(yMax+i,xMax+j);
      centerofmassx += j*c.GetCell(yMax+i,xMax+j);
      centerofmassy += i*c.GetCell(yMax+i,xMax+j);
      }
      }
      centerofmassx /= denom;
      centerofmassy /= denom;      
      }
      else  centerofmassx = centerofmassy = 0;
      
      //fill the list of cells that have already been found
      listofxmax.push_back(xMax);
      listofymax.push_back(yMax);
      
      // Find the lines equation
      c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
      //c.GetEquation(yMax, xMax, rho, theta);
      slope = -1./tan(theta);    
      intercept = (rho/sin(theta));
      //mf::LogVerbatim("HoughBaseAlg") << std::endl;
      //mf::LogVerbatim("HoughBaseAlg") << "slope: " << slope << " intercept: " << intercept << std::endl; 
      //mf::LogInfo("HoughBaseAlg") << "slope: " << slope << " intercept: " << intercept;
      double distance;
      /// \todo: the collection plane's characteristic hit width's are, 
      /// \todo: on average, about 5 time samples wider than the induction plane's. 
      /// \todo: this is hard-coded for now.
      if(sigt == geo::kInduction)
      indcolscaling = 5.;
      else
      indcolscaling = 0.;
      // What is this?   
      indcolscaling = 0;
      
      if(!std::isinf(slope) && !std::isnan(slope)){
      sequenceHolder.clear();
      hitTemp.clear();
      for(size_t i = 0; i < hit.size(); ++i){
      distance = (TMath::Abs(hit[i]->PeakTime()-slope*(double)(hit[i]->WireID().Wire)-intercept)/(std::sqrt(pow(xyScale*slope,2)+1)));
      
      if(distance < fMaxDistance+((hit[i]->EndTime()-hit[i]->StartTime())/2.)+indcolscaling  && skip[i]!=1){
      hitTemp.push_back(i);
      sequenceHolder.push_back(hit[i]->Channel());
      }
      
      }// end loop over hits
      
      if(hitTemp.size() < 2) continue;
      currentHits.clear();  
      lastHits.clear();
      int j; 
      currentHits.push_back(0);
      for(size_t i = 0; i + 1 < sequenceHolder.size(); ++i){  
      j = 1;
      while((chanFilt.BadChannel(sequenceHolder[i]+j)) == true) j++;
      if(sequenceHolder[i+1]-sequenceHolder[i] <= j + fMissedHits) currentHits.push_back(i+1);
      else if(currentHits.size() > lastHits.size()) {
      lastHits = currentHits;
      currentHits.clear();
      }
      else currentHits.clear();
      }
      
      
      if(currentHits.size() > lastHits.size()) lastHits = currentHits;
      
      
      
      
      // Check if lastHits has hits with big gaps in it
      uint32_t     channel = hit[0]->Wire()->RawDigit()->Channel();
      double wirePitch = geom->WirePitch(geom->View(channel));
      double wire_dist = wirePitch;
      double tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
      tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
      //std::cout << "New line" << std::endl;
      int missedHits=0;
      for(size_t i = 0; i < lastHits.size()-1; ++i) {
      //std::cout << hit[hitTemp[lastHits[i]]]->Wire()->RawDigit()->Channel() << std::endl;
      double pCorner0[2];
      pCorner0[0] = (hit[hitTemp[lastHits[i]]]->Wire()->RawDigit()->Channel())*wire_dist;
      pCorner0[1] = ((hit[hitTemp[lastHits[i]]]->StartTime()+hit[hitTemp[lastHits[i]]]->EndTime())/2.)*tickToDist;
      double pCorner1[2];
      pCorner1[0] = (hit[hitTemp[lastHits[i+1]]]->Wire()->RawDigit()->Channel())*wire_dist;
      pCorner1[1] = ((hit[hitTemp[lastHits[i+1]]]->StartTime()+hit[hitTemp[lastHits[i+1]]]->EndTime())/2.)*tickToDist;
      //std::cout << std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) << std::endl;
      if(std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) > fMissedHitsDistance             )
      missedHits++;
      }
      //std::cout << "missedHits " << missedHits << std::endl;
      //std::cout << "lastHits.size() " << lastHits.size() << std::endl;
      //std::cout << "missedHits/lastHits.size() " << (double)missedHits/(double)lastHits.size() << std::endl;
      if((double)missedHits/(double)lastHits.size() > fMissedHitsToLineSize)
      continue;
      
      
      
      
      
      clusterHits.clear();    
      double totalQ = 0.;
      if(lastHits.size() < 5) continue;
      
      
      for(size_t i = 0; i < lastHits.size(); ++i) {
      clusterHits.push_back(hit[hitTemp[lastHits[i]]]);
      totalQ += clusterHits.back()->Charge();
      skip[hitTemp[lastHits[i]]]=1;
      } 
      //protection against very steep uncorrelated hits
      if(std::abs(slope)>fMaxSlope 
      && std::abs((*clusterHits.begin())->Wire()->RawDigit()->Channel()-
      clusterHits[clusterHits.size()-1]->Wire()->RawDigit()->Channel())>=0
      )
      continue;
      
      
      
      unsigned int sw = (*clusterHits.begin())->WireID().Wire;
      unsigned int ew = (*(clusterHits.end()-1))->WireID().Wire;
      
      recob::Cluster cluster(sw, 0.,
      (*clusterHits.begin())->PeakTime(), 0.,
      ew, 0., 
      (clusterHits[clusterHits.size()-1])->PeakTime(), 0.,
      slope, 0., 
      -999., 0., 
      totalQ,
      geom->View((*clusterHits.begin())->Channel()),
      clusterID);	      
      
      ++clusterID;
      ccol.push_back(cluster);
      clusHitsOut.push_back(clusterHits);
      
      //Turn off hit sharing. T. Yang 9/14/12
      //	      //allow double assignment of first and last hits
      //	      for(size_t i = 0; i < lastHits.size(); ++i){ 
      //		if(skip[hitTemp[lastHits[i]]] ==1){
      //		  channel = hit[hitTemp[lastHits[i]]]->Wire()->RawDigit()->Channel();  
      //		  if( channel == sc || channel == ec) skip[i] = 0;
      //		}
      //	      }
      
      }// end if !std::isnan
      
      nLinesFound++;
      
      if(nLinesFound>(unsigned int)fMaxLines)
      break;
      
      
      }// end loop over hits*/
      
      std::vector<double> slopevec;std::vector<double> totalQvec;
      std::vector< art::PtrVector<recob::Hit> >   planeClusHitsOut;
      this->FastTransform(hit,planeClusHitsOut,slopevec,totalQvec );
      
      LOG_DEBUG("HoughBaseAlg") << "Made it through FastTransform" << planeClusHitsOut.size();

      for(size_t xx = 0; xx < planeClusHitsOut.size(); ++xx){
	const recob::Hit& FirstHit(**planeClusHitsOut.at(xx).begin());
	unsigned int sw = FirstHit.WireID().Wire;
	unsigned int ew = (*(planeClusHitsOut.at(xx).end()-1))->WireID().Wire;
	
	recob::Cluster cluster(sw, 0.,
			       (*planeClusHitsOut.at(xx).begin())->PeakTime(), 0.,
			       ew, 0., 
			       (planeClusHitsOut.at(xx).at(planeClusHitsOut.at(xx).size()-1))->PeakTime(), 0.,
			       slopevec.at(xx), 0., 
			       -999., 0., 
			       totalQvec.at(xx),
			       geom->View(FirstHit.Channel()),
			       clusterID,
			       FirstHit.WireID().planeID());
	
	++clusterID;
	ccol.push_back(cluster);
	clusHitsOut.push_back(planeClusHitsOut.at(xx));
      }
      
      
      hit.clear();
      //  lastHits.clear();
      if(clusterIter != clusIn.end()){
	clusterIter++;
	++cinctr;
      }
      // listofxmax.clear();
      // listofymax.clear();
    }//end loop over clusters
    
  }// end loop over views

  return ccol.size(); 

}

//------------------------------------------------------------------------------
  size_t cluster::HoughBaseAlg::FastTransform(std::vector < art::Ptr < recob::Hit > >                 & clusIn,
     	             std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut )
  {
   std::vector<double> slopevec; std::vector<double> totalQvec;
   return  FastTransform( clusIn, clusHitsOut, slopevec, totalQvec );
        
  }



//------------------------------------------------------------------------------
  size_t cluster::HoughBaseAlg::FastTransform(std::vector < art::Ptr < recob::Hit > >                 & clusIn,
     	             std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut, 
		     std::vector<double> &slopevec, std::vector<double> &totalQvec )
{
  std::vector<int> skip;  

  //art::FindManyP<recob::Hit> fmh(clusIn, evt, label);

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  filter::ChannelFilter chanFilt;
  HoughTransform c;

  // Get the random number generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine & engine = rng -> getEngine();
  CLHEP::RandFlat flat(engine);

  std::vector< art::Ptr<recob::Hit> > hit;

//   for(size_t cs = 0; cs < geom->Ncryostats(); ++cs){
//     for(size_t t = 0; t < geom->Cryostat(cs).NTPC(); ++t){
//       for(unsigned int p = 0; p < geom->Cryostat(cs).TPC(t).Nplanes(); ++p) {
// 	art::PtrVector<recob::Cluster>::const_iterator clusterIter = clusIn.begin();
// 	int clusterID = 0;//the unique ID of the cluster


	

  size_t cinctr = 0;
  //while(clusterIter != clusIn.end()) {
  hit.clear();
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator ii=clusIn.begin();ii!=clusIn.end();ii++)
    hit.push_back((*ii));   // this is new
  
  if(hit.size() == 0){ 
    if(fPerCluster){
      ++cinctr;
    }
    return -1;
  }
  
  unsigned int cs=hit.at(0)->WireID().Cryostat;
  unsigned int t=hit.at(0)->WireID().TPC;
  unsigned int p=hit.at(0)->WireID().Plane;   
  
  //	geo::View_t    view = geom->Cryostat(cs).TPC(t).Plane(p).View();
  geo::SigType_t sigt = geom->Cryostat(cs).TPC(t).Plane(p).SignalType();   
  
  // if(fPerCluster){
  //   if((*clusterIter)->View() == view) hit = fmh.at(cinctr);
  // }
  // else{   
  //  while(clusterIter != clusIn.end()){
  //    if( (*clusterIter)->View() == view ){
  
	//	hit = fmh.at(cinctr);
	 //     }// end if cluster is in correct view
	 //     clusterIter++;
	  //    ++cinctr;
	//    }//end loop over clusters
	//  }//end if not fPerCluster

  if(hit.size() == 0){ 
    if(fPerCluster){
      //      clusterIter++;
      ++cinctr;
    }
    return -1; //continue;
  }
  
  double wire_pitch[3];
  wire_pitch[0]= geom->WirePitch(0,1,0);
  wire_pitch[1]= geom->WirePitch(0,1,1);
  wire_pitch[2]= geom->WirePitch(0,1,2);

  //factor to make x and y scale the same units
  double xyScale[3];
  xyScale[2]= .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
  xyScale[0]=xyScale[2] * detprop->SamplingRate()/wire_pitch[0];
  xyScale[1]=xyScale[2] * detprop->SamplingRate()/wire_pitch[1];
  xyScale[2]*= detprop->SamplingRate()/wire_pitch[2];
  
  int x, y;
  int dx = geom->Cryostat(cs).TPC(t).Plane(p).Nwires();//number of wires 
  int dy = hit.at(0)->Wire()->NSignal();//number of time samples. 
  skip.clear();
  skip.resize(hit.size());
  std::vector<int> listofxmax;
  std::vector<int> listofymax;  
  std::vector<int> hitTemp;        //indecies ofcandidate hits
  std::vector<int> sequenceHolder; //channels of hits in list
  std::vector<int> currentHits;    //working vector of hits 
  std::vector<int> lastHits;       //best list of hits
  art::PtrVector<recob::Hit> clusterHits;
  float indcolscaling = 0.;       //a parameter to account for the different 
  //characteristic hit width of induction and collection plane
  float centerofmassx = 0;
  float centerofmassy = 0;
  float denom = 0; 
  float intercept=0.;
  float slope = 0.;
  //this array keeps track of the hits that have already been associated with a line. 
  int xMax = 0;
  int yMax = 0;
  float rho;
  float theta; 
  int accDx(0), accDy(0);
  
  
  //Init specifies the size of the two-dimensional accumulator 
  //(based on the arguments, number of wires and number of time samples). 
  //adds all of the hits (that have not yet been associated with a line) to the accumulator
  c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);
  
  // count is how many points are left to randomly insert
  unsigned int count = hit.size();
  std::vector<unsigned int> accumPoints;
  accumPoints.resize(hit.size());
  int nAccum = 0;
  unsigned int nLinesFound = 0;
  
  for( ; count > 0; count--){
    
    
    // The random hit we are examining
    unsigned int randInd = (unsigned int)(flat.fire()*hit.size());
    
    LOG_DEBUG("HoughBaseAlg") << "randInd=" << randInd << " and size is " << hit.size();

  // Skip if it's already in a line
    if(skip.at(randInd)==1)
      continue;
    
    // If we have already accumulated the point, skip it
    if(accumPoints.at(randInd))
      continue;
    accumPoints.at(randInd)=1;
    
    // zeroes out the neighborhood of all previous lines  
    for(unsigned int i = 0; i < listofxmax.size(); ++i){
      int yClearStart = listofymax.at(i) - fRhoZeroOutRange;
      if (yClearStart < 0) yClearStart = 0;
      
      int yClearEnd = listofymax.at(i) + fRhoZeroOutRange;
      if (yClearEnd >= accDy) yClearEnd = accDy - 1;
      
      int xClearStart = listofxmax.at(i) - fThetaZeroOutRange;
      if (xClearStart < 0) xClearStart = 0;
      
      int xClearEnd = listofxmax.at(i) + fThetaZeroOutRange;
      if (xClearEnd >= accDx) xClearEnd = accDx - 1;
      
      for (y = yClearStart; y <= yClearEnd; ++y){
	for (x = xClearStart; x <= xClearEnd; ++x){
	  c.SetCell(y,x,0);
	}
      }
    }// end loop over size of listxmax
  

    //find the weightiest cell in the accumulator.
    int maxCell = 0;
    unsigned int wireMax = hit.at(randInd)->WireID().Wire;
    
    // Add the randomly selected point to the accumulator
    std::vector<int> max = c.AddPointReturnMax(wireMax, 
					       (int)(hit.at(randInd)->PeakTime()));
    maxCell = max.at(0);
    xMax    = max.at(1);
    yMax    = max.at(2);
    nAccum++; 
    
    // Continue if the biggest maximum for the randomly selected point is smaller than fMinHits
    if (maxCell < fMinHits) 
      continue;
    
    // Find the center of mass of the 3x3 cell system (with maxCell at the center). 
    denom = centerofmassx = centerofmassy = 0;
    
    if(xMax > 0 && yMax > 0 && xMax + 1 < accDx && yMax + 1 < accDy){  
      for(int i = -1; i < 2; ++i){
	for(int j = -1; j < 2; ++j){
	  denom += c.GetCell(yMax+i,xMax+j);
	  centerofmassx += j*c.GetCell(yMax+i,xMax+j);
	  centerofmassy += i*c.GetCell(yMax+i,xMax+j);
	}
      }
      centerofmassx /= denom;
      centerofmassy /= denom;      
    }
    else  centerofmassx = centerofmassy = 0;
    
    //fill the list of cells that have already been found
    listofxmax.push_back(xMax);
    listofymax.push_back(yMax);
    
    // Find the lines equation
    c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
    //c.GetEquation(yMax, xMax, rho, theta);
    slope = -1./tan(theta);    
    intercept = (rho/sin(theta));
    //mf::LogVerbatim("HoughBaseAlg") << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "slope: " << slope << " intercept: " << intercept << std::endl; 
    //mf::LogInfo("HoughBaseAlg") << "slope: " << slope << " intercept: " << intercept;
    double distance;
    /// \todo: the collection plane's characteristic hit width's are, 
    /// \todo: on average, about 5 time samples wider than the induction plane's. 
    /// \todo: this is hard-coded for now.
    if(sigt == geo::kInduction)
      indcolscaling = 5.;
    else
      indcolscaling = 0.;
    // What is this?   
    //indcolscaling = 0;
    
    
    // note this doesn't work for wrapped wire planes!
    if(!std::isinf(slope) && !std::isnan(slope)){
      sequenceHolder.clear();
      hitTemp.clear();
      for(size_t i = 0; i < hit.size(); ++i){
	distance = (TMath::Abs(hit.at(i)->PeakTime()-slope*(double)(hit.at(i)->WireID().Wire)-intercept)/(std::sqrt(pow(xyScale[hit.at(i)->WireID().Plane]*slope,2)+1)));
	
	if(distance < fMaxDistance+((hit.at(i)->EndTime()-hit.at(i)->StartTime())/2.)+indcolscaling  && skip.at(i)!=1){
	  hitTemp.push_back(i);
	  sequenceHolder.push_back(hit.at(i)->Channel());
	}
	
      }// end loop over hits
      
      if(hitTemp.size() < 2) continue;
      currentHits.clear();  
      lastHits.clear();
      int j; 
      currentHits.push_back(0);
      for(size_t i = 0; i + 1 < sequenceHolder.size(); ++i){  
	j = 1;
	while((chanFilt.BadChannel(sequenceHolder.at(i)+j)) == true) j++;
	if(sequenceHolder.at(i+1)-sequenceHolder.at(i) <= j + fMissedHits) currentHits.push_back(i+1);
	else if(currentHits.size() > lastHits.size()) {
	  lastHits = currentHits;
	  currentHits.clear();
	}
	else currentHits.clear();
      }
      
      
      if(currentHits.size() > lastHits.size()) lastHits = currentHits;
      
      
      
      // Check if lastHits has hits with big gaps in it
      // lastHits[i] is ordered in increasing channel and then increasing peak time,
      // as a consequence, if the line has a negative slope and there are multiple hits in the line for a channel,
      // we have to go back to the first hit (in terms of lastHits[i]) of that channel to find the distance
      // between hits
      //std::cout << "New line" << std::endl;
      double tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
      tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
      int missedHits=0;
      int lastHitsChannel = 0;//lastHits.at(0);
      int nHitsPerChannel = 1;

      LOG_DEBUG("HoughBaseAlg") << "filling the pCorner arrays around here..."
				<< "\n but first, lastHits size is " << lastHits.size() 
				<< " and lastHitsChannel=" << lastHitsChannel;

      double pCorner0[2];
      double pCorner1[2];
      unsigned int lastChannel = hit.at(hitTemp.at(lastHits.at(0)))->Channel();

      for(size_t i = 0; i < lastHits.size()-1; ++i) {
	bool newChannel = false;
	if(slope < 0){
	  if(hit.at(hitTemp.at(lastHits.at(i+1)))->Channel() != lastChannel){
	    newChannel = true;
	  }
	  if(hit.at(hitTemp.at(lastHits.at(i+1)))->Channel() == lastChannel)
	    nHitsPerChannel++;
	}
	

	if(slope > 0 || (!newChannel && nHitsPerChannel <= 1)){

	  //std::cout << hits[hitsTemp[lastHits[i]]]->Wire()->RawDigit()->Channel() << " " << ((hits[hitsTemp[lastHits[i]]]->StartTime()+hits[hitsTemp[lastHits[i]]]->EndTime())/2.) << std::endl;
	  pCorner0[0] = (hit.at(hitTemp.at(lastHits.at(i)))->Channel())*wire_pitch[0];
	  pCorner0[1] = ((hit.at(hitTemp.at(lastHits.at(i)))->StartTime()+hit.at(hitTemp.at(lastHits.at(i)))->EndTime())/2.)*tickToDist;
	  pCorner1[0] = (hit.at(hitTemp.at(lastHits.at(i+1)))->Channel())*wire_pitch[0];
	  pCorner1[1] = ((hit.at(hitTemp.at(lastHits.at(i+1)))->StartTime()+hit.at(hitTemp.at(lastHits.at(i+1)))->EndTime())/2.)*tickToDist;
	  //std::cout << std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) << std::endl;
	  if(std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) > fMissedHitsDistance             )
	    missedHits++;
	}


	else if (slope < 0 && newChannel && nHitsPerChannel > 1){

	  //std::cout << hits[hitsTemp[lastHits[lastHitsChannel]]]->Wire()->RawDigit()->Channel() << " " << ((hits[hitsTemp[lastHits[lastHitsChannel]]]->StartTime()+hits[hitsTemp[lastHits[lastHitsChannel]]]->EndTime())/2.) << std::endl;
	  pCorner0[0] = (hit.at(hitTemp.at(lastHits.at(lastHitsChannel)))->Channel())*wire_pitch[0];
	  pCorner0[1] = ((hit.at(hitTemp.at(lastHits.at(lastHitsChannel)))->StartTime()+hit.at(hitTemp.at(lastHits.at(lastHitsChannel)))->EndTime())/2.)*tickToDist;
	  pCorner1[0] = (hit.at(hitTemp.at(lastHits.at(i+1)))->Channel())*wire_pitch[0];
	  pCorner1[1] = ((hit.at(hitTemp.at(lastHits.at(i+1)))->StartTime()+hit.at(hitTemp.at(lastHits.at(i+1)))->EndTime())/2.)*tickToDist;
	  //std::cout << std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) << std::endl;
	  if(std::sqrt( pow(pCorner0[0]-pCorner1[0],2) + pow(pCorner0[1]-pCorner1[1],2)) > fMissedHitsDistance             )
	    missedHits++;
	  lastChannel=hit.at(hitTemp.at(lastHits.at(i)))->Channel();
	  lastHitsChannel=i+1;
	  nHitsPerChannel=0;
	}
      }

      //std::cout << "missedHits " << missedHits << std::endl;
      //std::cout << "lastHits.size() " << lastHits.size() << std::endl;
      //std::cout << "missedHits/lastHits.size() " << (double)missedHits/((double)lastHits.size()-1) << std::endl;
      if((double)missedHits/((double)lastHits.size()-1) > fMissedHitsToLineSize)
	continue;
      
      
      
      
      clusterHits.clear();    
      double totalQ = 0.;
      if(lastHits.size() < 5) continue;
      
      for(size_t i = 0; i < lastHits.size(); ++i) {
	clusterHits.push_back(hit.at(hitTemp.at(lastHits.at(i))));
	totalQ += clusterHits.back()->Charge();
	skip.at(hitTemp.at(lastHits.at(i)))=1;
      } 
      //protection against very steep uncorrelated hits
      if(std::abs(slope)>fMaxSlope 
	 && std::abs((*clusterHits.begin())->Channel()-
		     clusterHits.at(clusterHits.size()-1)->Channel())>=0
	 )
	continue;
      
      
      
      // unsigned int sw = (*clusterHits.begin())->WireID().Wire;
      // unsigned int ew = (*(clusterHits.end()-1))->WireID().Wire;
      
      /*  recob::Cluster cluster(sw, 0.,
	  (*clusterHits.begin())->PeakTime(), 0.,
	  ew, 0., 
	  (clusterHits[clusterHits.size()-1])->PeakTime(), 0.,
	  slope, 0., 
	  -999., 0., 
	  totalQ,
	  geom->View((*clusterHits.begin())->Channel()),
	  clusterID, (*clusterHits.begin())->WireID().planeID()); */
      
      //   ++clusterID;
      //  ccol.push_back(cluster);
      clusHitsOut.push_back(clusterHits);
      slopevec.push_back(slope);
      totalQvec.push_back(totalQ);
      //Turn off hit sharing. T. Yang 9/14/12
      //	      //allow double assignment of first and last hits
      //	      for(size_t i = 0; i < lastHits.size(); ++i){ 
      //		if(skip[hitTemp[lastHits[i]]] ==1){
      //		  channel = hit[hitTemp[lastHits[i]]]->Wire()->RawDigit()->Channel();  
      //		  if( channel == sc || channel == ec) skip[i] = 0;
      //		}
      //	      }
      
    }// end if !std::isnan
    
    nLinesFound++;
    
    if(nLinesFound>(unsigned int)fMaxLines)
      break;
    
    
  }// end loop over hits
  
  // saves a bitmap image of the accumulator (useful for debugging), 
  // with scaling based on the maximum cell value
  if(fSaveAccumulator){   
    unsigned char *outPix = new unsigned char [accDx*accDy];
    //finds the maximum cell in the accumulator for image scaling
    int cell, pix = 0, maxCell = 0;
    for (y = 0; y < accDy; ++y){ 
      for (x = 0; x < accDx; ++x){
	cell = c.GetCell(y,x);
	if (cell > maxCell) maxCell = cell;
      }
    }
    for (y = 0; y < accDy; ++y){
      for (x = 0; x < accDx; ++x){ 
	//scales the pixel weights based on the maximum cell value     
	if(maxCell > 0)
	  pix = (int)((1500*c.GetCell(y,x))/maxCell);
	outPix[y*accDx + x] = pix;
      }
    }
    
    HLSSaveBMPFile("houghaccum.bmp", outPix, accDx, accDy);
    delete [] outPix;
  }// end if saving accumulator
  
  hit.clear();
  lastHits.clear();
  // 	  if(clusterIter != clusIn.end()){
  // 	    clusterIter++;
  // 	    ++cinctr;
  // 	  }
  listofxmax.clear();
  listofymax.clear();
  //}//end loop over clusters
  
  //       }//end loop over planes
  //     }// end loop over tpcs
  //   }// end loop over cryostats
  
  return clusHitsOut.size();
  
}





//------------------------------------------------------------------------------
size_t cluster::HoughBaseAlg::Transform(std::vector< art::Ptr<recob::Hit> > const& hits,
					double                                   & slope,
					double                                   & intercept)
{
  HoughTransform c;

  art::ServiceHandle<geo::Geometry> geom;
  int dx = geom->Nwires(0);               //number of wires 
  int dy = hits[0]->Wire()->NSignal();//number of time samples. 

  c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);

  for(unsigned int i=0;i < hits.size(); ++i){
    c.AddPointReturnMax(hits[i]->WireID().Wire, (int)(hits[i]->PeakTime()));
  }// end loop over hits

  //gets the actual two-dimensional size of the accumulator
  int accDx = 0;
  int accDy = 0;
  c.GetAccumSize(accDy, accDx);

  //find the weightiest cell in the accumulator.
  int xMax = 0;
  int yMax = 0;
  c.GetMax(yMax,xMax);

  //find the center of mass of the 3x3 cell system (with maxCell at the center). 
  float centerofmassx = 0.;
  float centerofmassy = 0.;
  float denom         = 0.;
    
  if(xMax > 0 && yMax > 0 && xMax+1 < accDx && yMax+1 < accDy){  
    for(int i = -1; i < 2; ++i){
      for(int j = -1; j < 2; ++j){
	denom         += c.GetCell(yMax+i,xMax+j);
	centerofmassx += j*c.GetCell(yMax+i,xMax+j);
	centerofmassy += i*c.GetCell(yMax+i,xMax+j);
      }
    }
    centerofmassx /= denom;
    centerofmassy /= denom;      
  }
  else  centerofmassx = centerofmassy = 0;

  float rho   = 0.;
  float theta = 0.;
  c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
  slope     = -1./tan(theta);    
  intercept = rho/sin(theta);
  
  ///\todo could eventually refine this method to throw out hits that are 
  ///\todo far from the hough line and refine the fit

  return hits.size();
}







