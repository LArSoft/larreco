/////////////////////////////////////////////////////////////////////
///
/// LSWMSBaseAlg class
///
/// Ben Carls, bcarls@fnal.gov
///
/// This is an implementation of LSWMS line finding by Marcos Nieto <marcos.nieto.doncel@gmail.com>.
/// Details can be found on his homepage at: http://sourceforge.net/projects/lswms/
//
///
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
#include "RecoAlg/LSWMSBaseAlg.h"
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

const double cluster::LSWMS::ANGLE_MARGIN = 22.5f;
const double cluster::LSWMS::MAX_ERROR = 0.19625f; // ((22.5/2)*CV_PI/180

void DIR_POINT::setTo14Quads()
{

	if(vx < 0)
	{
		vx = -vx;
		vy = -vy;
	}	
}


namespace cluster{
  const unsigned int kNO_CLUSTER    = UINT_MAX;
  const unsigned int kNOISE_CLUSTER = UINT_MAX-1;
}

//------------------------------------------------------------------------------
cluster::LSWMSBaseAlg::LSWMSBaseAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

//-----------------------------------------------------------------------------
double cluster::LSWMS::Gaussian(int x, int y, double sigma)
{
  double Norm=1./std::sqrt(2*TMath::Pi()*pow(sigma,2));
  double value=Norm*exp(-(pow(x,2)+pow(y,2))/(2*pow(sigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
double cluster::LSWMS::GaussianDerivativeX(int x,int y,double fGsigma)
{
  double Norm=1./(std::sqrt(2*TMath::Pi())*pow(fGsigma,3));
  double value=Norm*(-x)*exp(-(pow(x,2)+pow(y,2))/(2*pow(fGsigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
double cluster::LSWMS::GaussianDerivativeY(int x,int y,double fGsigma)
{
  double Norm=1./(std::sqrt(2*TMath::Pi())*pow(fGsigma,3));
  double value=Norm*(-y)*exp(-(pow(x,2)+pow(y,2))/(2*pow(fGsigma,2)));
  return value;
}

//------------------------------------------------------------------------------
cluster::LSWMSBaseAlg::~LSWMSBaseAlg()
{
}

//------------------------------------------------------------------------------
void cluster::LSWMSBaseAlg::reconfigure(fhicl::ParameterSet const& pset)
{ 
  fTimeBins      = pset.get< int    >("TimeBins");
  fUseWMS        = pset.get< int    >("UseWMS");
  fGsigma        = pset.get< double >("Gsigma");
  fR             = pset.get< double >("R");
  fVerbose       = pset.get< int    >("Verbose");
  return;
}

//------------------------------------------------------------------------------
cluster::LSWMS::LSWMS()
{  
  //m_accum=NULL;
}



//------------------------------------------------------------------------------
void cluster::LSWMS::Init( std::vector<art::Ptr<recob::Hit> > const& hits, 
    int numberwires, 
    int fTimeBins, 
    int numbertimesamples, 
    int fUseWMS, 
    int fVerbose, 
    double fGsigma, 
    double fR)
{




  m_useWMS = fUseWMS;
  m_verbose = fVerbose;
  m_R = fR;

  // Seeds
  m_seedsSize = 0;

  m_N = 2*m_R + 1;

  // Gradiant in the x direction
  //m_Gx.resize(numberwires);
  m_Gx.resize(fTimeBins, numberwires);
  // Gradiant in the y direction
  //m_Gy.resize(numberwires);
  m_Gy.resize(fTimeBins, numberwires);
  // Sum of Gradients in the x and y directions
  //m_G.resize(numberwires);
  m_G.resize(fTimeBins, numberwires);
  // Points that have been visited already
  //m_M.resize(numberwires);
  m_M.resize(fTimeBins,numberwires);
  m_A.resize(fTimeBins,numberwires);
  //m_hit_map(numberwires); //the map of hits 
  m_hit_map.resize(fTimeBins,numberwires); //the map of hits 
  m_hit_loc.resize(fTimeBins,numberwires); //the location of hits 

  //m_sampleIterator = std::vector<int>(fTimeBins*numberwires, 0);
  //for(unsigned int k=0; k<m_sampleIterator.size(); k++)
    //m_sampleIterator[k] = k;

  double wx[49] = {0.};
  double wy[49] = {0.};
  int ctr = 0;
  for(int i = -3; i < 4; ++i){
    for(int j = 3; j > -4; --j){
      wx[ctr] = GaussianDerivativeX(i,j,fGsigma);
      wy[ctr] = -GaussianDerivativeY(i,j,fGsigma);
      ++ctr;
    }
  }

  //zero the matrices
  m_hit_map = boost::numeric::ublas::scalar_matrix<double>(m_hit_map.size1(), m_hit_map.size2(), 0);
  m_hit_loc = boost::numeric::ublas::scalar_matrix<double>(m_hit_loc.size1(), m_hit_loc.size2(), 0);
  m_M = boost::numeric::ublas::scalar_matrix<double>(m_M.size1(), m_M.size2(), 0);
  m_G = boost::numeric::ublas::scalar_matrix<double>(m_G.size1(), m_G.size2(), 0);
  m_Gx = boost::numeric::ublas::scalar_matrix<double>(m_Gx.size1(), m_Gx.size2(), 0);
  m_Gy = boost::numeric::ublas::scalar_matrix<double>(m_Gy.size1(), m_Gy.size2(), 0);

  // Set m_A to NOT_A_VALID_ANGLE 
  m_A = boost::numeric::ublas::scalar_matrix<double>(m_A.size1(), m_A.size2(), NOT_A_VALID_ANGLE);	

  m_hit_loc = boost::numeric::ublas::scalar_matrix<int>(m_hit_loc.size1(), m_hit_loc.size2(), -1);	

  // Angular m_margin
  m_Margin = (double)(ANGLE_MARGIN*PI/180);

  unsigned int wire    = 0;
  double hit_AverageTime = 0;
  for(auto hitsItr = hits.begin(); hitsItr != hits.end(); ++hitsItr++){
    wire = (*hitsItr)->WireID().Wire;
    hit_AverageTime = ((*hitsItr)->StartTime()+(*hitsItr)->EndTime())/2.;
    for(int timebin = 1; timebin < fTimeBins-1; ++timebin){
      if( hit_AverageTime > (double)((double)timebin-0.5)*((double)numbertimesamples/(double)fTimeBins) 
	&& hit_AverageTime < (double)((double)timebin+0.5)*((double)numbertimesamples/(double)fTimeBins)){ 	       
	m_hit_loc(timebin,wire) = hitsItr-hits.begin();
	break;
      }
    }
    
    //pixelization using a Gaussian
    for(int j = 0; j <= (int)((*hitsItr)->EndTime()-(*hitsItr)->StartTime()+.5); ++j) {
      if((int)(((*hitsItr)->StartTime()+j)*((double)fTimeBins/(double)numbertimesamples)+.5) >= 0 &&
	(int)(((*hitsItr)->StartTime()+j)*((double)fTimeBins/(double)numbertimesamples)+.5) < fTimeBins)
	m_hit_map((int)(((*hitsItr)->StartTime()+j)*((double)fTimeBins/(double)numbertimesamples)+.5),wire) += Gaussian((int)(j-(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+.5),0,(*hitsItr)->EndTime()-(*hitsItr)->StartTime());
    }
  }









  // Gaussian derivative convolution 
  // Ben: I think this is the gradient
  double GxTotal = 0;
  double GyTotal = 0;
  int countG = 0;

  int windex = 0;//the wire index to make sure the end point finder does not fall off the edge of the hit map
  int tindex = 0;//the time index to make sure the end point finder does not fall off the edge of the hit map
  int n      = 0; //index of window cell. There are 49 cells in the 7X7 Gaussian and Gaussian derivative windows
  double m_hit_map_value = 0;
  double m_Gx_value = 0;
  double m_Gy_value = 0;
  double m_G_value  = 0;
  double wx_value=0;
  double wy_value=0;
  for(wire = 1; wire < (unsigned int)numberwires-1; ++wire){
    
    for(int timebin = 1; timebin < fTimeBins-1; ++timebin){

      m_Gx_value = m_Gx(timebin,wire);
      m_Gy_value = m_Gy(timebin,wire);
      m_G_value  = m_G(timebin,wire);
      n = 0;
      ++countG;

      //This if statement may change its behavior toward being less efficient
      if(m_hit_map(timebin,wire) > 0){
        for(int i = -3; i <= 3; ++i) {
          windex = wire+i;
          if(windex < 0 ) windex = 0;
          // this is ok, because the line before makes sure it's not negative
          else if (windex >= numberwires) windex = numberwires-1; 
          
          for(int j = -3; j <= 3; ++j){
            tindex = timebin+j;
            if(tindex < 0) tindex=0;
            else if(tindex >= fTimeBins) tindex = fTimeBins-1;
  
	    m_hit_map_value = m_hit_map(tindex,windex);
	    wx_value = wx[n];
	    wy_value = wy[n];


	    m_Gx_value += wx_value*m_hit_map_value;  
	    GxTotal += std::abs(wx_value*m_hit_map_value);
	    m_Gy_value += wy_value*m_hit_map_value; 
	    GyTotal += std::abs(wy_value*m_hit_map_value);
	    m_G_value += std::abs(wx_value*m_hit_map_value)+std::abs(wy_value*m_hit_map_value);  
	    ++n;
          } // end loop over j
        } // end loop over i
      }
      m_Gx(timebin,wire) += m_Gx_value;
      m_Gy(timebin,wire) += m_Gy_value;
      m_G(timebin,wire) += m_G_value;
    } // end loop over time bins
  }// end loop over wires


  // Ben: Find average gradient, may need to multiply this by 3
  if( !m_useWMS )
  	m_GMean = 3*(std::abs(GxTotal)+std::abs(GyTotal))/((double) countG);
  else
  	m_GMean = (std::abs(GxTotal)+std::abs(GyTotal))/((double) countG);

}



//------------------------------------------------------------------------------
size_t cluster::LSWMSBaseAlg::FindLineSegments(std::vector<art::Ptr<recob::Hit> > const& hits,
					std::vector<unsigned int>                 *fpointId_to_clusterId,
					int                                       *nClusters,
                                        std::vector<protoTrack>                   *linesFound
					)
{

  int nClustersTemp = *nClusters;
  
  LSWMS l;


  // Define the prototrack object
  protoTrack protoTrackToLoad;



  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<util::DetectorProperties> detp;

  //int windex = 0;//the wire index to make sure the end point finder does not fall off the edge of the hit map
  //int tindex = 0;//the time index to make sure the end point finder does not fall off the edge of the hit map
  //int n      = 0; //index of window cell. There are 49 cells in the 7X7 Gaussian and Gaussian derivative windows
  unsigned int numberwires = 0;
  double numbertimesamples = 0.;
  std::vector<double> Cornerness2;
  //gaussian window definitions. The cell weights are calculated here to help the algorithm's speed
      
  geo::WireID wid = hits[0]->WireID();
  numberwires = geom->Cryostat(wid.Cryostat).TPC(wid.TPC).Plane(wid.Plane).Nwires();
  numbertimesamples = hits[0]->Wire()->NSignal();
  mf::LogInfo("LSWMS") << " --- endpoints check " 
      		       << numberwires << " " 
      		       << numbertimesamples << " " 
      		       << fTimeBins;


  l.Init(hits,numberwires,fTimeBins,numbertimesamples,fUseWMS,fVerbose,fGsigma,fR);




  //Ben: find the first hit that's above the mean gradient, GMean
  double error = 0;
  unsigned int wire    = 0;
  unsigned int timebin = 0; 
  for(auto hitsItr = hits.begin(); hitsItr != hits.end(); ++hitsItr++){
    wire = (*hitsItr)->WireID().Wire;
    double hit_AverageTime = ((*hitsItr)->StartTime()+(*hitsItr)->EndTime())/2.;
    for(int itimebin = 1; itimebin < fTimeBins-1; ++itimebin){
      if( hit_AverageTime > (double)((double)itimebin-0.5)*((double)numbertimesamples/(double)fTimeBins) 
	&& hit_AverageTime < (double)((double)itimebin+0.5)*((double)numbertimesamples/(double)fTimeBins)){ 	       
	timebin=itimebin;
	break;
      }
    }
    if(l.GetM(timebin,wire)>0)
      continue;
    if(l.GetG(timebin,wire) > l.GetGMean()){
      DIR_POINT dpOrig(openCVPoint(wire,timebin),l.GetGx(timebin,wire),l.GetGy(timebin,wire));
      mf::LogVerbatim("LSWMSBaseAlg") << "-------------------------------";
     
      //if(fVerbose) { printf("Try dpOrig=(%d,%d,%.2f,%.2f)...\n", dpOrig.pt.wire, dpOrig.pt.timebin, dpOrig.vx, dpOrig.vy); }
      mf::LogVerbatim("LSWMSBaseAlg") << "Try dpOrig=("<<dpOrig.pt.wire<<","<<dpOrig.pt.timebin<<","<<dpOrig.vx<<","<<dpOrig.vy<<")..."; 
    
      l.LineSegmentGeneration(dpOrig, &protoTrackToLoad, hits, error, nClustersTemp, fpointId_to_clusterId);
      linesFound->push_back(protoTrackToLoad);
    }
    l.SetM(timebin,wire,255);
  }



  //double error = 0;
  //// Loop over the image
  //int x0, y0;
  //int kIterator = 0;
  //double imgSize=numberwires*fTimeBins;
  //while(true)
  //{				
    //if (kIterator == imgSize)
      //break;
    //x0 = l.GetsampleIterator(kIterator)%numberwires;
    //y0 = l.GetsampleIterator(kIterator)/numberwires;
    //kIterator++;
    //if(l.GetM(y0,x0)>0)
      //continue;
    //l.SetM(y0,x0,255);
    ////if(l.GetG(timebin,wire) > 0)
    //if(l.GetG(y0,x0) > l.GetGMean()){
      //DIR_POINT dpOrig(openCVPoint(x0,y0),l.GetGx(y0,x0),l.GetGy(y0,x0));
      //if(fVerbose) { printf("-------------------------------\n"); }
      //if(fVerbose) { printf("Try dpOrig=(%d,%d,%.2f,%.2f)...\n", dpOrig.pt.wire, dpOrig.pt.timebin, dpOrig.vx, dpOrig.vy); }
      //l.LineSegmentGeneration(dpOrig, &protoTrackToLoad, hits, error, nClustersTemp, fpointId_to_clusterId);
      //linesFound->push_back(protoTrackToLoad);
    //}
	//for(unsigned int j=y0-l.GetR(); j<=y0+l.GetR(); ++j)
	//{
		//for(unsigned int i=x0-l.GetR(); i<=x0+l.GetR(); ++i){		
			//if(j >= l.GetMsize1() || i >= l.GetMsize2() || j < 0 || i < 0)
			  //continue;
			//l.SetM(j,i,255);				
		//}
	//} 
  //}//End while




  ////Ben: find the first hit that's above the mean gradient, GMean
  //double error = 0;
  //for(unsigned int wire = 1; wire < numberwires-1; ++wire){
    //for(int timebin = 1; timebin < fTimeBins-1; ++timebin){
      //if(l.GetM(timebin,wire)>0)
        //continue;
      //l.SetM(timebin,wire,255);

      //if(l.GetG(timebin,wire) > l.GetGMean()){
        //DIR_POINT dpOrig(openCVPoint(wire,timebin),l.GetGx(timebin,wire),l.GetGy(timebin,wire));
	//if(fVerbose) { printf("-------------------------------\n"); }
	//if(fVerbose) { printf("Try dpOrig=(%d,%d,%.2f,%.2f)...\n", dpOrig.pt.wire, dpOrig.pt.timebin, dpOrig.vx, dpOrig.vy); }
        //l.LineSegmentGeneration(dpOrig, &protoTrackToLoad, hits, error, nClustersTemp, fpointId_to_clusterId);
        //linesFound->push_back(protoTrackToLoad);
      //}
    //}
  //}




  *nClusters=nClustersTemp;

  return 1; 

  
  
  
  
}




//------------------------------------------------------------------------------
int cluster::LSWMS::LineSegmentGeneration(const DIR_POINT& _dpOrig, protoTrack* protoTrackToAdd, std::vector<art::Ptr<recob::Hit> > const& hits, double& _error, int &nClustersTemp,std::vector<unsigned int>                 *fpointId_to_clusterId)
{
	// **********************************************
	// Starts at dpOrig and generates lSeg
	// 
	// Args:
	// 	-> dpOrig - starting DIR_POINT 
	//	<- lSeg - detected line segment
	// Ret:
	// 	RET_OK - lSeg created
	// 	RET_ERROR - lSeg not created
	// **********************************************


	// Check input data
	if(_dpOrig.pt.wire < 0 || _dpOrig.pt.wire >= m_G.size2() || _dpOrig.pt.timebin<0 || _dpOrig.pt.timebin >= m_G.size1())
		return RET_ERROR;

	// Find best candidate with Mean-Shift
	// -----------------------------------------------------
	DIR_POINT dpCentr = _dpOrig;	
	if( m_useWMS )
	{
		
		{

			 mf::LogVerbatim("LSWMSBaseAlg") << " Mean-Shift(Centr): from("
			   <<_dpOrig.pt.wire<<","<<_dpOrig.pt.timebin<<","<<_dpOrig.vx<<","<<_dpOrig.vy<<") to..."; 

		}
		int retMSC = weightedMeanShift(_dpOrig, dpCentr, m_M); /// Since it is receiveing m_M, it knows if it has been visited or no
		mf::LogVerbatim("LSWMSBaseAlg") << " Mean-Shift(Centr): from("
		   <<dpCentr.pt.wire<<","<<dpCentr.pt.timebin<<","<<dpCentr.vx<<","<<dpCentr.vy<<") to..."; 
		
		if(retMSC == RET_ERROR)	
		{
			// { printf("\tMean-Shift reached not a valid point\n"); }
			mf::LogVerbatim("LSWMSBaseAlg") << "  Mean-Shift reached not a valid point"; 
			return RET_ERROR;	
		}
	}
	
	// Grow in two directions from dpCentr
	// -----------------------------------------------------
	// { printf("\tGROW 1:"); fflush(stdout); }
	mf::LogVerbatim("LSWMSBaseAlg")<< "  GROW 1:"; 
	openCVPoint pt1;
	double retG1 = grow(dpCentr, pt1, 1);
	double d1 = (double)((dpCentr.pt.wire - pt1.wire)*(dpCentr.pt.wire - pt1.wire) + (dpCentr.pt.timebin - pt1.timebin)*(dpCentr.pt.timebin - pt1.timebin));
	// { printf("\tpt1(%d,%d), dist = %.2f, error=%.4f\n", pt1.wire, pt1.timebin, d1, retG1); }
	mf::LogVerbatim("LSWMSBaseAlg") << "  pt1("<<pt1.wire<<","<<pt1.timebin << "), dist = "<< d1<<  ", error="<< retG1; 
	
	// { printf("\tGROW 2:"); fflush(stdout); }
	mf::LogVerbatim("LSWMSBaseAlg")<< "  GROW 2:"; 
	openCVPoint pt2;
	double retG2 = grow(dpCentr, pt2, 2);
	double d2 = (double)((dpCentr.pt.wire - pt2.wire)*(dpCentr.pt.wire - pt2.wire) + (dpCentr.pt.timebin - pt2.timebin)*(dpCentr.pt.timebin - pt2.timebin));
	// { printf("\tpt2(%d,%d), dist = %.2f, error=%.4f\n", pt2.wire, pt2.timebin, d2, retG2); }
	mf::LogVerbatim("LSWMSBaseAlg") << "  pt2("<<pt2.wire<<","<<pt2.timebin << "), dist = "<< d2<<  ", error="<< retG2; 

	if(retG1 == -1 && retG2 == -1)
		return RET_ERROR;

	// Select the most distant extremum
	if(d1<d2)
	{
		pt1 = pt2;
		_error = retG2;
		//{ printf("Longest dir is 2\n"); }
		mf::LogVerbatim("LSWMSBaseAlg") << "Longest dir is 2"; 
	}
	else
	{
		_error = retG1;
		 //{ printf("Longest dir is 1\n"); }
		 mf::LogVerbatim("LSWMSBaseAlg") << "Longest dir is 1";
	}

	// Grow to the non-selected direction, with the new orientation
	double dirX = (double)(dpCentr.pt.wire - pt1.wire);
	double dirY = (double)(dpCentr.pt.timebin - pt1.timebin);
	double norm = sqrt(dirX*dirX + dirY*dirY);
	
	if(norm>0)
	{
		dirX = dirX/norm;
		dirY = dirY/norm;
		
		DIR_POINT dpAux(dpCentr.pt, -(-dirY), -dirX);	// DIR_POINT must be filled ALWAYS with gradient vectors
		double retG = grow(dpAux, pt2, 1);
		_error = retG;
	}
	else
	{
		pt2 = dpCentr.pt;	
	}
	
	// Check
	dirX = (double)(pt1.wire -pt2.wire);
	dirY = (double)(pt1.timebin -pt2.timebin);
	if( sqrt(dirX*dirX + dirY*dirY) < m_N)
	{
		 //{ printf("Line segment not generated: Too short.\n"); }
		 mf::LogVerbatim("LSWMSBaseAlg") << "Line segment not generated: Too short.";
		return RET_ERROR;
	}
	
	// Output line segment
	 //{ printf("LSeg = (%d,%d)-(%d,%d)\n", pt2.wire, pt2.timebin, pt1.wire, pt1.timebin); }
	mf::LogVerbatim("LSWMSBaseAlg") << "LSeg = ("<<pt2.wire<<","<< pt2.timebin <<")-("<<pt1.wire <<","<< pt1.timebin<<")"; 
	//_lSeg.clear();
	////_lSeg.push_back(cv::Point(pt2.x - 2*m_R, pt2.y - 2*m_R));
	////_lSeg.push_back(cv::Point(pt1.x - 2*m_R, pt1.y - 2*m_R));
	//_lSeg.push_back(openCVPoint(pt2.wire, pt2.timebin));
	//_lSeg.push_back(openCVPoint(pt1.wire, pt1.timebin));
	
	// Update visited positions matrix
	updateMask(pt1,pt2, protoTrackToAdd, hits, nClustersTemp, fpointId_to_clusterId);

               

        return RET_OK;
}



double cluster::LSWMS::grow(const DIR_POINT& _dpOrig, openCVPoint& _ptDst, int _dir)

{


	// **********************************************
	// Finds end-point ptDst starting from dpOrig
	// 
	// Args:
	// 	-> dpOrig - starting DIR_POINT 
	//	<- ptDst - end-point
	//	-> dir - growing direction (1(+) or 2(-))
	// Ret:
	// 	error - error of line segment
	//
	// Called from lineSegmentGeneration
	// **********************************************	
	
	openCVPoint ptEnd1, ptEnd2; //auxiliar
	DIR_POINT dpEnd, dpRef;	  // auxiliar

	// Init output
	_ptDst = _dpOrig.pt;

	// Starting gradient vector and director vector
	double gX, gY;
	if(_dir == 1)
	{	
		gX = _dpOrig.vx;
		gY = _dpOrig.vy;		
	}
	else if(_dir == 2)
	{
		gX = -_dpOrig.vx;
		gY = -_dpOrig.vy;
	}
	else return RET_ERROR;

	// Compute currentAngle in 1º4º
	double error1 = 0;
	double growAngle, auxAngle, minAngle, maxAngle, diffAngle;
	//if(gX < 0)	// In this case, direction must not be fliped to 1º4º, because otherwise the sense is lost for the second grow procedure...
	//{
	//	// Move to 1º4º
	//	gX = -gX;
	//	gY = -gY;
	//}
	growAngle = atan2(gY, gX);

	// Starting point and angle - Bresenham
	openCVPoint pt1 = _dpOrig.pt;
	openCVPoint pt2(pt1.wire + (int)(1000*(-gY)), pt1.timebin + (int)(1000*(gX)));
	//cv::clipLine(m_imPadSize, pt1, pt2);	
	
	// Loop	- Bresenham
	int k1=0;
	
	int x1 = pt1.wire, x2 = pt2.wire, y1 = pt1.timebin, y2 = pt2.timebin;
	int dx = (int)std::abs((double)x2-(double)x1);
	int dy = (int)std::abs((double)y2-(double)y1);
	int sx, sy, err, e2;

	 //{ printf("From (%d,%d) to (%d,%d)...", x1, y1, x2, y2); fflush(stdout); }
	mf::LogVerbatim("LSWMSBaseAlg") << "From = ("<<x1<<","<< y1<<") to ("<<x2 <<","<< y2<<")..."; 

	if(x1 < x2) sx = 1; else sx = -1;
	if(y1 < y2) sy = 1; else sy = -1;
	err = dx-dy;
	
	int maxNumZeroPixels = 2*m_R, countZeroPixels=0;
	while(true)
	{		
		// Current value is (x1,y1)	
		// { printf("\n\tBresenham(%d,%d)", x1, y1); fflush(stdout); }
		// -------------------------------
		// Do...
		// Check if angle has been computed
		if((unsigned int) y1 >= m_A.size1() || (unsigned int) y1 < 0 || (unsigned int) x1 >= m_A.size2() || (unsigned int) x1 < 0)
		  break;

		if(m_A(y1,x1) != NOT_A_VALID_ANGLE)
			auxAngle = m_A(y1,x1);
		else
		{
			auxAngle = atan2((double)m_Gy(y1,x1), (double)m_Gx(y1,x1));
			m_A(y1,x1) = auxAngle;
		}		

		// Check early-termination of Bresenham		
		if(m_G(y1,x1) == 0) 
		{
			// printf("Zero-pixel num. %d\n", countZeroPixels);
		
			countZeroPixels++;
			if(countZeroPixels >= maxNumZeroPixels)
				break;		// No gradient point
			
		}
		
		// Check angular limits
		if(growAngle - m_Margin < -PI_2)		    // e.g. currentAngle = -80º, margin = 20º
		{
			minAngle = growAngle - m_Margin + (double)PI; // e.g. -80 -20 +180 = 80º
			maxAngle = growAngle + m_Margin;	    // e.g. -80 +20      =-60º	

			if( auxAngle < 0)
			{
				if( auxAngle > maxAngle ) break; // printf("Early-termination: auxAngle(%.2f) > maxAngle(%.2f) && auxAngle < 0\n", auxAngle, maxAngle);
				diffAngle = std::abs(growAngle - auxAngle);
			}
			else // auxAngle > 0
			{
				if( auxAngle < minAngle) break; // printf("Early-termination: auxAngle(%.2f) < minAngle(%.2f) && auxAngle > 0\n", auxAngle, minAngle);
				diffAngle = std::abs(growAngle - (auxAngle - (double)PI));
			}			
		}
		else if(growAngle + m_Margin > PI_2)		    // e.g. currentAngle = 80º, margin = 20º
		{
			minAngle = growAngle - m_Margin;	    // e.g.  80 -20      = 60º
			maxAngle = growAngle + m_Margin - (double)PI; // e.g.  80 +20 -180 = -80º

			if( auxAngle > 0 )
			{
				if( auxAngle < minAngle) break;	// printf("Early-termination: auxAngle(%.2f) < minAngle(%.2f) && auxAngle > 0\n", auxAngle, minAngle);
				diffAngle = std::abs(growAngle - auxAngle);
			}
			else // auxAngle < 0
			{
				if( auxAngle > maxAngle) break; // printf("Early-termination: auxAngle(%.2f) > maxAngle(%.2f) && auxAngle < 0\n", auxAngle, maxAngle);
				diffAngle = std::abs(growAngle - (auxAngle + (double)PI));
			}
		}
		else
		{	
			minAngle = growAngle - m_Margin;
			maxAngle = growAngle + m_Margin;
			if(auxAngle < minAngle || auxAngle > maxAngle) break; // printf("Early-termination: auxAngle(%.2f) < minAngle(%.2f) || > maxAngle(%.2f)\n", auxAngle, minAngle, maxAngle);
			
			diffAngle = std::abs(growAngle - auxAngle);
		}
		
		// If arrived here, the point is valid (inside the angular limits, and with G!=0)
		//error1 += std::abs(std::abs(m_Gx.at<short>(y1,x1)) - std::abs(gX)) +
		//	  std::abs(std::abs(m_Gy.at<short>(y1,x1)) - std::abs(gY));
		//error1 += std::abs(auxAngle - growAngle); // OJO, SI HA HABIDO DISCONTINUIDAD, ESTO NO ES CORRECTO...
		error1 += diffAngle;
		ptEnd1 = openCVPoint(x1,y1);
		k1++;

		// -------------------------------
		// Check end
		if (x1 == x2 && y1 == y2) break;

		// Update position for next iteration
		e2 = 2*err;
		if(e2 > -dy) { err = err - dy; x1 = x1 + sx;}
		if(e2 < dx)  { err = err + dx; y1 = y1 + sy;} 		
	}

	// "k1": how many points have been visited
	// "ptEnd": last valid point
	if( k1==0 ) // this means that even the closest point has not been accepted	
	{
		ptEnd1 = _dpOrig.pt;
		error1 = (double)PI;
	}	
	else error1 /= k1;
	 //{ printf(", Arrived to (%d,%d), error=%.2f", ptEnd1.wire, ptEnd1.timebin, error1); fflush(stdout); }
	 mf::LogVerbatim("LSWMSBaseAlg") << ", Arrived to (" << ptEnd1.wire<<","<<ptEnd1.timebin<<"), error="<<error1; 
	
	// Set ptDst
	_ptDst = ptEnd1;

	// Apply Mean-Shift to refine the end point	
	// printf("Check grow movement: From (%d,%d) to (%d,%d)\n", dpOrig.pt.x, dpOrig.pt.y, ptEnd1.x, ptEnd1.y);	
	 //{ printf(", Dist = (%d,%d)\n", (int)std::abs((double)ptEnd1.wire - (double)_dpOrig.pt.wire), (int)std::abs((double)ptEnd1.timebin - (double)_dpOrig.pt.timebin)); }
	 mf::LogVerbatim("LSWMSBaseAlg")<< "Dist = ("<<(int)std::abs((double)ptEnd1.wire - (double)_dpOrig.pt.wire)<<","<<(int)std::abs((double)ptEnd1.timebin - (double)_dpOrig.pt.timebin)<<")"; 
	if((int)std::abs((double)ptEnd1.wire - (double)_dpOrig.pt.wire) > m_R || (int)std::abs((double)ptEnd1.timebin - (double)_dpOrig.pt.timebin) > m_R) // ONLY IF THERE HAS BEEN (SIGNIFICANT) MOTION FROM PREVIOUS MEAN-SHIFT MAXIMA
	{		
		int counter = 0;
		while(true)
		{
			 //{ printf("\tMean-Shift(Ext): from (%d,%d,%.2f,%.2f) to...", ptEnd1.wire, ptEnd1.timebin, gX, gY); fflush(stdout); }
			mf::LogVerbatim("LSWMSBaseAlg")<<"Mean-Shift(Ext): from ("<<ptEnd1.wire<<","<<ptEnd1.timebin<<","<<gX<<","<<gY<<") to...";
			counter++;

			// Mean-Shift on the initial extremum
			// -------------------------------------------------------------
			dpEnd.pt = ptEnd1; dpEnd.vx = gX; dpEnd.vy = gY;	// gX and gY have been update in the last iter
			dpRef.pt = ptEnd1; dpRef.vx = gX; dpRef.vy = gY;

			if( m_useWMS )
			{
				int retMSExt = weightedMeanShift(dpEnd, dpRef,m_M);

				 //{ printf("(%d,%d,%.2f,%.2f)\n", dpRef.pt.wire, dpRef.pt.timebin, dpRef.vx, dpRef.vy); }
				mf::LogVerbatim("LSWMSBaseAlg")<< "("<<dpRef.pt.wire<<","<<dpRef.pt.timebin<<","<<dpRef.vx<<","<<dpRef.vy<<")";

				if(retMSExt == RET_ERROR)
				{
					// The refinement gave and incorrect value, keep last Bresenham value
					_ptDst = ptEnd1;
					return RET_OK;
				}
			}
			else
				return RET_OK;

			// Check motion caused by Mean-Shift
			if(dpRef.pt.wire == dpEnd.pt.wire && dpRef.pt.timebin == dpEnd.pt.timebin)
			{
				_ptDst = dpRef.pt;
				return RET_OK;
			}

			// Check displacement from dpOrig
			gX = (double)(dpRef.pt.timebin - _dpOrig.pt.timebin);	// 	double dX = dpRef.x - dpOrig.x; and gX = dY;
			gY = (double)(_dpOrig.pt.wire - dpRef.pt.wire);	//	double dY = dpRef.y - dpOrig.y; and gY = -dX;
			if(gX == 0 && gY == 0)
			{	
				_ptDst = dpRef.pt;
				return RET_OK;
			}
			double norm = sqrt(gX*gX + gY*gY);
			gX /= norm;
			gY /= norm;	

			// New Bresenham procedure	
			if(gX < 0)
			{
				// MOve to 1º4º
				gX = -gX;
				gY = -gY;
			}
			growAngle = atan2(gY, gX);

			int k2=0;
			double error2 = 0;

			pt2.wire = pt1.wire + (int)(1000*(-gY)); pt2.timebin = pt1.timebin + (int)(1000*(gX));

			x1 = pt1.wire; x2 = pt2.wire; y1 = pt1.timebin; y2 = pt2.timebin;
			dx = (int)std::abs((double)x2-(double)x1);	
			dy = (int)std::abs((double)y2-(double)y1);

			if(x1 < x2) sx = 1; else sx = -1;
			if(y1 < y2) sy = 1; else sy = -1,
			err = dx-dy;

			mf::LogVerbatim("LSWMSBaseAlg") << "  Refined GROW: From ("<<x1<<","<<y1<<") to ("<<x2<<","<<y2<<")...";
			while(true)
			{
				// Current value is (x1,y1)	
				// { printf("\n\tBresenham(%d,%d)", x1, y1); fflush(stdout); }
				
				// -------------------------------
				// Do...
				// Check if angle has been computed
				if(m_A(y1,x1) != NOT_A_VALID_ANGLE)
					auxAngle = m_A(y1,x1);
				else
				{
					auxAngle = atan2((double)m_Gy(y1,x1), (double)m_Gx(y1,x1));
					m_A(y1,x1) = auxAngle;
				}

				// Check early-termination of Bresenham		
				if(m_G(y1,x1) == 0) 
				{
					// printf("Zero-pixel num. %d\n", countZeroPixels);
			
					countZeroPixels++;
					if(countZeroPixels >= maxNumZeroPixels)
						break;		// No gradient point			
				}
		
				// Check angular limits
				if(growAngle - m_Margin < -PI_2)		    // e.g. currentAngle = -80º, margin = 20º
				{
					minAngle = growAngle - m_Margin + (double)PI; // e.g. -80 -20 +180 = 80º
					maxAngle = growAngle + m_Margin;	    // e.g. -80 +20      =-60º	

					if( auxAngle < 0)
					{
						if( auxAngle > maxAngle ) break; // printf("Early-termination: auxAngle(%.2f) > maxAngle(%.2f) && auxAngle < 0\n", auxAngle, maxAngle);
						diffAngle = std::abs(growAngle - auxAngle);
					}
					else // auxAngle > 0
					{
						if( auxAngle < minAngle) break; // printf("Early-termination: auxAngle(%.2f) < minAngle(%.2f) && auxAngle > 0\n", auxAngle, minAngle);
						diffAngle = std::abs(growAngle - (auxAngle - (double)PI));
					}			
				}
				else if(growAngle + m_Margin > PI_2)		    // e.g. currentAngle = 80º, margin = 20º
				{
					minAngle = growAngle - m_Margin;	    // e.g.  80 -20      = 60º
					maxAngle = growAngle + m_Margin - (double)PI; // e.g.  80 +20 -180 = -80º

					if( auxAngle > 0 )
					{
						if( auxAngle < minAngle) break;	// printf("Early-termination: auxAngle(%.2f) < minAngle(%.2f) && auxAngle > 0\n", auxAngle, minAngle);
						diffAngle = std::abs(growAngle - auxAngle);
					}
					else // auxAngle < 0
					{
						if( auxAngle > maxAngle) break; // printf("Early-termination: auxAngle(%.2f) > maxAngle(%.2f) && auxAngle < 0\n", auxAngle, maxAngle);
						diffAngle = std::abs(growAngle - (auxAngle + (double)PI));
					}
				}
				else
				{	
					minAngle = growAngle - m_Margin;
					maxAngle = growAngle + m_Margin;
					if(auxAngle < minAngle || auxAngle > maxAngle) break; // printf("Early-termination: auxAngle(%.2f) < minAngle(%.2f) || > maxAngle(%.2f)\n", auxAngle, minAngle, maxAngle);
			
					diffAngle = std::abs(growAngle - auxAngle);
				}
				
				error2 += diffAngle;
				ptEnd2 = openCVPoint(x1,y1);
				k2++;
				
				// -------------------------------

				// Check end
				if (x1 == x2 && y1 == y2) break;

				// Update position for next iteration
				e2 = 2*err;
				if(e2 > -dy) { err = err - dy; x1 = x1 + sx;}
				if(e2 < dx)  { err = err + dx; y1 = y1 + sy;} 		
			} // Bresenham while

			
			// "k2": how many points have been visited
			// "ptEnd2": last valid point
			if( k2==0 ) // this means that even the closest point has not been accepted	
			{
				ptEnd2 = _dpOrig.pt;
				error2 = (double)PI;
			}	
			else error2 = error2 / k2;			

			fflush(stdout);	// Don't really know why, but this is necessary to avoid dead loops...
			 mf::LogVerbatim("LSWMSBaseAlg") <<", Arrived to ("<<ptEnd2.wire<<","<<ptEnd2.timebin<<"), error="<<error2; 
			 mf::LogVerbatim("LSWMSBaseAlg") <<", Dist = ("<<(int)std::abs((double)ptEnd2.wire - (double)_dpOrig.pt.wire)<<","<<(int)std::abs((double)ptEnd1.timebin - (double)_dpOrig.pt.timebin)<<")";

			// Compare obtained samples
			if(error1 <= error2)
			{
				_ptDst = ptEnd1;
				return error1;
			}
			else
			{
				// Update ptEnd1 with ptEnd2 because it is better and iterate
				ptEnd1 = ptEnd2;
				k1 = k2;
				error1 = error2;
			}		
		} // Mean-Shift while
	}
	//else 
	//	printf("Not enough movement\n");

	//return RET_OK;
	return error1;
}












//------------------------------------------------------------------------------
cluster::LSWMS::~LSWMS()
{
}
int cluster::LSWMS::weightedMeanShift(const DIR_POINT& _dpOrig, DIR_POINT& _dpDst, const boost::numeric::ublas::matrix<int> & _M )
//int weightedMeanShift(const DIR_POINT& _dpOrig, DIR_POINT& _dpDst, const boost::numeric::ublas::matrix<int> & _M = boost::numeric::ublas::matrix<int>);
{
	// **********************************************
	// Refines dpOrig and creates dpDst
	// 
	// Args:
	// 	-> dpOrig - starting DIR_POINT 
	//	<- dpDst - refined DIR_POINT
	// Ret:
	// 	RET_OK - dpDst created
	// 	RET_ERROR - dpDst not found
	//
	// Called from "lineSegmentGeneration"
	// **********************************************

	// MAIN LOOP: loop until MS generates no movement (or dead-loop)
	m_seeds.clear();
	m_seedsSize = 0;
	DIR_POINT dpCurr = _dpOrig;	// The initial dp is in 1º4º
	_dpDst = _dpOrig;

	while(true)
	{

		// Check point
		if(dpCurr.pt.wire < 0 || dpCurr.pt.wire >= m_G.size2() || dpCurr.pt.timebin<0 || dpCurr.pt.timebin >= m_G.size1())
			return RET_ERROR;

		// Check direction
		if(dpCurr.vx==0 && dpCurr.vy == 0)	
			return RET_ERROR;

		// Convert to 1º4º (maybe not needed)
		//setTo14Quads(dpCurr);
		dpCurr.setTo14Quads();

		
		// Check already visited
		//if(!_M.empty())
		//{
			if(m_M(dpCurr.pt.timebin, dpCurr.pt.wire) == 255) 
			{
				return RET_ERROR;
			}
		//}

		// Check if previously used as seed for this MS-central (this is to avoid dead-loops)
		for(unsigned int i=0; i<m_seedsSize; i++)
		{
			if(m_seeds[i].wire == dpCurr.pt.wire && m_seeds[i].timebin == dpCurr.pt.timebin)
			{
				_dpDst = dpCurr;
				return RET_ERROR;
			}
		}		

		// Define bounds
		int xMin = dpCurr.pt.wire - (int)m_R;
		int yMin = dpCurr.pt.timebin - (int)m_R;
		int xMax = dpCurr.pt.wire + (int)m_R;
		int yMax = dpCurr.pt.timebin + (int)m_R;
		int offX = m_R;
		int offY = m_R;

		if( xMin < 0 || yMin < 0 || xMax+1 >= (int)m_G.size2() || yMax+1 >= (int)m_G.size1())
			return RET_ERROR;
		
		m_seeds.push_back(dpCurr.pt);
		m_seedsSize++;

		// Define rois
		//Rect roi(xMin, yMin, xMax-xMin+1, yMax-yMin+1);
		//boost::numeric::ublas::matrix<double>  gBlock = Mat(m_G, roi);
                boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<double> >  gBlock(m_G,boost::numeric::ublas::range(yMin,yMax+1), boost::numeric::ublas::range(xMin,xMax+1));
		//boost::numeric::ublas::matrix<double>  gXBlock = Mat(m_Gx, roi);
                boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<double> >  gXBlock(m_Gx,boost::numeric::ublas::range(yMin,yMax+1), boost::numeric::ublas::range(xMin,xMax+1));
		//boost::numeric::ublas::matrix<double>  gYBlock = Mat(m_Gy, roi);
                boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<double> >  gYBlock(m_Gy,boost::numeric::ublas::range(yMin,yMax+1), boost::numeric::ublas::range(xMin,xMax+1));
		//boost::numeric::ublas::matrix<double>  aBlock = Mat(m_A, roi);
                boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<double> >  aBlock(m_A,boost::numeric::ublas::range(yMin,yMax+1), boost::numeric::ublas::range(xMin,xMax+1));
		boost::numeric::ublas::matrix<double> insideBlock( gBlock.size1(), gBlock.size2()); // 0: outside, 1:inside
		
                insideBlock = boost::numeric::ublas::scalar_matrix<double>(insideBlock.size1(), insideBlock.size2(),1);	
		
                // Update angles (this is to compute angles only once)
		for(size_t j=0; j<aBlock.size1(); ++j)
		{
			for(size_t i=0; i<aBlock.size2(); ++i)
			{
				// This is guaranteed to be in 1º and 4º quadrant
				aBlock(j,i) = atan2((double)gYBlock(j,i), (double)gXBlock(j,i));
			}
		}

		// printf("dpCurr(%d,%d)(%.2f,%.2f)\n", dpCurr.pt.x, dpCurr.pt.y, dpCurr.vx, dpCurr.vy);
				

		// ----------------------------------
		// Angle analysis
		double currentAngle = atan2(dpCurr.vy, dpCurr.vx);	// output is between (-PI/2, PI/2)
		// printf("currentAngle = %.2f\n", currentAngle);
		// ----------------------------------

		//double angleShift = 0;
		int outsideCounter = 0;
		if(currentAngle - m_Margin < -PI_2)
		{
			// Shift angles according to currentAngle to avoid discontinuities			
			// printf("shift angles since %.2f - %.2f < %.2f\n", currentAngle, m_Margin, -PI_2);
			//angleShift = currentAngle;
                        aBlock -= boost::numeric::ublas::scalar_matrix<double>(aBlock.size1(), aBlock.size2(),currentAngle);	
			currentAngle = 0;
			double minAngle = currentAngle - m_Margin;
			double maxAngle = currentAngle + m_Margin;

			for(size_t j=0; j<aBlock.size1(); j++)
			{
				for(size_t i=0; i<aBlock.size2(); i++)
				{
					if(aBlock(j,i) < -PI_2) aBlock(j,i) += (double)PI;
					if(aBlock(j,i) > PI_2) aBlock(j,i) -= (double)PI;
					if(aBlock(j,i) < minAngle || aBlock(j,i) > maxAngle) 
					{
						//ptRowGBlock[i] = -1;
						insideBlock(j,i) = 0;
						outsideCounter++;
					}	
				}
			}		
			//for(size_t j=0; j<aBlock.size1(); j++)
			//{
				//double *ptRowABlock = aBlock.ptr<double>(j);
				//unsigned char *ptRowGBlock = gBlock.ptr<unsigned char>(j);
				//for(int i=0; i<aBlock.size2(); i++)
				//{
					//if(ptRowABlock[i] < -PI_2) ptRowABlock[i] += (double)PI;
					//if(ptRowABlock[i] > PI_2) ptRowABlock[i] -= (double)PI;
					//if(ptRowABlock[i] < minAngle || ptRowABlock[i] > maxAngle) 
					//{
						////ptRowGBlock[i] = -1;
						//insideBlock(j,i) = 0;
						//outsideCounter++;
					//}	
				//}
			//}		
			// Restore
                        aBlock += boost::numeric::ublas::scalar_matrix<double>(aBlock.size1(), aBlock.size2(),currentAngle);	


		}
		else if(currentAngle + m_Margin > PI_2)
		{
			// Shift angles according to currentAngle to avoid discontinuities
			// printf("shift angles since %.2f + %.2f > %.2f\n", currentAngle, m_Margin, PI_2);
			//angleShift = currentAngle;
                        aBlock -= boost::numeric::ublas::scalar_matrix<double>(aBlock.size1(), aBlock.size2(),currentAngle);	
			currentAngle = 0;

			double minAngle = currentAngle - m_Margin;
			double maxAngle = currentAngle + m_Margin;

			for(size_t j=0; j<aBlock.size1(); j++)
			{
				for(size_t i=0; i<aBlock.size2(); i++)
				{
					if(aBlock(j,i) < -PI_2) aBlock(j,i) += (double)PI;
					if(aBlock(j,i) > PI_2) aBlock(j,i) -= (double)PI;
					if(aBlock(j,i) < minAngle || aBlock(j,i) > maxAngle) 
					{
						//ptRowGBlock[i] = -1;
						insideBlock(j,i) = 0;
						outsideCounter++;
					}	
				}
			}
			//for(size_t j=0; j<aBlock.size1(); j++)
			//{
				//double *ptRowABlock = aBlock.ptr<double>(j);
				//uchar *ptRowGBlock = gBlock.ptr<uchar>(j);
				//for(size_t i=0; i<aBlock.size2(); i++)
				//{
					//if(ptRowABlock[i] < -PI_2) ptRowABlock[i] += (double)PI;
					//if(ptRowABlock[i] > PI_2) ptRowABlock[i] -= (double)PI;
					//if(ptRowABlock[i] < minAngle || ptRowABlock[i] > maxAngle) 
					//{
						////ptRowGBlock[i] = -1;
						//insideBlock(j,i) = 0;
						//outsideCounter++;
					//}	
				//}
			//}
			// Restore
                        aBlock += boost::numeric::ublas::scalar_matrix<double>(aBlock.size1(), aBlock.size2(),currentAngle);	
		}
		else
		{
			//angleShift = 0;
			double minAngle = currentAngle - m_Margin;
			double maxAngle = currentAngle + m_Margin;
			for(size_t j=0; j<aBlock.size1(); j++)
			{
				for(size_t i=0; i<aBlock.size2(); i++)
				{					
					if(aBlock(j,i) < minAngle || aBlock(j,i) > maxAngle) 
					{
						//ptRowGBlock[i] = -1;
						insideBlock(j,i) = 0;
						outsideCounter++;
					}	
				}
			}
			//for(size_t j=0; j<aBlock.size1(); j++)
			//{
				//double *ptRowABlock = aBlock.ptr<double>(j);
				//uchar *ptRowGBlock = gBlock.ptr<uchar>(j);
				//for(size_t i=0; i<aBlock.size2(); i++)
				//{					
					//if(ptRowABlock[i] < minAngle || ptRowABlock[i] > maxAngle) 
					//{
						////ptRowGBlock[i] = -1;
						//insideBlock(j,i) = 0;
						//outsideCounter++;
					//}	
				//}
			//}
		}


		// Check number of samples inside the bandwidth
		if(outsideCounter == (2*m_R+1)*(2*m_R+1))
			return RET_ERROR;

		// New (Circular) Mean angle (weighted by G)
		double sumWeight = 0;
		double foffX = 0;
		double foffY = 0;
		double meanAngle = 0;
		
		for(size_t j=0; j<gBlock.size1(); j++)
		{			
			for(size_t i=0; i<gBlock.size2(); i++)
			{
				//if(ptRowGBlock[i] != -1)
				if(insideBlock(j,i) != 0)
				{
					// This sample is inside the Mean-Shift bandwidth
					// Weighted mean of positons
					foffX += (double)(i+1)*gBlock(j,i);	  // This cannot be precomputed...
					foffY += (double)(j+1)*gBlock(j,i);
							
					// Weighted mean of angle
					meanAngle += aBlock(j,i)*gBlock(j,i);
					sumWeight += gBlock(j,i);		  
				}
			}
		}
		//for(size_t j=0; j<gBlock.size1(); j++)
		//{			
			//unsigned char *ptRowGBlock = gBlock.ptr<unsigned char>(j);
			//double *ptRowABlock = aBlock.ptr<double>(j);

			//for(size_t i=0; i<gBlock.size2(); i++)
			//{
				////if(ptRowGBlock[i] != -1)
				//if(insideBlock(j,i) != 0)
				//{
					//// This sample is inside the Mean-Shift bandwidth
					//// Weighted mean of positons
					//foffX += (double)(i+1)*ptRowGBlock[i];	  // This cannot be precomputed...
					//foffY += (double)(j+1)*ptRowGBlock[i];
							
					//// Weighted mean of angle
					//meanAngle += ptRowABlock[i]*ptRowGBlock[i];
					//sumWeight += ptRowGBlock[i];		  
				//}
			//}
		//}
		foffX /= sumWeight; foffX--;
		foffY /= sumWeight; foffY--;
		meanAngle /= sumWeight;

		// printf("meanAngle = %.2f\n", meanAngle);
				
		// Check convergence (movement with respect to the center)
		//if(cvRound(foffX) == offX && cvRound(foffY) == offY)
		if(round(foffX) == offX && round(foffY) == offY)
		{
			// Converged. Assign and return.
			_dpDst = DIR_POINT(dpCurr.pt, cos(meanAngle), sin(meanAngle));
			//setTo14Quads(_dpDst);
			_dpDst.setTo14Quads();
			return RET_OK;
		}
		else
		{
			// Not converged: update dpCurr and iterate
			dpCurr.pt.wire += round(foffX) - offX;
			dpCurr.pt.timebin += round(foffY) - offY;
			dpCurr.vx = cos(meanAngle);
			dpCurr.vy = sin(meanAngle);
		}

	}	

	return RET_OK;
}

void cluster::LSWMS::updateMask(openCVPoint _pt1, openCVPoint _pt2, protoTrack* protoTrackToAdd, std::vector<art::Ptr<recob::Hit> > const& hits, int &nClustersTemp, std::vector<unsigned int> *fpointId_to_clusterId)
{	

        art::ServiceHandle<geo::Geometry> geom;
        art::ServiceHandle<util::LArProperties> larprop;
        art::ServiceHandle<util::DetectorProperties> detprop;
        
        uint32_t     channel = hits[0]->Wire()->RawDigit()->Channel();

        double wirePitch = geom->WirePitch(geom->View(channel));
        //double xyScale  = 0.001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature())*detprop->SamplingRate()/wirePitch;
        //xyScale        *= detprop->SamplingRate()/wirePitch;
      
        double wire_dist = wirePitch;
        double tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature())*1.e-3 * detprop->SamplingRate();


        double pCornerMin[2];
        double pCornerMax[2];

        unsigned int fMaxWire = 0;
        int iMaxWire = 0;
        unsigned int fMinWire = 99999999;
        int iMinWire = -1;
        unsigned int wire;
        double totalQ = 0;

	int m_hit_loc_val = -1;

        std::vector< art::Ptr<recob::Hit> > lineHits;
        nClustersTemp++;

	// Bresenham from one extremum to the other
	int x1 = _pt1.wire, x2 = _pt2.wire, y1 = _pt1.timebin, y2 = _pt2.timebin;
	int dx = (int)std::abs((double)x2-(double)x1);
	int dy = (int)std::abs((double)y2-(double)y1);
	int sx, sy, err, e2;

	bool foundHit = false;

	if(x1 < x2) sx = 1; else sx = -1;
	if(y1 < y2) sy = 1; else sy = -1;
	err = dx-dy;
	while(true)
	{


		// Current value is (x1,y1)
		// -------------------------------
		// Do...
		// Set window to "visited=255"
		for(unsigned int j=y1-m_R; j<=y1+m_R; ++j)
		{

			for(unsigned int i=x1-m_R; i<=x1+m_R; ++i){		
				if(j >= m_M.size1() || i >= m_M.size2() || j < 0 || i < 0)
				  continue;
				m_M(j,i) = 255;				
				//unsigned char* ptRowM = m_M.ptr<uchar>(j);
				//for(int i=x1-m_R; i<=x1+m_R; ++i)		
				//ptRowM[i] = 255;		
                                 
			        if(j >= m_hit_loc.size1() || i >= m_hit_loc.size2() )
			          continue;

				m_hit_loc_val = m_hit_loc(j,i);

			        if(m_hit_loc_val<0)
			          continue;
                               
			        foundHit = true;
			        fpointId_to_clusterId->at(m_hit_loc_val) = nClustersTemp-1;
                                totalQ += hits[m_hit_loc_val]->Charge();
                                wire = hits[m_hit_loc_val]->WireID().Wire;
			       
                                
			        lineHits.push_back(hits[m_hit_loc_val]);
                                
                                if(wire < fMinWire){
                                  fMinWire = wire;
                                  iMinWire = m_hit_loc_val;
                                }
                                if(wire > fMaxWire){
                                  fMaxWire = wire;
                                  iMaxWire = m_hit_loc_val;
                                }
                       
                            
                        
                        }
                                
                            
		}
		// -------------------------------
		
		// Check end
		if (x1 == x2 && y1 == y2) break;

		// Update position for next iteration
		e2 = 2*err;
		if(e2 > -dy) { err = err - dy; x1 = x1 + sx;}
		if(e2 < dx)  { err = err + dx; y1 = y1 + sy;} 	
	}

	if(!foundHit){
	  nClustersTemp--;
	  return;
	}


        pCornerMin[0] = (hits[iMinWire]->Wire()->RawDigit()->Channel())*wire_dist;
        pCornerMin[1] = ((hits[iMinWire]->StartTime()+hits[iMinWire]->EndTime())/2.)*tickToDist;
        pCornerMax[0] = (hits[iMaxWire]->Wire()->RawDigit()->Channel())*wire_dist;
        pCornerMax[1] = ((hits[iMaxWire]->StartTime()+hits[iMaxWire]->EndTime())/2.)*tickToDist;

        double slope = (pCornerMax[1]-pCornerMin[1])/(pCornerMax[0]-pCornerMin[0]);
        double intercept = pCornerMax[1] - slope*pCornerMax[0];

        protoTrackToAdd->Init(nClustersTemp-1,
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

}
