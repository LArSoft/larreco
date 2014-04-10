////////////////////////////////////////////////////////////////////////
// LSWMSBaseAlg.h
//
// LSWMSBaseAlg class
//
// Ben Carls (bcarls@fnal.gov)
//
////////////////////////////////////////////////////////////////////////
#ifndef LSWMSBASEALG_H
#define LSWMSBASEALG_H

#include "TMath.h"
#include <vector>
#include <unordered_map>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "boost/numeric/ublas/matrix.hpp" 
#include "boost/numeric/ublas/matrix_proxy.hpp" 
#include "boost/numeric/ublas/io.hpp" 

//#include "RecoBase/Hit.h"

#define PI_2 TMath::Pi()/2

enum {RET_OK, RET_ERROR};

namespace recob { 
  class Hit;
  class Cluster; 
}


    struct openCVPoint
    {
      unsigned int wire;
      unsigned int timebin;
      openCVPoint(int _wire, int _timebin)
      {
        wire=_wire;
        timebin=_timebin;
      }
      openCVPoint()
      {
        wire=0;
        timebin=0;
      }
    };


    struct DIR_POINT
    {
      openCVPoint pt;
      double vx;
      double vy;

      DIR_POINT( openCVPoint _pt, double _vx, double _vy)
      {
        pt=_pt;
        vx=_vx;
        vy=_vy;
      }
      DIR_POINT()
      {
        pt=openCVPoint(0,0);
        vx=0;
        vy=0;
      }
      void setTo14Quads();

    };



    struct protoTrackLSWMS
    {
      int clusterNumber=999999;
      int oldClusterNumber=999999;
      double clusterSlope=999999;
      double clusterIntercept=999999;
      double totalQ=-999999;
      double pMin0=999999;
      double pMin1=999999;
      double pMax0=-999999;
      double pMax1=-999999;
      double iMinWire=999999;
      double iMaxWire=-999999;
      double minWire=999999;
      double maxWire=-999999;
      double isolation=-999999;
      double showerLikeness=-999999;
      bool merged=false;
      bool showerMerged=false;
      bool mergedLeft=false;
      bool mergedRight=false;
      std::vector<art::Ptr<recob::Hit>> hits;
      protoTrackLSWMS(){
      }
      
      void Init(unsigned int num=999999, 
          double slope=999999, 
          double intercept=999999,
          double totalQTemp=-999999,
          double Min0=999999, 
          double Min1=999999, 
          double Max0=-999999, 
          double Max1=-999999,
          int    iMinWireTemp=999999,
          int    iMaxWireTemp=-999999,
          int    minWireTemp=999999,
          int    maxWireTemp=-999999,
          std::vector<art::Ptr<recob::Hit>> hitsTemp=std::vector<art::Ptr<recob::Hit>>())
      {
        clusterNumber = num;
        oldClusterNumber = num;
        clusterSlope = slope;
        clusterIntercept = intercept;
        totalQ=totalQTemp;
        pMin0 = Min0;
        pMin1 = Min1;
        pMax0 = Max0;
        pMax1 = Max1;
        iMinWire = iMinWireTemp;
        iMaxWire = iMaxWireTemp;
        minWire = minWireTemp;
        maxWire = maxWireTemp;
        merged = false;
        showerMerged = false;
        showerLikeness = 0;
        hits.swap(hitsTemp);
      }
    };


namespace cluster {
   
  class LSWMS{
  public:
    
    LSWMS();
    ~LSWMS();
     
    void Init( std::vector<art::Ptr<recob::Hit> > const& hits, int numberwires, int fTimeBins, int numbertimesamples, int fUseWMS, int fVerbose, double fGsigma, double fR);
    std::vector<int>  AddPointReturnMax(int x, int y);

    void reconfigure(fhicl::ParameterSet const& pset);


    void  SetM(int row, int col, int _m)    { m_M(row,col) = _m; }

    int   GetsampleIterator(int row)        {return m_sampleIterator[row];}

    int   GetM(int row, int col)            { return m_M(row,col); }
    double GetG(int row, int col)            { return m_G(row,col); }
    double GetGx(int row, int col)           { return m_Gx(row,col); }
    double GetGy(int row, int col)           { return m_Gy(row,col); }

    unsigned int   GetMsize1()            { return m_M.size1(); }
    unsigned int   GetMsize2()            { return m_M.size2(); }


    double GetGMean()            { return m_GMean; }

    double GetR() {return m_R;}

    int LineSegmentGeneration(const DIR_POINT& _dpOrig, protoTrackLSWMS *protoTrackLSWMSToAdd, std::vector<art::Ptr<recob::Hit> > const& hits, double& _error, int &nClustersTemp, std::vector<unsigned int>                 *fpointId_to_clusterId);
    
    // Growing and Mean-Shift
    double grow(const DIR_POINT& _dpOrig, openCVPoint& _ptDst, int _dir);
    int weightedMeanShift(const DIR_POINT& _dpOrig, DIR_POINT& _dpDst, const boost::numeric::ublas::matrix<int> & _M);
    //int weightedMeanShift(const DIR_POINT& _dpOrig, DIR_POINT& _dpDst, const cv::Mat& _M = cv::Mat() );

    void updateMask(openCVPoint _pt1, openCVPoint _pt2, protoTrackLSWMS* protoTrackLSWMSToAdd, std::vector<art::Ptr<recob::Hit> > const& hits, int &nClustersTemp,std::vector<unsigned int>                 *fpointId_to_clusterId);

    private:

    // Gradiant in the x direction
    //std::vector < std::vector < double > > m_Gx;
    boost::numeric::ublas::matrix<double> m_Gx;
    // Gradiant in the y direction
    //std::vector < std::vector < double > > m_Gy;
    boost::numeric::ublas::matrix<double> m_Gy;
    // Sum of Gradients in the x and y directions
    //std::vector < std::vector < double > > m_G;
    boost::numeric::ublas::matrix<double> m_G;
    // Points that have been visited already
    //std::vector < std::vector < int > > m_M;
    boost::numeric::ublas::matrix<double> m_M;
    //std::vector < std::vector < double > > m_hit_map; //the map of hits 
    boost::numeric::ublas::matrix<double> m_hit_map;
    boost::numeric::ublas::matrix<double> m_A;// Map of angles
    boost::numeric::ublas::matrix<int> m_hit_loc;

    // Weighted-mean shift
    bool m_useWMS;

    double m_R;
    double m_N;

	// Thresholds and variables
    double m_GMean;
    std::vector<int> m_sampleIterator;
    double m_Margin;

    // Mean-Shift central
    std::vector<openCVPoint> m_seeds;
    std::size_t m_seedsSize;

    //Verbose
    bool m_verbose;	

    double Gaussian(int x, int y, double sigma);
    double GaussianDerivativeX(int x, int y, double fGsigma);
    double GaussianDerivativeY(int x, int y, double fGsigma);

    // Constants
    static const int NOT_A_VALID_ANGLE = 5;
    static const double ANGLE_MARGIN;
    static const double MAX_ERROR;


  }; // class LSWMS





  class LSWMSBaseAlg {
    
  public:
    
    LSWMSBaseAlg(fhicl::ParameterSet const& pset); 
    virtual ~LSWMSBaseAlg();



    size_t FindLineSegments(std::vector<art::Ptr<recob::Hit> > const& hits,
                     std::vector<unsigned int>     *fpointId_to_clusterId,
                     int *nClusters,
                     std::vector<protoTrackLSWMS> *protoTrackLSWMSs);


    virtual void reconfigure(fhicl::ParameterSet const& pset);





  protected:

    void HLSSaveBMPFile(char const*, unsigned char*, int, int);

  private:

    int    fMaxLines;                      ///< Max number of lines that can be found 
    int    fMinHits;                       ///< Min number of hits in the accumulator to consider 
                                           ///< (number of hits required to be considered a line).
    int    fSaveAccumulator;               ///< Save bitmap image of accumulator for debugging?
    int    fNumAngleCells;                 ///< Number of angle cells in the accumulator 
                                           ///< (a measure of the angular resolution of the line finder). 
                                           ///< If this number is too large than the number of votes 
                                           ///< that fall into the "correct" bin will be small and consistent with noise.
    double  fMaxDistance;                   ///< Max distance that a hit can be from a line to be considered part of that line
    double  fMaxSlope;                      ///< Max slope a line can have
    int    fRhoZeroOutRange;               ///< Range in rho over which to zero out area around previously found lines in the accumulator
    int    fThetaZeroOutRange;             ///< Range in theta over which to zero out area around previously found lines in the accumulator
    double  fRhoResolutionFactor;           ///< Factor determining the resolution in rho
    int    fPerCluster;                    ///< Tells the original LSWMS algorithm to look at clusters individually, or all hits
                                           ///< at once
    int    fMissedHits;                    ///< Number of wires that are allowed to be missed before a line is broken up into
                                           ///< segments
    double  fMissedHitsDistance;            ///< Distance between hits in a hough line before a hit is considered missed
    double  fMissedHitsToLineSize;          ///< Ratio of missed hits to line size for a line to be considered a fake

    int    fUseWMS;
    int    fVerbose;
    int    fR;


    int          fTimeBins;
    int          fMaxCorners;
    double       fGsigma;
    int          fWindow;
    double       fThreshold;
    int          fSaveVertexMap;



    void VSSaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy);


  protected:

    friend class LSWMSTransformClus;
  };
  
  
}// namespace

#endif // LSWMSBASEALG_H
