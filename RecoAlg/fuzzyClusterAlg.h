/////////////////////////////////////////////////////////////////
/*===================================================================
  The standard implementation of FLAME data clustering algorithm.

  FLAME (Fuzzy clustering by Local Approximation of MEmberships)
  was first described in:
  "FLAME, a novel fuzzy clustering method for the analysis of DNA
  microarray data", BMC Bioinformatics, 2007, 8:3.
  Available from: http://www.biomedcentral.com/1471-2105/8/3
  
  Copyright(C) 2007, Fu Limin (phoolimin@gmail.com).
  All rights reserved.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
  2. The origin of this software must not be misrepresented; you must 
     not claim that you wrote the original software. If you use this 
     software in a product, an acknowledgment in the product 
     documentation would be appreciated but is not required.
  3. Altered source versions must be plainly marked as such, and must
     not be misrepresented as being the original software.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
===================================================================*/
////////////////////////////////////////////////////////////////////
#ifndef fuzzyClusterALG_H
#define fuzzyClusterALG_H
#include <vector>
#include <cmath>
#include <iostream>
#include <stdint.h>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

#include "RecoAlg/HoughBaseAlg.h"
#include "RecoAlg/DBScanAlg.h"
#include "Geometry/Geometry.h"

#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVectorF.h"
#include "TVector.h"
#include "TH1.h"


/* Since data for clustering are usually noisy,
 * so it is not very necessary to have EPSILON extremely small.
 */
#define EPSILON 1E-9

typedef double (*DistFunction)( std::vector<double> *x, std::vector<double> *y, int m );
    
extern DistFunction basicDistFunctions[];

class TH1F;





    
    //struct IntArray
    //{
            ////int *array;
      //std::vector<int> array;
      //int  size;
      //int  bufsize;
    //};
    
    /* For sorting and storing the orignal indices. */
    struct Indexdouble
    {
    	int   index;
    	double value;
    };
    

    typedef struct Indexdouble Indexdouble;
    //typedef struct IntArray IntArray;


        enum DistSimTypes
        {
        	DST_USER = 0,
        	DST_EUCLID ,
        	DST_COSINE ,
        	DST_PEARSON ,
        	DST_UC_PEARSON ,
        	DST_SQ_PEARSON ,
        	DST_DOT_PROD ,
        	DST_COVARIANCE ,
        	DST_MANHATTAN ,
        	DST_NULL
        };



    enum FlameObjectTypes
    {
    	OBT_NORMAL ,
    	OBT_SUPPORT ,
    	OBT_OUTLIER
    };





namespace recob { class Hit; }

namespace cluster{


  // This stores information about a showerlike cluster
  class showerCluster
    {
      public:
        int clusterNumber=-999999;
        std::vector<protoTrack> clusterProtoTracks;
        showerCluster (protoTrack protoTrackTemp)
        {
          clusterNumber=protoTrackTemp.clusterNumber;
          clusterProtoTracks.push_back(protoTrackTemp);
        }

        void addProtoTracks(std::vector<protoTrack> tracksToAdd){
          
          for(auto tracksToAddItr = tracksToAdd.begin(); tracksToAddItr != tracksToAdd.end(); tracksToAddItr++)
            tracksToAddItr->clusterNumber = clusterNumber;
          clusterProtoTracks.insert(clusterProtoTracks.end(),tracksToAdd.begin(),tracksToAdd.end());
        }
        
        void clearProtoTracks(){
          clusterProtoTracks.clear();
        }

    };

  // This stores information about a tracklike cluster
  class trackCluster
    {
      public:
        int clusterNumber=-999999;
        std::vector<protoTrack> clusterProtoTracks;
        trackCluster (protoTrack protoTrackTemp)
        {
          clusterNumber=protoTrackTemp.clusterNumber;
          clusterProtoTracks.push_back(protoTrackTemp);
        }

        void addProtoTracks(std::vector<protoTrack> tracksToAdd){

          for(auto tracksToAddItr = tracksToAdd.begin(); tracksToAddItr != tracksToAdd.end(); tracksToAddItr++)
            tracksToAddItr->clusterNumber = clusterNumber;
          clusterProtoTracks.insert(clusterProtoTracks.end(),tracksToAdd.begin(),tracksToAdd.end());
        }
        
        void clearProtoTracks(){
          clusterProtoTracks.clear();
        }

    };




    class Flame {
      public:



        /* Sort until the smallest "part" items are sorted. */
        void PartialQuickSort( std::vector<Indexdouble> *data, int first, int last, int part );


        /* Free allocated memory, and set all fields to zero. */
        //void Flame_Clear( Flame *self );
        
        /* Set a NxM data matrix, and compute distances of type T.
         * 
         * If T==DST_USER or T>=DST_NULL, and Flame::distfunc member is set,
         * then Flame::distfunc is used to compute the distances;
         * Otherwise, Flame_Euclidean() is used. */
        void Flame_SetDataMatrix( std::vector<std::vector<double>> *data, int N, int M, int T );
        
        /* Set a pre-computed NxN distance matrix. */
        void Flame_SetDistMatrix( std::vector<std::vector<double>> *data, int N );
        
        /* Define knn-nearest neighbors for each object 
         * and the Cluster Supporting Objects (CSO). 
         * 
         * The actual number of nearest neighbors could be large than knn,
         * if an object has neighbors of the same distances.
         *
         * Based on the distances of the neighbors, a density can be computed
         * for each object. Objects with local maximum density are defined as
         * CSOs. The initial outliers are defined as objects with local minimum
         * density which is less than mean( density ) + thd * stdev( density );
         */
        void Flame_DefineSupports( int knn, double thd );
        
        /* Local Approximation of fuzzy memberships.
         * Stopped after the maximum steps of iterations;
         * Or stopped when the overall membership difference between
         * two iterations become less than epsilon. */
        void Flame_LocalApproximation( int steps, double epsilon );
        
        /* Construct clusters.
         * If 0<thd<1:
         *   each object is assigned to all clusters in which
         *   it has membership higher than thd; if it can not be assigned
         *   to any clusters, it is then assigned to the outlier group.
         * Else:
         *   each object is assigned to the group (clusters/outlier group)
         *   in which it has the highest membership. */
        void Flame_MakeClusters( double thd );

        void Flame_SetMatrix( std::vector<std::vector<double>> *data, int n, int m );


    	int simtype;
    
    	/* Number of objects */
    	int N;
    
    	/* Number of K-Nearest Neighbors */
    	int K;
    
    	/* Upper bound for K defined as: sqrt(N)+10 */
    	int KMAX;
    
    	/* Stores the KMAX nearest neighbors instead of K nearest neighbors
    	 * for each objects, so that when K is changed, weights and CSOs can be
    	 * re-computed without referring to the original data.
    	 */
            //int   **graph;
        std::vector<std::vector<int>> graph;
        /* Distances to the KMAX nearest neighbors. */
            //double **dists;
        std::vector<std::vector<double>> dists;
    
    	/* Nearest neighbor count.
    	 * it can be different from K if an object has nearest neighbors with
    	 * equal distance. */
            //int    *nncounts;
        std::vector<double> nncounts;
            //double **weights;
        std::vector<std::vector<double>> weights;
    
	/* Number of identified Cluster Supporting Objects */
	int cso_count;
	//char *obtypes;
        std::vector<char> obtypes;

	//double **fuzzyships;
        std::vector<std::vector<double>> fuzzyships;
	
	/* Number of clusters including the outlier group */
	int count;
	/* The last one is the outlier group. */
        //std::vector<IntArray> clusters;
        std::vector<std::vector<int>> clusters;
	
	DistFunction distfunc;


    



      private:



    };/// End Flame class





  //--------------------------------------------------------------- 
  class fuzzyClusterAlg {
  public:
    
    
    fuzzyClusterAlg(fhicl::ParameterSet const& pset);
    virtual ~fuzzyClusterAlg();
    
    void reconfigure(fhicl::ParameterSet const& p);
    void InitFuzzy(std::vector<art::Ptr<recob::Hit> >& allhits, std::set<uint32_t> badChannels);
    // Three differnt version of the clustering code
    void run_fuzzy_cluster(std::vector<art::Ptr<recob::Hit> >& allhits);     
   

    std::vector<std::vector<unsigned int> > fclusters;               ///< collection of something
    std::vector<std::vector<double> >       fps;                     ///< the collection of points we are working on     
    std::vector<unsigned int>               fpointId_to_clusterId;   ///< mapping point_id -> clusterId     
    std::vector<std::vector<double> >       fsim;                    ///<
    std::vector<std::vector<double> >       fsim2;            	     ///<
    std::vector<std::vector<double> >       fsim3;            	     ///<
    double fMaxWidth;

    //Needed for Ben's FLAME cluster
    int NNumOfRows;
    int MNumOfCols;

    // Get functions and structures from HoughBaseAlg
    //friend class HoughBaseAlg;

   

    static double Flame_Euclidean( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_Cosine( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_Pearson( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_UCPearson( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_SQPearson( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_DotProduct( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_Covariance( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_Manhattan( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_CosineDist( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_PearsonDist( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_UCPearsonDist( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_SQPearsonDist( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_DotProductDist( std::vector<double> *x, std::vector<double> *y, int m );
    static double Flame_CovarianceDist( std::vector<double> *x, std::vector<double> *y, int m );



   //double **data = NULL;
   std::vector<std::vector<double>> data;




  private:
   

    // The distance metric chosen for clustering, you likely do not need to modify this
    int fDistanceMetric;
    // The number of hits to be considered per k-nearest neighbors cluster 
    int fKNN;
    // The maximum number of iterations to try, needed for FLAME clustering
    int fIterations;
    // The limit in the difference between memberships when FLAME clustering stops
    double fEpsilon;
    // Sets the threshold parameter in FLAME cluster, it effectively sets a lower limit on hit density for whether a hit is considered part of a cluster or an outlier
    double fThreshold;


    // Run the Hough line finder?
    bool    fRunHough;


    bool    fDoFuzzyRemnantMerge;           ///< Tell the algorithm to merge fuzzy cluster remnants into showers or tracks (0-off, 1-on)
    double  fFuzzyRemnantMergeCutoff;       ///< cut off on merging the fuzzy cluster remnants into the nearest shower or track 

    bool    fDoTrackClusterMerge;           ///< Turn on cut on product of charge asymmetry and sin of angle between slopes of lines
    double  fTrackClusterMergeCutoff;          ///< Max distance between Hough lines before two lines are merged (muon tracks), 
    double  fChargeAsymAngleCut;            ///< Cut on product of charge asymmetry and sin of angle between slopes of lines
    double  fSigmaChargeAsymAngleCut;       ///< Cut on product of charge asymmetry and sin of angle between slopes of lines
  
    bool    fDoShowerClusterMerge;          ///< Turns on shower Hough line merging (0-off, 1-on)
    double  fShowerClusterMergeAngle;       ///< Max angle between slopes before two lines are merged, for lines in shower line regions
    double  fShowerClusterMergeCutoff;    ///< Max distance between Hough lines before two lines are merged (electron showers),

    bool    fDoShowerTrackClusterMerge;     ///< Turn on cut on product of charge asymmetry and sin of angle between slopes of lines
    double  fShowerTrackClusterMergeCutoff;    ///< Max distance between Hough lines before two lines are merged (electron showers),
    double  fShowerTrackClusterMergeAngle;  ///< Max angle between slopes before two lines are merged, for lines in shower line regions
    
    double  fShowerLikenessCut;             ///< Cut on shower likeness (larger the more shower like, smaller the less shower like)

    int    fMaxVertexLines;                ///< Max number of line end points allowed in a Hough line merge region for a merge to happen
    double  fVertexLinesCutoff;             ///< Size of the vertex region to count up lines for fMaxVertexLines 







    void mergeHoughLinesBySegment(unsigned int k,
        std::vector<protoTrack> *protoTracks, 
        double xyScale,
        int mergeStyle,
        double wire_dist,
        double tickToDist);

    bool mergeTrackClusters(unsigned int k,
        std::vector<trackCluster> *trackClusters, 
        double xyScale,
        double wire_dist,
        double tickToDist);

    bool mergeShowerClusters(unsigned int k,
        std::vector<showerCluster> *showerClusters, 
        double xyScale,
        double wire_dist,
        double tickToDist);
    
    bool mergeShowerTrackClusters(showerCluster *showerClusterI, 
        trackCluster *trackClusterJ, 
        double xyScale,
        double wire_dist,
        double tickToDist);

    //std::vector<lineSlope> linesFound;
    double HoughLineDistance(double p0MinLine1, 
        double p1MinLine1, 
        double p0MaxLine1, 
        double p1MaxLine1, 
        double p0MinLine2, 
        double p1MinLine2, 
        double p0MaxLine2, 
        double p1MaxLine2);
    bool   HoughLineIntersect(double x11,
        double  y11,
        double  x12,
        double  y12,
        double  x21,
        double  y21,
        double  x22,
        double  y22);
    double PointSegmentDistance(double px,
        double  py,
        double  x1,
        double  y1,
        double  x2,
        double  y2);

    double DistanceBetweenHits(
        art::Ptr<recob::Hit> hit0,
        art::Ptr<recob::Hit> hit1,
        double wire_dist,
        double tickToDist);

    // noise vector
    std::vector<bool>      fnoise;	
    std::vector<bool>      fvisited;					     
    std::vector<double>    fWirePitch;     ///< the pitch of the wires in each plane
    std::set<uint32_t>     fBadChannels;   ///< set of bad channels in this detector
    std::vector<uint32_t>  fBadWireSum;    ///< running total of bad channels. Used for fast intervening 
                                           ///< dead wire counting ala fBadChannelSum[m]-fBadChannelSum[n]. 

    


    //void IntArray_Push( IntArray *self, int value );


    // Object used for Hough transforms
    HoughBaseAlg fHBAlg;        

    // Object used for DBScan
    DBScanAlg fDBScan;        

    art::ServiceHandle<geo::Geometry> fGeom; ///< handle to geometry service

  }; // class fuzzyClusterAlg
    




} // namespace

#endif // ifndef fuzzyClusterALG_H














