////////////////////////////////////////////////////////////////////////
// HoughBaseAlg.h
//
// HoughBaseAlg class
//
// Ben Carls (bcarls@fnal.gov)
// Some optimization by Gianluca Petrillo (petrillo@fnal.gov)
//
// Hough transform algorithm
// ============================================================================
// 
// This implementation is a pattern recognition algorithm trying to detect
// straight lines in a image.
// In our application, the image is a hitmap on a wire plane, represented by
// wire number as one coordinate and drift time as the other.
// For each hit coordinate, all straight lines through that coordinate are
// recorded on counters, one for each line. If a counter rises to two, it
// mesans that there are two hits who can cast that line, i.e. there are two
// hits on that line.
// A line can be represented by two independent parameters, therefore we need a
// two-dimension parameter space to represent all of them (two degrees of
// freedom).
// Since there are infinite straight lines, each counter represents not just
// one line but a small region of parameter space.
// Passing by one point constraints one of the two parameters, therefore the
// parameter set of a generic line passing by a given point has one degree of
// freedom.
// We follow the custom of choosing as parameters of an arbitrary line passing
// through a point (x, y) the angle (from "x" axis) of the line, and the
// distance of the line from the chosen point. Fixing the point, we can assign
// freedom to the angle (that has a limited range by definition) and then
// the distance will be defined by some functional form d(angle).
// The shape (d vs. a) of the set of parameters of all lines passing by a point
// is a sinusoidal curve. We can define the angle to be between 0 and pi
// (half a period) and the distance potentially covers the whole real axis
// (but actually it will be never larger than the distance of the selected
// point from the origin).
// 
// The implementation of the algorithm is based on a "accumulator", the set
// of counters of how many hits are passed by a given line (represented by
// its tow parameters). The accumulator is a two-dimensional container of
// counters, with the first dimension given by the angle and the second by the
// distance.
// In this scenario, all angles are sampled, therefore each angle will have
// at least one count per hit (somewhere at some distance d).
// Each angle will see one entry per hit.
// We choose to have a large number of sampled angles (the current standard
// configuration says 10800), therefore most of the counters will be empty.
// The sinusoidal shape can be steep enough that got example the angle a(1)
// has d(a1)=30 and the next angle (a2) has d(a2)=50. In this case we cover
// also the distances between 31 and 50, (assigning them, as an approximation,
// to a2). In this way we don't leave gaps and make sure that each two
// sinusoidal curves cross in at least one point, or, equivalently, that we
// can always find at least one straight line throw any two points.
// This translates in having very often the counts at a given angle clustered
// around some distances.
// 
// We need to discretize the angles and the distances. Tests show that for a
// "natural" plane size of 9600 (TDC counts) x ~3000 (wires), we need to use
// O(10000) angles and to oversample the distance, that typically goes in the
// range [0-10000], by some factor (5 in the default parameters).
// The oversampling factor is just an artifact of the fact that we use integral
// distances but we want to have a resolution better than "1" -- we just
// multiply the distance by e.g. 5 and we get that the distace "1" actually
// means 1/5 = 0.2.
// Note that the algorithm is usually applied to subsets of that plane rather
// than on the full plane, making the typical distance range somehow smaller.
// 
// The fastest data structure for the accumulator is a two-dimensional array.
// The proper size of this array would be some gigabyte, that makes this
// approach unfeasable. Since the dimension of angles has a very clear number
// of "bins" (covering always a fixed 0-pi range), but for each given angle
// the counters actually used are a few, a sparse structure (associative
// container, or "map") is used to describe all the counters at a given angle,
// with key the parameter r (discretized).
// Since all the angles always have data, the outer container is a vector
// whose index is the parameter a (also discretized).
// Therefore, for a given line (a, r), we can find how many hits pass though
// it by finding the associative container of the angle a from a vector,
// and in there looking for an entry with d as key: accum[a][d].
// 
// Optimization
// ----------------------------------------------------------------------------
// 
// Given the constraints and data structure described above, it turns out that
// a lot of time is spent creating new counters. On each hit, and each angle,
// one or more elements of the associative containers must be created.
// In total, millions of counter instances ar used at once.
// The standard C++ implementation of it, std::map, dynamically a new node
// each time a new counter is required, and in the end it frees them one by
// one. Both operations are very time-demanding.
// We use here a custom memory allocator (BulkAllocator) that prepares memory
// for chunks of nodes and then returns a preallocated space at each request
// for a new node. The allocator is designed to be fast, giving up features:
// nodes are never really freed, that saves a lot of book-keeping.
// A lot of tricky aspects of it are documented in its own header file.
// Also freeing the memory is fast, since it implies the deletion of a few
// (very large) chunks rather than millions.
// The other aspect taking a lot of time is the lookup of a counter: where is
// counter (a,d)? the location of a is fast (constant time, from a vector).
// The location of d is the next best thing, a binary tree with log(N) access
// time. Unfortunately a balanced tree needs to be rebalanced often, and that
// takes a lot of time. Also, N may be "small" (O(1000)), but there are still
// million insertions and look ups.
// We use here a replacement of the plain map. Since it often happens that,
// to "fill the gaps", sequential counters are allocated and increased,
// we have blocks of couners allocated together (in the current version, 64
// of them). This can save memory in case of crowded spaces: the overhead
// for each node is 40 bytes, whether it is a single counter or a block of
// 64). Look up within the same block becomes constant-time.
// Also, to access sequences of counters, special code is used so that we
// take advantage of the result from the previous look up to perform the next
// one, including also the insertion of a new counter block after an existing
// one.
// Finally, the algorithm acts when it finds a relatively small number of
// aligned hits, afterward removing the hits already clustered. The hit count
// keeps small all the time: we use a signed char as basic data type for
// the counters, allowing a maximum of 127 aligned hits. This saves a lot of
// memory, at the cost of a small slow-down for large (e.g. 64-bit) bus
// architectures. No check is performed for overflow; that can also be
// implemented at a small cost.
//
//
////////////////////////////////////////////////////////////////////////
#ifndef HOUGHBASEALG_H
#define HOUGHBASEALG_H

#include <vector>
#include <array>
#include <map>
#include <utility> // std::pair<>

#include <TMath.h>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 

#include "art/Persistency/Common/PtrVector.h" 
#include "lardata/Utilities/BulkAllocator.h"
#include "lardata/Utilities/CountersMap.h"

namespace art { class Event; }

namespace recob { 
  class Hit;
  class Cluster; 
}

    struct houghCorner
    {
      double strength=0;
      double p0=0;
      double p1=0;
      houghCorner(double strengthTemp=0,
          double p0Temp=0,
          double p1Temp=0)
      {
        strength=strengthTemp;
        p0=p0Temp;
        p1=p1Temp;
      }

      bool operator < (const houghCorner& houghCornerComp) const
      {
        return (strength < houghCornerComp.strength);
      }
    };


    // This stores information about merged lines
    struct mergedLines
    {
      double totalQ=0;
      double pMin0=0;
      double pMin1=0;
      double pMax0=0;
      double pMax1=0;
      int clusterNumber=-999999;
      double showerLikeness=0;
      mergedLines (double totalQTemp=0,
          double pMin0Temp=0,
          double pMin1Temp=0,
          double pMax0Temp=0,
          double pMax1Temp=0,
          double clusterNumberTemp=-999999,
          double showerLikenessTemp=0)
      {
        totalQ=totalQTemp;
        pMin0=pMin0Temp;
        pMin1=pMin1Temp;
        pMax0=pMax0Temp;
        pMax1=pMax1Temp;
        clusterNumber=clusterNumberTemp;
        showerLikeness=showerLikenessTemp;
      }
    };



    struct protoTrack
    {
      int clusterNumber=999999;
      int planeNumber=999999;
      int oldClusterNumber=999999;
      float clusterSlope=999999;
      float clusterIntercept=999999;
      float totalQ=-999999;
      float pMin0=999999;
      float pMin1=999999;
      float pMax0=-999999;
      float pMax1=-999999;
      float iMinWire=999999;
      float iMaxWire=-999999;
      float minWire=999999;
      float maxWire=-999999;
      float isolation=-999999;
      float showerLikeness=-999999;
      bool merged=false;
      bool showerMerged=false;
      bool mergedLeft=false;
      bool mergedRight=false;
      std::vector<art::Ptr<recob::Hit>> hits;
      protoTrack(){
      }
      
      void Init(unsigned int num=999999, 
	  unsigned int pnum=999999,
          float slope=999999, 
          float intercept=999999,
          float totalQTemp=-999999,
          float Min0=999999, 
          float Min1=999999, 
          float Max0=-999999, 
          float Max1=-999999,
          int    iMinWireTemp=999999,
          int    iMaxWireTemp=-999999,
          int    minWireTemp=999999,
          int    maxWireTemp=-999999,
          std::vector<art::Ptr<recob::Hit>> hitsTemp=std::vector<art::Ptr<recob::Hit>>())
      {
        clusterNumber = num;
        planeNumber = pnum;
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
   
   
  /**
   * @brief CountersMap with access optimized for Hough Transform algorithm
   * @param KEY the type of the key of the counters map
   * @param COUNTER the type of a basic counter (can be signed or unsigned)
   * @param BLOCKSIZE the number of counters in a cluster
   * @param ALLOC allocator for the underlying STL map
   * @param SUBCOUNTERS split each counter in subcounters (not implemented yet)
   * @see CountersMap
   *
   * In addition to the standard CountersMap interface, a special range
   * increment, increment with max detection and decrement methods are provided.
   */
  template <
    typename KEY,
    typename COUNTER,
    size_t SIZE,
    typename ALLOC = std::allocator<std::pair<KEY,std::array<COUNTER, SIZE>>>,
    unsigned int SUBCOUNTERS=1
    >
  class HoughTransformCounters:
    public lar::CountersMap<KEY, COUNTER, SIZE, ALLOC, SUBCOUNTERS>
  {
      public:

    /// This class
    using CounterMap_t
      = HoughTransformCounters<KEY, COUNTER, SIZE, ALLOC, SUBCOUNTERS>;
    /// Base class
    using Base_t = lar::CountersMap<KEY, COUNTER, SIZE, ALLOC, SUBCOUNTERS>;
    
    // import useful types
    using BaseMap_t = typename Base_t::BaseMap_t;
    using Allocator_t = typename Base_t::Allocator_t;
    using Key_t = typename Base_t::Key_t;
    using SubCounter_t = typename Base_t::SubCounter_t;
    using CounterBlock_t = typename Base_t::CounterBlock_t;
    using const_iterator = typename Base_t::const_iterator;
    
    /// Pair identifying a counter and its current value
    using PairValue_t = std::pair<const_iterator, SubCounter_t>;

    /// Default constructor (empty map)
    HoughTransformCounters(): Base_t() {}
    
    /// Constructor, specifies an allocator
    HoughTransformCounters(Allocator_t alloc): Base_t(alloc) {}
    
    
    /**
     * @brief Sets the specified counter to a count value
     * @param key key of the counter to be set
     * @param value the count value
     * @return new value of the counter
     */
    SubCounter_t set(Key_t key, SubCounter_t value)
      { return Base_t::set(key, value); }
    
    /**
     * @brief Increments by 1 the specified counter
     * @param key key of the counter to be increased
     * @return new value of the counter
     */
    SubCounter_t increment(Key_t key) { return Base_t::increment(key); }
    
    /**
     * @brief Decrements by 1 the specified counter
     * @param key key of the counter to be decreased
     * @return new value of the counter
     */
    SubCounter_t decrement(Key_t key) { return Base_t::decrement(key); }
    
    
    /**
     * @brief Sets the specified range of counters to a count value
     * @param key_begin key of the first counter to be set
     * @param key_end key of the first counter not to be set
     * @param value the count value
     * @return new value of all the counters
     * @see increment(), decrement(), increment_and_get_max()
     */
    SubCounter_t set(Key_t key_begin, Key_t key_end, SubCounter_t value)
      { return unchecked_set_range(key_begin, key_end, value); }
    
    /**
     * @brief Increments by 1 the specified range of counters
     * @param key_begin key of the first counter to be increased
     * @param key_end key of the first counter not to be increased
     * @see decrement(), increment_and_get_max()
     */
    void increment(Key_t key_begin, Key_t key_end);
    
    
    /**
     * @brief Increments by 1 the specified counters and returns the maximum
     * @param key_begin key of the first counter to be increased
     * @param key_end key of the first counter not to be increased
     * @return pair with an iterator to the largest counter and its value
     * @see increment(Key_t, Key_t)
     *
     * This method works like the corresponding increment() method, and in
     * addition it returns the location of the counter with the largest count
     * (after the increase).
     * The return value consist of a pair: the second is the largest counter
     * value in the range after the increase, the first member is the constant
     * iterator pointing to the first (lowest key) counter with that (get its
     * key with const_iterator::key() method).
     * 
     * If no maximum is found, the maximum in the return value is equal to
     * current_max, while the iterator points to the end of the map (end()).
     * Note that if all the counters are at the minimum possible value, no
     * maximum will be returned.
     */
    PairValue_t increment_and_get_max(Key_t key_begin, Key_t key_end)
      { return unchecked_add_range_max(key_begin, key_end, +1); }
    
    
    /**
     * @brief Increments by 1 the specified counters and returns the maximum
     * @param key_begin key of the first counter to be increased
     * @param key_end key of the first counter not to be increased
     * @param current_max only counters larger than this will be considered
     * @return pair with an iterator to the largest counter and its value
     * @see increment(Key_t, Key_t), increment_and_get_max(Key_t, Key_t)
     *
     * This method works like increment_and_get_max() method, except that it
     * does not update the maximum if it's not (strictly) larger than
     * current_max. If no such a maximum is found, the maximum in the return
     * value is equal to current_max, while the iterator points to the end of
     * the map (end()).
     */
    PairValue_t increment_and_get_max
      (Key_t key_begin, Key_t key_end, SubCounter_t current_max)
      { return unchecked_add_range_max(key_begin, key_end, +1, current_max); }
    
    
    /**
     * @brief Decrements by 1 the specified range of counters
     * @param key_begin key of the first counter to be increased
     * @param key_end key of the first counter not to be increased
     * @see increment()
     */
    void decrement(Key_t key_begin, Key_t key_end);
    
    
    /**
     * @brief Returns the largest counter
     * @param current_max only counters larger than this will be considered
     * @return pair with an iterator to the largest counter and its value
     * @see get_max(), increment_and_get_max(Key_t, Key_t)
     *
     * All counters are parsed, and the first one with the largest count
     * is returned.
     * The return value consist of a pair: the second is the largest counter
     * value, the first member is the constant iterator pointing to the first
     * (lowest key) counter with that (get its key with const_iterator::key()
     * method).
     * 
     * This method does not update the maximum if it's not (strictly) larger
     * than current_max. If no such a maximum is found, the maximum in the
     * return value is equal to current_max, while the iterator points to the
     * end of the map (end()).
     */
    PairValue_t get_max(SubCounter_t current_max) const;
    
    /**
     * @brief Increments by 1 the specified counters and returns the maximum
     * @param current_max only counters larger than this will be considered
     * @return pair with an iterator to the largest counter and its value
     * @see get_max(SubCounter_t), increment_and_get_max(Key_t, Key_t)
     *
     * This method works like get_max() method, except that it
     * does not update the maximum if it's not (strictly) larger than
     * current_max. If no such a maximum is found, the maximum in the return
     * value is equal to current_max, while the iterator points to the end of
     * the map (end()).
     */
    PairValue_t get_max() const;
    
    
      protected:
    
    using CounterKey_t = typename Base_t::CounterKey_t;
    
    
      private:
    
    SubCounter_t unchecked_set_range(
      Key_t key_begin, Key_t key_end, SubCounter_t value,
      typename BaseMap_t::iterator start
      );
    SubCounter_t unchecked_set_range
      (Key_t key_begin, Key_t key_end, SubCounter_t value);
    PairValue_t unchecked_add_range_max(
      Key_t key_begin, Key_t key_end, SubCounter_t delta,
      typename BaseMap_t::iterator start,
      SubCounter_t min_max = std::numeric_limits<SubCounter_t>::min()
      );
    PairValue_t  unchecked_add_range_max(
      Key_t key_begin, Key_t key_end, SubCounter_t delta,
      SubCounter_t min_max = std::numeric_limits<SubCounter_t>::min()
      );
    
  }; // class HoughTransformCounters
  
  
#define FC_DEVELOP 0
  
  class HoughTransform {
  public:
    
    HoughTransform();
    ~HoughTransform();
     
    void Init
      (unsigned int dx, unsigned int dy, float rhores, unsigned int numACells);
    std::array<int,3> AddPointReturnMax(int x, int y);
    bool SubtractPoint(int x, int y);
    int  GetCell(int row, int col) const;
    void SetCell(int row, int col, int value) { m_accum[row].set(col, value); }
    void GetAccumSize(int &numRows, int &numCols) 
    { 
      numRows = m_accum.size();
      numCols  = (int) m_rowLength;
    }
    int NumAccumulated()                      { return m_numAccumulated; }
    void GetEquation( float row, float col, float &rho, float &theta) const;
    int GetMax(int & xmax, int & ymax) const;

    void reconfigure(fhicl::ParameterSet const& pset);

  private:
    
    /// rho -> # hits (for convenience)
    typedef HoughTransformCounters<int, signed char, 64> BaseMap_t;
    
    /// Special allocator for large chunks of pairs (turns out map won't use it)
    typedef lar::BulkAllocator<BaseMap_t::allocator_type::value_type>
      BulkPairAllocator_t;
    
    /// Type of map distance (discretized) =># hits,
    /// #hits stored in counters allocated in blocks
    typedef HoughTransformCounters<int, signed char, 64, BulkPairAllocator_t>
      DistancesMap_t;
    
    /// Type of the Hough transform (angle, distance) map with custom allocator
    typedef std::vector<DistancesMap_t> HoughImage_t;
    
    
    unsigned int m_dx;
    unsigned int m_dy;
    unsigned int m_rowLength;
    unsigned int m_numAngleCells;
    float m_rhoResolutionFactor;
    // Note, m_accum is a vector of associative containers,
    // the vector elements are called by rho, theta is the container key,
    // the number of hits is the value corresponding to the key
    HoughImage_t m_accum;  ///< column (map key)=rho, row (vector index)=theta
    int m_numAccumulated;
    std::vector<double> m_cosTable;
    std::vector<double> m_sinTable;
    
    std::array<int,3> DoAddPointReturnMax(int x, int y, bool bSubtract = false);


  }; // class HoughTransform  





  class HoughBaseAlg {
    
  public:
    
    /// Data structure collecting charge information to be filled in cluster
    struct ChargeInfo_t {
      float integral         = 0.0F;
      float integral_stddev  = 0.0F;
      float summedADC        = 0.0F;
      float summedADC_stddev = 0.0F;
      
      ChargeInfo_t(float in, float in_stdev, float sum, float sum_stdev):
        integral(in), integral_stddev(in_stdev),
        summedADC(sum), summedADC_stddev(sum_stdev)
        {}
    }; // ChargeInfo_t
    
    
    HoughBaseAlg(fhicl::ParameterSet const& pset); 
    virtual ~HoughBaseAlg();

    size_t FastTransform(const std::vector<art::Ptr<recob::Cluster> >         & clusIn,
			 std::vector<recob::Cluster>                    & ccol,  
			 std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
			 art::Event                                const& evt,
			 std::string                               const& label);

    size_t Transform(std::vector<art::Ptr<recob::Hit> > const& hits,
                     std::vector<unsigned int>     *fpointId_to_clusterId,
                     unsigned int clusterId, // The id of the cluster we are examining
                     unsigned int *nClusters,
                     std::vector<protoTrack> *protoTracks);
    
    
    // interface to look for lines only on a set of hits,without slope and totalQ arrays
    size_t FastTransform(
      std::vector<art::Ptr<recob::Hit>>      & clusIn,
      std::vector<art::PtrVector<recob::Hit>>& clusHitsOut
      );
    
    // interface to look for lines only on a set of hits
    size_t FastTransform(
      std::vector<art::Ptr<recob::Hit>>      & clusIn,
      std::vector<art::PtrVector<recob::Hit>>& clusHitsOut,
      std::vector<double>                    & slope,
      std::vector<ChargeInfo_t>              & totalQ
      );
    

    size_t Transform(std::vector<art::Ptr<recob::Hit> > const& hits);

    size_t Transform(std::vector< art::Ptr<recob::Hit> > const& hits,
		     double                                   & slope,
		     double                                   & intercept);

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
    float  fMaxDistance;                   ///< Max distance that a hit can be from a line to be considered part of that line
    float  fMaxSlope;                      ///< Max slope a line can have
    int    fRhoZeroOutRange;               ///< Range in rho over which to zero out area around previously found lines in the accumulator
    int    fThetaZeroOutRange;             ///< Range in theta over which to zero out area around previously found lines in the accumulator
    float  fRhoResolutionFactor;           ///< Factor determining the resolution in rho
    int    fPerCluster;                    ///< Tells the original Hough algorithm to look at clusters individually, or all hits
                                           ///< at once
    int    fMissedHits;                    ///< Number of wires that are allowed to be missed before a line is broken up into
                                           ///< segments
    float  fMissedHitsDistance;            ///< Distance between hits in a hough line before a hit is considered missed
    float  fMissedHitsToLineSize;          ///< Ratio of missed hits to line size for a line to be considered a fake

  protected:

    friend class HoughTransformClus;
  };
  
  
}// namespace

#endif // HOUGHBASEALG_H
