////////////////////////////////////////////////////////////////////////
//
// LineMerger class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// biagio.rossi@lhep.unibe.ch
// msoderbe@syr.edu
// joshua.spitz@yale.edu
//
// This algorithm is designed to merge 2D lines with similar slope and endpoints 
//  
////////////////////////////////////////////////////////////////////////

#include <string>
#include <cmath> // std::abs(), std::sqrt()
#include <iomanip>
#include <vector>
#include <array>
#include <memory> // std::unique_ptr<>
#include <utility> // std::move()


//Framework includes:
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

//LArSoft includes:
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "larreco/ClusterFinder/ClusterCreator.h"



//#ifndef LINEMERGER_H
//#define LINEMERGER_H


namespace cluster {
  
  class LineMerger : public art::EDProducer {
    
  public:
    
    explicit LineMerger(fhicl::ParameterSet const& pset); 
    ~LineMerger();
    
    void produce(art::Event& evt);
    void beginJob();
    
  private:
        
    std::string     fClusterModuleLabel;
    double          fSlope; // tolerance for matching angles between two lines (in units of radians) 
    double          fEndpointWindow; // tolerance for matching endpoints (in units of time samples) 
   
    bool SlopeCompatibility(double slope1,double slope2);
    int  EndpointCompatibility(
      float sclstartwire, float sclstarttime,
      float sclendwire,   float sclendtime,
      float cl2startwire, float cl2starttime,
      float cl2endwire,   float cl2endtime
      );
    
  protected: 
    
  }; // class LineMerger

}

//#endif // LINEMERGER_H



namespace cluster{
  
  
  /// Class merging clusters: recomputes start and end position and hit list
  class ClusterMerger {
      public:
    // NOTE if you feel like copying this class, move it into its own header
    // instead, and if you need it for a different hit or hit pointer, make it
    // a template instead
    using HitPtr_t = art::Ptr<recob::Hit>; ///< type of pointer to hits
    using HitVector_t = std::vector<HitPtr_t>; ///< vector of pointers to hits
    
    using ID_t = recob::Cluster::ID_t ; ///< Type of cluster ID
    using ClusterEnds_t = recob::Cluster::ClusterEnds_t;
    
    /*
      typedef enum {
        clStart,       ///< Represents the most likely start of the cluster
        clEnd,         ///< Represents the end, or the alternative start, of the cluster
        NEnds,         ///< End count
        clFirstEnd = 0 ///< Just an alias for loops
      } ClusterEnds_t; ///< Used to decide which end to use
    */
    
    ClusterMerger() = default;
    
    ClusterMerger(recob::Cluster const& cluster):
      ClusterMerger()
      { Add(cluster); }
    
    /**
     * @brief Merges a single cluster into this object
     * @param cluster the cluster to be merged
     * @return whether the addition was successful
     * 
     * The two ends of the cluster are merged into this one, that gets extended.
     * 
     * The new cluster must have the same view as the prevopus ones and must lay
     * on the same plane.
     * If the new cluster has invalid plane, the current one is kept; if the
     * current plane is invalid, it is overwritten (that means that if both are
     * invalid, the merged cluster will also have an invalid plane).
     * 
     * Note that this code is crap unless the cluster is track-like.
     */
    bool Add(recob::Cluster const& cluster);
    
    
    /// @{
    /// @name Accessors
    
    /// Returns the wire coordinate of the start of the cluster
    float StartWire() const { return fEndWires[ClusterEnds_t::clStart]; }
    
    /// Returns the tick coordinate of the start of the cluster
    float StartTick() const { return fEndTicks[ClusterEnds_t::clStart]; }
    
    /// Returns the uncertainty on wire coordinate of the start of the cluster
    float SigmaStartWire() const { return fSigmaEndWires[ClusterEnds_t::clStart]; }
    
    // Returns the uncertainty on tick coordinate of the start of the cluster
    float SigmaStartTick() const { return fSigmaEndTicks[ClusterEnds_t::clStart]; }
    
    /// Returns the wire coordinate of the end of the cluster
    float EndWire() const { return fEndWires[ClusterEnds_t::clEnd]; }
    
    /// Returns the tick coordinate of the end of the cluster
    float EndTick() const { return fEndTicks[ClusterEnds_t::clEnd]; }
    
    /// Returns the uncertainty on wire coordinate of the end of the cluster
    float SigmaEndWire() const { return fSigmaEndWires[ClusterEnds_t::clEnd]; }
    
    /// Returns the uncertainty on tick coordinate of the end of the cluster
    float SigmaEndTick() const { return fSigmaEndTicks[ClusterEnds_t::clEnd]; }
    
    /// Returns the wire coordinate of one of the end sides of the cluster
    float WireCoord(ClusterEnds_t side) const { return fEndWires[side]; }
  //  float WireCoord(unsigned int side) const { return fEndWires[side]; }
    
    /// Returns the tick coordinate of one of the end sides of the cluster
    float TickCoord(ClusterEnds_t side) const { return fEndTicks[side]; }
  //  float TickCoord(unsigned int side) const { return fEndTicks[side]; }
    
    /// Returns the uncertainty on wire coordinate of one of the end sides of the cluster
    float SigmaWireCoord(ClusterEnds_t side) const { return fSigmaEndWires[side]; }
  //  float SigmaWireCoord(unsigned int side) const { return fSigmaEndWires[side]; }
    
    /// Returns the uncertainty on tick coordinate of one of the end sides of the cluster
    float SigmaTickCoord(ClusterEnds_t side) const { return fSigmaEndTicks[side]; }
  //  float SigmaTickCoord(unsigned int side) const { return fSigmaEndTicks[side]; }
    
    
    /// Returns the charge on the first wire of the cluster
    float StartCharge() const { return fEndCharges[ClusterEnds_t::clStart]; }
    
    /// Returns the starting angle of the cluster
    float StartAngle() const { return fAngles[ClusterEnds_t::clStart]; }
    
    /// Returns the opening angle at the start of the cluster
    float StartOpeningAngle() const { return fOpeningAngles[ClusterEnds_t::clStart]; }
    
    /// Returns the charge on the last wire of the cluster
    float EndCharge() const { return fEndCharges[ClusterEnds_t::clEnd]; }
    
    /// Returns the ending angle of the cluster
    float EndAngle() const { return fAngles[ClusterEnds_t::clEnd]; }
    
    /// Returns the opening angle at the end of the cluster
    float EndOpeningAngle() const { return fOpeningAngles[ClusterEnds_t::clEnd]; }
    
    /// Returns the charge on the first or last wire of the cluster
    float EdgeCharge(ClusterEnds_t side) const { return fEndCharges[side]; }
  //  float EdgeCharge(unsigned int side) const { return fEndCharges[side]; }
    
    /// Returns the angle at either end of the cluster
    float Angle(ClusterEnds_t side) const { return fAngles[side]; }
  //  float Angle(unsigned int side) const { return fAngles[side]; }
    
    /// Returns the opening angle at either end of the cluster
    float OpeningAngle(ClusterEnds_t side) const { return fOpeningAngles[side]; }
  //  float OpeningAngle(unsigned int side) const { return fOpeningAngles[side]; }
    
    
    /// A measure of the cluster width, in homogenized units.
    float Width() const { return fWidth; }
    
    /// Returns the view for this cluster
    geo::View_t View() const { return fView; }
    
    /// Returns the plane ID this cluster lies on
    geo::PlaneID Plane() const { return fPlaneID; }
    
    /// Returns whether geometry plane is valid
    bool hasPlane() const { return Plane().isValid; }
    
    
      protected:
    
    /// Data referring to start and end of the cluster
    float fEndWires[ClusterEnds_t::NEnds];
      
    /// Uncertainty on wire coordinate of the start and end of the cluster
    float fSigmaEndWires[ClusterEnds_t::NEnds];
    
    /// Tick coordinate of the start and end of the cluster
    float fEndTicks[ClusterEnds_t::NEnds];
    
    /// Uncertainty on tick coordinate of the start and end of the cluster
    float fSigmaEndTicks[ClusterEnds_t::NEnds];
    
    /// Charge on the start and end wire of the cluster.
    float fEndCharges[ClusterEnds_t::NEnds];
    
    /// Angle of the start and end of the cluster, defined in [-pi,pi]
    float fAngles[ClusterEnds_t::NEnds];
    
    /// Opening angle of the cluster shape at the start and end of the cluster.
    float fOpeningAngles[ClusterEnds_t::NEnds];
    
    /// A measure of the cluster width, in homogenized units.
    float fWidth;
    
    geo::View_t fView; ///< View for this cluster
    
    geo::PlaneID fPlaneID; ///< Location of the start of the cluster
    
    unsigned int n_clusters = 0; ///< number of clusters added so far
    
    /// Imports all the member of the corresponding end
    void AdoptEnd(recob::Cluster const& cluster, ClusterEnds_t iEnd);
    
    template <typename T>
    static void top(T& var, T value) { if (value > var) var = value; }
    template <typename T>
    static void bot(T& var, T value) { if (value < var) var = value; }
  }; // class ClusterMerger
  
  void ClusterMerger::AdoptEnd
    (recob::Cluster const& cluster, ClusterEnds_t iSrcEnd)
  {
    const ClusterEnds_t iDestEnd = iSrcEnd;
    fEndWires[iDestEnd]      = cluster.WireCoord(iSrcEnd);
    fSigmaEndWires[iDestEnd] = cluster.SigmaWireCoord(iSrcEnd);
    fEndTicks[iDestEnd]      = cluster.TickCoord(iSrcEnd);
    fSigmaEndTicks[iDestEnd] = cluster.SigmaTickCoord(iSrcEnd);
    fEndCharges[iDestEnd]    = cluster.EdgeCharge(iSrcEnd);
    fAngles[iDestEnd]        = cluster.Angle(iSrcEnd);
    fOpeningAngles[iDestEnd] = cluster.OpeningAngle(iSrcEnd);
  } // ClusterMerger::AdoptEnd()
  
  bool ClusterMerger::Add(recob::Cluster const& cluster) {
    if (!cluster.isValid()) return false;
    
    if (n_clusters == 0) { // special case: we are still empty
      AdoptEnd(cluster, ClusterEnds_t::clStart);
      AdoptEnd(cluster, ClusterEnds_t::clEnd);
      fWidth           = cluster.Width();
      fView            = cluster.View();
      fPlaneID         = cluster.Plane();
      ++n_clusters;
      return true;
    } // if empty
    
    if (cluster.View() != View()) return false;
    
    if (cluster.hasPlane() && hasPlane() && (cluster.Plane() != Plane()))
      return false;
    
    // this code has been moved here from the old recon::Cluster::operator+
    // of recob::Cluster v13.
    if (cluster.StartWire() < StartWire()) { // adopt the new start
      AdoptEnd(cluster, ClusterEnds_t::clStart);
    }
    if (cluster.EndWire() < EndWire()) { // adopt the new end
      AdoptEnd(cluster, ClusterEnds_t::clEnd);
    }
    
    top(fWidth, cluster.Width()); // extend width
    
    if (!hasPlane()) fPlaneID = cluster.Plane();
    
    return true;
  } // ClusterMerger::Add(Cluster)
  
  
  /// Class merging clusters: recomputes start and end position and hit list
  class ClusterAndHitMerger: public ClusterMerger {
      public:
    // NOTE if you feel like copying this class, move it into its own header
    // instead, and if you need it for a different hit or hit pointer, make it
    // a template instead
    using HitPtr_t = art::Ptr<recob::Hit>; ///< type of pointer to hits
    using HitVector_t = std::vector<HitPtr_t>; ///< vector of pointers to hits
    
    ClusterAndHitMerger() = default;
    
    ClusterAndHitMerger
      (recob::Cluster const& cluster, HitVector_t const& cluster_hits)
      { Add(cluster, cluster_hits); }
    
    /**
     * @brief Merges a single cluster into this object
     * @param cluster the cluster to be merged
     * @param cluster_hits the list of hits in this cluster
     * @param prepend if true, hits are inserted at the beginning of the list
     * @return whether the addition was successful
     * @see ClusterMerger::Add()
     * 
     * The two ends of the cluster are merged into this one, that gets extended.
     * Hit lists are merged too: no check on existing hits nor double addition.
     * 
     * Note that this code is crap unless the cluster is track-like.
     */
    bool Add(
      recob::Cluster const& cluster, HitVector_t const& cluster_hits,
      bool prepend = false
      );
    
    /// @{
    /// @name Accessors
    
    /// Returns a constant reference to the current list of hits
    HitVector_t const& Hits() const { return hits; }
    
    /// Number of hits in the cluster
    unsigned int NHits() const { return hits.size(); }
    
    ///@}
    
      protected:
    
    HitVector_t hits; ///< hits in the cluster
    
    void AddHits(HitVector_t const& cluster_hits, bool prepend)
      {
        hits.insert(prepend? hits.begin(): hits.end(),
          cluster_hits.begin(), cluster_hits.end());
      } // AddHits()
    
  }; // class ClusterAndHitMerger
  
  bool ClusterAndHitMerger::Add(
    recob::Cluster const& cluster, HitVector_t const& cluster_hits,
    bool prepend /* = false */
  ) {
    if (!ClusterMerger::Add(cluster)) return false;
    
    AddHits(cluster_hits, prepend);
    return true;
  } // ClusterAndHitMerger::Add()
  
  
  //-------------------------------------------------
  LineMerger::LineMerger(fhicl::ParameterSet const& pset) 
    : fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel"))
    , fSlope             (pset.get<double     >("Slope"))
    , fEndpointWindow    (pset.get<double     >("EndpointWindow"))
  {
    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  }

  //-------------------------------------------------
  LineMerger::~LineMerger()
  {
  }

  //-------------------------------------------------
  void LineMerger::beginJob()
  {
    //this doesn't do anything now, but it might someday
  }
    
  //------------------------------------------------------------------------------------//
  void LineMerger::produce(art::Event& evt)
  { 
    // Get a Handle for the input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    evt.getByLabel(fClusterModuleLabel,clusterVecHandle);

    constexpr size_t nViews = 3; // number of views we map
    
    //one vector for each view in the geometry (holds the index of the cluster)
    std::array< std::vector<size_t>, nViews > ClsIndices;

    //vector with indicators for whether a cluster has been merged already
    std::array< std::vector<int>, nViews > Cls_matches;

    // loop over the input Clusters
    for(size_t i = 0; i < clusterVecHandle->size(); ++i){
      
      //get a art::Ptr to each Cluster
      art::Ptr<recob::Cluster> cl(clusterVecHandle, i);
      
      size_t view = 0;
      switch(cl->View()){
      case geo::kU :
        view = 0;
        break;
      case geo::kV :
        view = 1;
        break;
      case geo::kZ :
        view = 2;
        break;
      default :
        continue; // ignore this cluster and process the next one
      }// end switch on view
      
      Cls_matches[view].push_back(0);
      ClsIndices[view].push_back(i);
    }// end loop over input clusters

    std::unique_ptr<std::vector<recob::Cluster> >             SuperClusters(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here; configuration would happen here,
    // but we are using the default configuration for that algorithm
    ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;
    
    art::FindManyP<recob::Hit> fmh(clusterVecHandle, evt, fClusterModuleLabel);
    
    for(size_t i = 0; i < nViews; ++i){

      int clustersfound = 0; // how many merged clusters found in each plane
      int clsnum1       = 0;

      for(size_t c = 0; c < ClsIndices[i].size(); ++c){
        if(Cls_matches[i][clsnum1] == 1){
          ++clsnum1;
          continue;
        }

        // make a new cluster to put into the SuperClusters collection 
        // because we want to be able to adjust it later;
        // use the hits associated with the current cluster
        recob::Cluster const& StartingCluster
          = clusterVecHandle->at(ClsIndices[i][c]);
        ClusterAndHitMerger cl1(StartingCluster, fmh.at(ClsIndices[i][c]));
        const recob::Cluster::ID_t clusterID = StartingCluster.ID();

        Cls_matches[i][clsnum1] = 1; 
        ++clustersfound;
        
        int clsnum2 = 0;
        for(size_t c2 = 0; c2 < ClsIndices[i].size(); ++c2){

          if(Cls_matches[i][clsnum2] == 1){
            ++clsnum2;
            continue;
          }

          const recob::Cluster& cl2( clusterVecHandle->at(ClsIndices[i][c2]) );

          
          // check that the slopes are the same
          // added 13.5 ticks/wirelength in ArgoNeuT. 
          // \todo need to make this detector agnostic
          // would be nice to have a LArProperties function that returns ticks/wire.
          bool sameSlope = SlopeCompatibility(cl1.StartAngle(), cl2.EndAngle())
            || SlopeCompatibility(cl1.EndAngle(), cl2.StartAngle());
          
          // check that the endpoints fall within a circular window of each other 
          // done in place of intercept matching
          int sameEndpoint = EndpointCompatibility(
            cl1.StartWire(), cl1.StartTick(),
            cl1.EndWire(),   cl1.EndTick(),
            cl2.StartWire(), cl2.StartTick(),
            cl2.EndWire(),   cl2.EndTick()
            );
          
          // if the slopes and end points are the same, combine the clusters
          // note that after 1 combination cl1 is no longer what we started 
          // with
          if(sameSlope && (sameEndpoint != 0)) {
            // combine the hit collections too
            // (find the hits associated with this second cluster);
            // take into account order when merging hits from two clusters: doc-1776
            // if sameEndpoint is 1, the new hits come first
            cl1.Add(cl2, fmh.at(ClsIndices[i][c2]), sameEndpoint == 1);
            Cls_matches[i][clsnum2] = 1;
          }
          
          ++clsnum2;
        }// end loop over second cluster iterator

        // now add the final version of cl1 to the collection of SuperClusters
        // and create the association between the super cluster and the hits
        ClusterParamAlgo.ImportHits(cl1.Hits());
        
        // create the recob::Cluster directly in the vector
        SuperClusters->emplace_back(
          cl1.StartWire(),                            // start_wire
          cl1.SigmaStartWire(),                       // sigma_start_wire
          cl1.StartTick(),                            // start_tick
          cl1.SigmaStartTick(),                       // sigma_start_tick
          cl1.StartCharge(),                          // start_charge
          cl1.StartAngle(),                           // start_angle
          cl1.StartOpeningAngle(),                    // start_opening
          cl1.EndWire(),                              // end_wire
          cl1.SigmaEndWire(),                         // sigma_end_wire
          cl1.EndTick(),                              // end_time
          cl1.SigmaEndTick(),                         // sigma_end_tick
          cl1.EndCharge(),                            // end_charge
          cl1.EndAngle(),                             // end_angle
          cl1.EndOpeningAngle(),                      // end_opening
          ClusterParamAlgo.Integral().value(),        // integral
          ClusterParamAlgo.IntegralStdDev().value(),  // integral_stddev
          ClusterParamAlgo.SummedADC().value(),       // summedADC
          ClusterParamAlgo.SummedADCStdDev().value(), // summedADC_stddev
          ClusterParamAlgo.NHits(),                   // n_hits
          ClusterParamAlgo.MultipleHitDensity(),        // multiple_hit_density
          cl1.Width(),                                // width
          clusterID,                                  // ID
          cl1.View(),                                 // view
          cl1.Plane(),                                // planeID
          recob::Cluster::Sentry                      // sentry
          );
        
        util::CreateAssn(*this, evt, *(SuperClusters.get()), cl1.Hits(), *(assn.get()));
        ++clsnum1;

      }// end loop over first cluster iterator
    }// end loop over planes

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "LineMerger Summary:";
    for(size_t i = 0; i < SuperClusters->size(); ++i) 
      mf::LogVerbatim("Summary") << SuperClusters->at(i);

    evt.put(std::move(SuperClusters));
    evt.put(std::move(assn));

    return;

  }

  //------------------------------------------------------------------------------------//
  //checks the difference between angles of the two lines
  bool LineMerger::SlopeCompatibility(double slope1, double slope2)
  { 
    double sl1 = atan(slope1);
    double sl2 = atan(slope2);

    //the units of fSlope are radians
    bool comp  = std::abs(sl1-sl2) < fSlope ? true : false;

    return comp;
  }
  //------------------------------------------------------------------------------------//
  int LineMerger::EndpointCompatibility(
    float sclstartwire, float sclstarttime,
    float sclendwire,   float sclendtime,
    float cl2startwire, float cl2starttime,
    float cl2endwire,   float cl2endtime
  ) {

    /// \todo 13.5 ticks/wire. need to make this detector agnostic--spitz
    float distance = std::sqrt((pow(sclendwire-cl2startwire,2)*13.5) + pow(sclendtime-cl2starttime,2));

    //not sure if this line is necessary--spitz
    float distance2 = std::sqrt((pow(sclstartwire-cl2endwire,2)*13.5) + pow(sclstarttime-cl2endtime,2));
    
//    bool comp = (distance  < fEndpointWindow ||
//                 distance2 < fEndpointWindow) ? true : false;

    //determine which way the two clusters should be merged. TY
    int comp = 0;
    if (distance < fEndpointWindow) 
      comp = 1;
    else if (distance2 < fEndpointWindow)
      comp = -1;
    return comp;
  }



} // end namespace





namespace cluster{

  DEFINE_ART_MODULE(LineMerger)
  
} // end namespace 

