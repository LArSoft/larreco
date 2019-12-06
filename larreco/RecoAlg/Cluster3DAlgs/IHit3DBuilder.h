/**
 *  @file   IHit3DBuilder.h
 *
 *  @brief  This provides an art tool interface definition for tools which construct 3D hits used in 3D clustering
 *          and outputs a new hit collection based on those 3D hits
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IHit3DBuilder_h
#define IHit3DBuilder_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace art
{
    class ProducesCollector;
}

namespace lar_cluster3d
{
/**
 *  @brief  IHit3DBuilder interface class definiton
 */
class IHit3DBuilder
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IHit3DBuilder() noexcept = default;

    /**
     *  @brief The space point building should output the hit collection
     *         for those hits which combine to form space points - a nice noise filter!
     */
    virtual void produces(art::ProducesCollector&) = 0;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    /**
     *  @brief Defines a structure mapping art representation to internal
     */
    using RecobHitToPtrMap = std::unordered_map<const recob::Hit*, art::Ptr<recob::Hit>>;

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    virtual void Hit3DBuilder(art::Event&, reco::HitPairList&, RecobHitToPtrMap&) = 0;

    /**
     *  @brief enumerate the possible values for time checking if monitoring timing
     */
    enum TimeValues {COLLECTARTHITS   = 0,
                     BUILDTHREEDHITS  = 1,
                     BUILDNEWHITS     = 2,
                     NUMTIMEVALUES
    };

    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    virtual float getTimeToExecute(TimeValues index) const = 0;

};

} // namespace lar_cluster3d
#endif
