/**
 *  @file   IHitReader.h
 *
 *  @brief  This provides an art tool interface definition for reading hits into the SpacePointSolver universe
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IHitReader_h
#define IHitReader_h

// Algorithm includes
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"

#include <vector>

//------------------------------------------------------------------------------------------------------------------------------------------
namespace reco3d {
  /**
 *  @brief  IHitReader interface class definiton
 */
  class IHitReader {
  public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IHitReader() noexcept = default;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    /**
     *  @brief Scan an input collection of clusters and modify those according
     *         to the specific implementing algorithm
     *
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    virtual bool readHits(const std::vector<art::Ptr<recob::Hit>>&,      // input hits
                          std::vector<art::Ptr<recob::Hit>>&,            // output hits plane 0
                          std::vector<art::Ptr<recob::Hit>>&,            // output hits plane 1
                          std::vector<art::Ptr<recob::Hit>>&) const = 0; // output hits plane 2
  };

} // namespace lar_cluster3d
#endif
