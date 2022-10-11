/**
 *  @file   INormalizeCharge.h
 *
 *  @brief  This is an interface for an art Tool which scales charge by some
 *          factor given information about its associated hit.
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef INormalizeCharge_h
#define INormalizeCharge_h

// Framework Includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

/**
 *  @brief  INormalizeCharge interface class definiton
 */
class INormalizeCharge {
public:
  /**
     *  @brief  Virtual Destructor
     */
  virtual ~INormalizeCharge() noexcept = default;

  virtual void configure(const fhicl::ParameterSet&) = 0;
  virtual double Normalize(double dQdx,
                           const art::Event& e,
                           const recob::Hit& h,
                           const geo::Point_t& location,
                           const geo::Vector_t& direction,
                           double t0) = 0;
};

#endif
