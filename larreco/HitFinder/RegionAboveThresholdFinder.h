#ifndef REGIONABOVETHRESHOLDFINDER_H
#define REGIONABOVETHRESHOLDFINDER_H

/*!
 * Title:   RegionAboveThresholdFinder Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Class that finds a region above threshold in which to do
 *              hit-finding.
 *
 * Input:  Vector of floats (like a recob::Wire vector)
 * Output: Vector of begin times, and vector of end times.
*/

#include <vector>

namespace hit{

  class RegionAboveThresholdFinder {

  public:
    RegionAboveThresholdFinder(float threshold) { fThreshold = threshold; }

    void FillStartAndEndTicks(const std::vector<float>& signal,
			      std::vector<unsigned int>& start_ticks,
			      std::vector<unsigned int>& end_ticks);

  private:

    float fThreshold;

  };

}

#endif
