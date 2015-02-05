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

#include "RegionAboveThresholdFinder.h"
#include <stdexcept>

void hit::RegionAboveThresholdFinder::FillStartAndEndTicks(const std::vector<float>& signal,
							   std::vector<unsigned int>& start_ticks,
							   std::vector<unsigned int>& end_ticks)
{

  start_ticks.clear(); end_ticks.clear();

  bool in_RAT = false;
  for(unsigned int i_tick=0; i_tick<signal.size(); i_tick++){

    if(!in_RAT && signal[i_tick]>=fThreshold){
      start_ticks.push_back(i_tick);
      in_RAT = true;
    }
    else if(in_RAT && signal[i_tick]<fThreshold){
      end_ticks.push_back(i_tick);
      in_RAT = false;
    }

  }

  if(in_RAT)
    end_ticks.push_back(signal.size());

  if(end_ticks.size()!=start_ticks.size())
    throw std::runtime_error("ERROR in RegionAboveThresholdFinder: start and end tick vectors not equal.");

}
