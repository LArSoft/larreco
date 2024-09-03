////////////////////////////////////////////////////////////////////////
//
//  Cosmic removal tools
//
//  Tingjun Yang (tjyang@fnal.gov)
//
////////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERTCCR_H
#define TRAJCLUSTERTCCR_H

namespace detinfo {
  class DetectorClocksData;
}

namespace tca {

  struct PFPStruct;
  struct TCSlice;

  void SaveCRInfo(detinfo::DetectorClocksData const& clockData,
                  TCSlice& tcs,
                  PFPStruct& ms,
                  bool fIsRealData);
  int GetOrigin(detinfo::DetectorClocksData const& clockData, TCSlice& tcs, PFPStruct& ms);
  void ClearCRInfo(TCSlice& tcs);
}

#endif
