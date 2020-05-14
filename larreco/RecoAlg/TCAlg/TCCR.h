////////////////////////////////////////////////////////////////////////
//
//  Cosmic removal tools
//
//  Tingjun Yang (tjyang@fnal.gov)
//
////////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERTCCR_H
#define TRAJCLUSTERTCCR_H

namespace tca {

  struct PFPStruct;
  struct TCSlice;

  void SaveCRInfo(TCSlice& tcs, PFPStruct& ms, bool prt, bool fIsRealData);
  int GetOrigin(TCSlice& tcs, PFPStruct& ms);
  void ClearCRInfo(TCSlice& tcs);
}

#endif
