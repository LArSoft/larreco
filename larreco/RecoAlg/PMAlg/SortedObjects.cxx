/**
 *  @file   SortedObjects.cxx
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Base classes for chains and branched chains of nodes and segments.
 *          See PmaTrack3D.h file for details.
 */

#include "larreco/RecoAlg/PMAlg/SortedObjects.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

//***********************  SortedObjectBase  ***********************
pma::SortedObjectBase::SortedObjectBase(pma::SortedObjectBase* prevElement,
                                        pma::SortedObjectBase* nextElement)
  : next(0), prev(0)
{
  if (prevElement) prevElement->AddNext(this);
  if (nextElement) AddNext(nextElement);
}

void pma::SortedObjectBase::Disconnect(void)
{
  if (prev) prev->RemoveNext(this);
  if (next) RemoveNext(next);
}

bool pma::SortedObjectBase::AddNext(pma::SortedObjectBase* nextElement)
{
  if (!nextElement || (nextElement == this)) return false;

  if (next && (next->prev == this)) next->prev = 0;

  if (nextElement->prev && (nextElement->prev != this)) nextElement->prev->RemoveNext(nextElement);

  next = nextElement;
  next->prev = this;
  return true;
}

int pma::SortedObjectBase::RemoveNext(pma::SortedObjectBase* nextElement)
{
  if (nextElement && (next == nextElement)) {
    if (next->prev == this)
      next->prev = 0;
    else {
      mf::LogError("pma::SortedObjectBase") << "Object structure is broken.";
    }

    next = 0;
    return 0;
  }
  else
    return -1;
}
//******************************************************************

//***********************  SortedBranchBase  ***********************
void pma::SortedBranchBase::Disconnect(void)
{
  while (next_vector.size())
    RemoveNext(next_vector.front());
  if (prev) prev->RemoveNext(this);
}

bool pma::SortedBranchBase::AddNext(pma::SortedObjectBase* nextElement)
{
  if (!nextElement) {
    mf::LogError("pma::SortedBranchBase") << "Next == 0.";
    return false;
  }

  if (nextElement == this) {
    mf::LogWarning("pma::SortedBranchBase") << "Next == This.";
    return false;
  }

  bool present = false;
  for (size_t i = 0; i < next_vector.size(); i++) {
    if (next_vector[i] == nextElement) {
      mf::LogWarning("pma::SortedBranchBase") << "Contained.";
      present = true;
      break;
    }
  }
  if (!present) {
    if (nextElement->prev) // && (nextElement->prev != this)
      nextElement->prev->RemoveNext(nextElement);

    next = nextElement;
    next->prev = this;
    next_vector.push_back(next);
    return true;
  }
  else
    return false;
}

int pma::SortedBranchBase::RemoveNext(pma::SortedObjectBase* nextElement)
{
  if (!nextElement || (nextElement == this)) return -1;

  int index = -1;
  for (unsigned int i = 0; i < next_vector.size(); i++) {
    if (next_vector[i] == nextElement) {
      index = i;
      break;
    }
  }
  if (index >= 0) {
    if (next_vector[index]->prev != this) {
      mf::LogError("pma::SortedBranchBase") << "Object structure is broken.";
      return -1;
    }
    next_vector[index]->prev = 0;
    next_vector.erase(next_vector.begin() + index);

    if (next_vector.size() > 0)
      next = next_vector.back();
    else
      next = 0;
  }
  return index;
}
//******************************************************************
