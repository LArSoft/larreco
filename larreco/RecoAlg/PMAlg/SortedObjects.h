/**
 *  @file   SortedObjects.h
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Base classes for chains and branched chains of nodes and segments.
 *          See PmaTrack3D.h file for details.
 */

#ifndef SortedObjects_h
#define SortedObjects_h

#include <vector>

namespace pma
{
	class SortedObjectBase;
	class SortedBranchBase;
}

class pma::SortedObjectBase
{
	friend class pma::SortedBranchBase;

public:
	SortedObjectBase(void) : next(0), prev(0) {}
	SortedObjectBase(pma::SortedObjectBase* prevElement, pma::SortedObjectBase* nextElement);

	/// Note: copy constructor does not preserve connections.
	SortedObjectBase(const pma::SortedObjectBase& src) : next(0), prev(0) {}

	virtual ~SortedObjectBase(void) { Disconnect(); }

	virtual void Disconnect(void);

	virtual bool AddNext(pma::SortedObjectBase* nextElement);
	virtual int RemoveNext(pma::SortedObjectBase* nextElement);

	virtual bool IsFirst(void) const { return !prev; }
	virtual bool IsLast(void) const { return !next; }

	virtual pma::SortedObjectBase* Prev(void) const { return prev; }
	virtual pma::SortedObjectBase* Next(unsigned int index = 0) const { return next; }
	virtual unsigned int NextCount(void) const
	{
		if (next) return 1;
		else return 0;
	}

protected:
	pma::SortedObjectBase* next;
	pma::SortedObjectBase* prev;
};


/// Base for classes, where a single object is assigned to Prev()
/// and many objects may be assigned to Next().
class pma::SortedBranchBase : public pma::SortedObjectBase
{
public:
	SortedBranchBase(void) : pma::SortedObjectBase() {}
	SortedBranchBase(pma::SortedObjectBase* prevElement, pma::SortedObjectBase* nextElement = 0) :
		pma::SortedObjectBase(prevElement, nextElement)
	{
		if (nextElement) next_vector.push_back(next);
	}
	/// Note: copy constructor does not preserve connections.
	SortedBranchBase(const pma::SortedBranchBase& src) : pma::SortedObjectBase() {}

	virtual ~SortedBranchBase(void) { Disconnect(); }

	virtual void Disconnect(void);

	virtual bool AddNext(pma::SortedObjectBase* nextElement);
	virtual int RemoveNext(pma::SortedObjectBase* nextElement);

	virtual pma::SortedObjectBase* Next(unsigned int index = 0) const
	{
		if (next_vector.size()) return next_vector[index];
		else return 0;
	}
	virtual unsigned int NextCount(void) const { return next_vector.size(); }
	virtual bool IsLast(void) const { return !(next_vector.size()); }

protected:
	std::vector< pma::SortedObjectBase* > next_vector;
};

#endif

