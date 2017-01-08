//////////////////////////////////////////////////////////////////////////////
// \version 
//
// \brief Helper functions for MVAReader and MVAWriter wrappers
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_MVAWRAPPERBASE_H
#define ANAB_MVAWRAPPERBASE_H

#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/AnalysisBase/MVAOutput.h"

#include <typeinfo>
#include <functional>
#include <string>

namespace anab {

/// Helper functions for MVAReader and MVAWriter wrappers.
class MVAWrapperBase {
public:

protected:

    std::string getProductName(std::type_info const & ti) const;
    size_t getProductHash(std::type_info const & ti) const { return ti.hash_code(); }

    template <class T, size_t N>
    std::array<float, N> pAccumulate(
        std::vector< art::Ptr<T> > const & items,
        std::vector< FeatureVector<N> > const & outs) const;

    template <class T, size_t N>
    std::array<float, N> pAccumulate(
        std::vector< art::Ptr<T> > const & items, std::vector<float> const & weights,
        std::vector< FeatureVector<N> > const & outs) const;

    template <class T, size_t N>
    std::array<float, N> pAccumulate(
        std::vector< art::Ptr<T> > const & items, std::function<float (T const &)> fweight,
        std::vector< FeatureVector<N> > const & outs) const;

    template <class T, size_t N>
    std::array<float, N> pAccumulate(
        std::vector< art::Ptr<T> > const & items, std::function<float (art::Ptr<T> const &)> fweight,
        std::vector< FeatureVector<N> > const & outs) const;
};

} // namespace anab

//----------------------------------------------------------------------------
// MVAReader functions.
//
template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
    std::vector< art::Ptr<T> > const & items,
    std::vector< anab::FeatureVector<N> > const & outs) const
{
    std::array<double, N> acc;
    acc.fill(0);

	float pmin = 1.0e-6, pmax = 1.0 - pmin;
	float log_pmin = log(pmin), log_pmax = log(pmax);

	for (auto const & ptr : items)
	{
		auto const & vout = outs[ptr.key()];
		for (size_t i = 0; i < vout.size(); ++i)
		{
		    float v;
			if (vout[i] < pmin) v = log_pmin;
			else if (vout[i] > pmax) v = log_pmax;
			else v = log(vout[i]);

			acc[i] += v;
		}
	}

	if (!items.empty())
	{
		double totp = 0.0;
		for (size_t i = 0; i < N; ++i)
		{
			acc[i] = exp(acc[i] / items.size());
			totp += acc[i];
		}
		for (size_t i = 0; i < N; ++i)
		{
			acc[i] /= totp;
		}
	}
	else std::fill(acc.begin(), acc.end(), 1.0 / N);


    std::array<float, N> result;
    for (size_t i = 0; i < N; ++i) result[i] = acc[i];
    return result;
}
//----------------------------------------------------------------------------

template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
    std::vector< art::Ptr<T> > const & items, std::vector<float> const & weights,
    std::vector< anab::FeatureVector<N> > const & outs) const
{
    std::array<double, N> acc;
    acc.fill(0);

	float pmin = 1.0e-6, pmax = 1.0 - pmin;
	float log_pmin = log(pmin), log_pmax = log(pmax);
	double totw = 0.0;

	for (size_t k = 0; k < items.size(); ++k)
	{
	    auto const & ptr = items[k];
		float w = weights[k];

		if (w == 0) continue;

		auto const & vout = outs[ptr.key()];
		for (size_t i = 0; i < vout.size(); ++i)
		{
		    float v;
			if (vout[i] < pmin) v = log_pmin;
			else if (vout[i] > pmax) v = log_pmax;
			else v = log(vout[i]);

			acc[i] += w * v;
		}
		totw += w;
	}

	if (!items.empty())
	{
		double totp = 0.0;
		for (size_t i = 0; i < N; ++i)
		{
			acc[i] = exp(acc[i] / totw);
			totp += acc[i];
		}
		for (size_t i = 0; i < N; ++i)
		{
			acc[i] /= totp;
		}
	}
	else std::fill(acc.begin(), acc.end(), 1.0 / N);


    std::array<float, N> result;
    for (size_t i = 0; i < N; ++i) result[i] = acc[i];
    return result;
}
//----------------------------------------------------------------------------

template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
    std::vector< art::Ptr<T> > const & items, std::function<float (T const &)> fweight,
    std::vector< anab::FeatureVector<N> > const & outs) const
{
    std::array<double, N> acc;
    acc.fill(0);

	float pmin = 1.0e-6, pmax = 1.0 - pmin;
	float log_pmin = log(pmin), log_pmax = log(pmax);
	double totw = 0.0;

	for (size_t k = 0; k < items.size(); ++k)
	{
	    auto const & ptr = items[k];
		float w = fweight(*ptr);

		if (w == 0) continue;

		auto const & vout = outs[ptr.key()];
		for (size_t i = 0; i < vout.size(); ++i)
		{
		    float v;
			if (vout[i] < pmin) v = log_pmin;
			else if (vout[i] > pmax) v = log_pmax;
			else v = log(vout[i]);

			acc[i] += w * v;
		}
		totw += w;
	}

	if (!items.empty())
	{
		double totp = 0.0;
		for (size_t i = 0; i < N; ++i)
		{
			acc[i] = exp(acc[i] / totw);
			totp += acc[i];
		}
		for (size_t i = 0; i < N; ++i)
		{
			acc[i] /= totp;
		}
	}
	else std::fill(acc.begin(), acc.end(), 1.0 / N);


    std::array<float, N> result;
    for (size_t i = 0; i < N; ++i) result[i] = acc[i];
    return result;
}
//----------------------------------------------------------------------------

template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
    std::vector< art::Ptr<T> > const & items, std::function<float (art::Ptr<T> const &)> fweight,
    std::vector< anab::FeatureVector<N> > const & outs) const
{
    std::array<double, N> acc;
    acc.fill(0);

	float pmin = 1.0e-6, pmax = 1.0 - pmin;
	float log_pmin = log(pmin), log_pmax = log(pmax);
	double totw = 0.0;

	for (size_t k = 0; k < items.size(); ++k)
	{
	    auto const & ptr = items[k];
		float w = fweight(ptr);

		if (w == 0) continue;

		auto const & vout = outs[ptr.key()];
		for (size_t i = 0; i < vout.size(); ++i)
		{
		    float v;
			if (vout[i] < pmin) v = log_pmin;
			else if (vout[i] > pmax) v = log_pmax;
			else v = log(vout[i]);

			acc[i] += w * v;
		}
		totw += w;
	}

	if (!items.empty())
	{
		double totp = 0.0;
		for (size_t i = 0; i < N; ++i)
		{
			acc[i] = exp(acc[i] / totw);
			totp += acc[i];
		}
		for (size_t i = 0; i < N; ++i)
		{
			acc[i] /= totp;
		}
	}
	else std::fill(acc.begin(), acc.end(), 1.0 / N);


    std::array<float, N> result;
    for (size_t i = 0; i < N; ++i) result[i] = acc[i];
    return result;
}
//----------------------------------------------------------------------------

#endif //ANAB_MVAWRAPPERBASE

