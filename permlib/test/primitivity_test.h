// ---------------------------------------------------------------------------
//
//  This file is part of PermLib.
//
// Copyright (c) 2009-2011 Thomas Rehn <thomas@carmen76.de>
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//    derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ---------------------------------------------------------------------------


#ifndef PRIMITIVITY_TEST_H_
#define PRIMITIVITY_TEST_H_

#include <permlib/prime_helper.h>

#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <vector>
#include <list>

namespace permlib {

/// Tests a transitive group is availble for primitivity.
/**
 * If group is not primitive, it can compute a minimal block.
 * Note that PrimitivitySGSTest may be faster if a strong generating set is known.
 * 
 * This class implements the algorithm described in
 * Holt, Eick, O'Brien: Handbook of Computational Group Theory, 2005. Chapter 4.3
 */
template<typename PERM>
class PrimitivityTest {
public:
	/**
	 * Sets up the test
	 * 
	 * @param n number of elements in the group domain
	 * @param genBegin iterator<PERM::ptr> begin for group generators
	 * @param genEnd iterator<PERM::ptr> end for group generators
	 */
	template<typename InputIterator>
	PrimitivityTest(const unsigned int n, InputIterator genBegin, InputIterator genEnd);
	
	/**
	 * @param minimalBlock If not null, this vector will be filled a (non-sorted) minimal block for the group.
	 *                     If the group is primitive, the vector will contain all elements of the domain.
	 * @return true iff group is primitive
	 */
	bool blockOfImprimitivity(std::vector<dom_int>* minimalBlock) const;
	
	/**
	 * @return true iff group is primitive
	 */
	bool isPrimitive() const { return blockOfImprimitivity(NULL); }
	
private:
	const unsigned int m_n;
	unsigned int m_primeLimit;
	const std::list<typename PERM::ptr> m_generators;
	
	bool fillTrivialBlock(std::vector<dom_int>* minimalBlock) const;
	
	static dom_int rep(dom_int kappa, std::vector<dom_int>& p);
	
	bool merge(dom_int kappa, dom_int lambda, std::vector<dom_int>& c, std::vector<dom_int>& p, std::vector<dom_int>& q, unsigned int& l) const;
};



template<typename PERM>
template<typename InputIterator>
PrimitivityTest<PERM>::PrimitivityTest(const unsigned int n, InputIterator genBegin, InputIterator genEnd)
	: m_n(n), m_primeLimit(m_n), m_generators(genBegin, genEnd)
{
	for (const unsigned int* p = PrimeHelper::firstPrime(); p != PrimeHelper::lastPrime(); ++p) {
		if (m_n % (*p) == 0) {
			m_primeLimit = m_n / (*p);
			break;
		}
	}
}


template<typename PERM>
bool PrimitivityTest<PERM>::blockOfImprimitivity(std::vector<dom_int>* minimalBlock) const {
	std::vector<dom_int> alphas(2);
	alphas[0] = 0;
	
	for (dom_int a = 1; a < m_n; ++a) {
		alphas[1] = a;
		
		const unsigned int k = alphas.size();
		unsigned int l = k - 1;
		std::vector<dom_int> p(m_n);
		std::vector<dom_int> q(m_n);
		std::vector<dom_int> c(m_n);
		
		for (unsigned int i = 0; i < m_n; ++i) {
			c[i] = 1;
			p[i] = i;
		}
		
		for (unsigned int i = 0; i < k - 1; ++i) {
			p[alphas[i+1]] = alphas[0];
			q[i] = alphas[i+1];
		}
		
		bool tryNextAlpha = false;
		c[alphas[0]] = k;
		for (unsigned int i = 0; i < l; ++i) {
			const dom_int gamma = q[i];
			BOOST_FOREACH(const typename PERM::ptr& x, m_generators) {
				const dom_int delta = rep(gamma, p);
				if (merge(x->at(gamma), x->at(delta), c, p, q, l)) {
					tryNextAlpha = true;
					goto TRY_NEXT_ALPHA;
				}
			}
		}
TRY_NEXT_ALPHA:
		if (tryNextAlpha)
			continue;
		
		for (unsigned int i = 0; i < m_n; ++i)
			rep(i, p);
		
		const unsigned int minBlockSize = c[rep(alphas[0], p)];
		if (minBlockSize < m_n) {
			if (minimalBlock) {
				minimalBlock->clear();
				minimalBlock->reserve(minBlockSize);
				for (unsigned int i = 0; i < m_n; ++i)
					if (p[i] == p[alphas[0]])
						minimalBlock->push_back(i);
			}
			return false;
		}
	}
	
	return fillTrivialBlock(minimalBlock);
}

template<typename PERM>
bool PrimitivityTest<PERM>::fillTrivialBlock(std::vector<dom_int>* minimalBlock) const {
	if (minimalBlock) {
		minimalBlock->clear();
		minimalBlock->resize(m_n);
		for (unsigned int i = 0; i < m_n; ++i)
			minimalBlock->at(i) = i;
	}
	return true;
}

template<typename PERM>
dom_int PrimitivityTest<PERM>::rep(dom_int kappa, std::vector<dom_int>& p) {
	dom_int lambda = kappa;
	dom_int rho = p[lambda];
	while (rho != lambda) {
		lambda = rho;
		rho = p[lambda];
	}
	
	dom_int mu = kappa;
	rho = p[mu];
	while (rho != lambda) {
		p[mu] = lambda;
		mu = rho;
		rho = p[mu];
	}
	
	return lambda;
}

template<typename PERM>
bool PrimitivityTest<PERM>::merge(dom_int kappa, dom_int lambda, std::vector<dom_int>& c, std::vector<dom_int>& p, std::vector<dom_int>& q, unsigned int& l) const {
	dom_int phi = rep(kappa, p);
	dom_int psi = rep(lambda, p);
	
	if (phi != psi) {
		dom_int mu, nu;
		if (c[phi] >= c[psi]) {
			mu = phi;
			nu = psi;
		} else {
			mu = psi;
			nu = phi;
		}
		p[nu] = mu;
		c[mu] += c[nu];
		if (c[mu] > m_primeLimit)
			return true;
		
		q[l] = nu;
		++l;
	}
	
	return false;
}

}

#endif
