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


#ifndef DADE_INVARIANTS_H_
#define DADE_INVARIANTS_H_

#include <permlib/transversal/orbit_set.h>
#include <permlib/invariant/linear_form_list.h>

namespace permlib {

/// computes some invariants of a permutation group with Dade's algorithm
template<class BSGSIN>
class DadeInvariants {
public:
	typedef typename BSGSIN::PERMtype PERM;
	typedef typename BSGSIN::TRANStype TRANS;

	/// constructor
	DadeInvariants(const BSGSIN& bsgs);

	/// destructor
	virtual ~DadeInvariants(){}

	/// computes some algebraically independent invariants, but not necessarily all generators of the invariant ring
	/**
	 * @param invariantList list to store the invariants (LinearFormList)
	 * @param maximalDegree maximal degree that a constructed invariant may have, or zero if unbounded
	 */
	void invariants(std::list<LinearFormList>& invariantList, unsigned int maximalDegree = 0) const;
private:
	const BSGSIN& m_bsgs;
};

//
// IMPLEMENTATION
//


template<class BSGSIN>
DadeInvariants<BSGSIN>::DadeInvariants(const BSGSIN& bsgs) 
	: m_bsgs(bsgs)
{  }

template<class BSGSIN>
void DadeInvariants<BSGSIN>::invariants(std::list<LinearFormList>& invariantList, unsigned int maximalDegree) const
{
	const unsigned long n = m_bsgs.n;
	
	std::map<unsigned long,OrbitSet<PERM,unsigned long> > orbits;
	boost::dynamic_bitset<> checked(n);
	unsigned long alpha = 0;
	while (true) {
		while (checked[alpha] && alpha < n) {
			++alpha;
		}
		if (alpha >= n)
			break;
		
		OrbitSet<PERM, unsigned long> orbit;
		orbit.orbit(alpha, m_bsgs.S, typename Transversal<PERM>::TrivialAction());
		BOOST_FOREACH(const unsigned long& beta, std::make_pair(orbit.begin(), orbit.end())) {
			checked.set(beta, 1);
		}
		orbits.insert(std::make_pair(alpha, orbit));
	}
	
	typedef std::pair<unsigned long,OrbitSet<PERM,unsigned long> > pair_t;
	std::set<unsigned long>::const_iterator setIt;
	BOOST_FOREACH(const pair_t& orbit, orbits) {
		unsigned long l = 1;
		while (l < orbit.second.size()) {
			LinearForm form(n);
			setIt = orbit.second.begin();
			for (unsigned int i = 0; i < orbit.second.size()-l; ++i) {
				form.set(*setIt, 1);
				++setIt;
			}
			
			OrbitSet<PERM, LinearForm> formOrbit;
			formOrbit.orbit(form, m_bsgs.S, LinearFormAction<PERM>());
			LinearFormList list;
			BOOST_FOREACH(const LinearForm& lform, std::make_pair(formOrbit.begin(), formOrbit.end())) {
				list.add(lform);
			}
			
			if (!maximalDegree || list.size() < maximalDegree)
				invariantList.push_back(list);
			l <<= 1;
		}
	}
}

}

#endif // -- DADE_INVARIANTS_H_
