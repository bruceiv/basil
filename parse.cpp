/** Implements Basil input file parsing from parse.hpp.
 *
 *  @author Aaron Moss
 */

#include <istream>
#include <ostream>
#include <vector>

#include <boost/make_shared.hpp>

#include "automorphism.hpp"
#include "basil.hpp"
#include "parse.hpp"
#include "permUtils.hpp"

#include "lrs/cobasis.hpp"

#include "permlib/permlib_api.h"


namespace basil {
	
	std::ostream& operator<< (std::ostream& o, parse_results const& p) {
		
// 		permutation_group_ptr g;
		
		std::ostream& (*endl)(std::ostream&) = std::endl;
		
		/* print name */
		if ( ! p.name.empty() ) o << p.name << endl;
		
		/* print representation */
		switch ( p.rep ) {
		case halfspace: o << "H-representation"; break;
		case vertex: o << "V-representation"; break;
		case arrangement: o << "A-representation"; break;
		}
		o << endl;
		
		/* print pre-matrix unparsed lines */
		for (uind i = 0; i < p.preLines.size(); ++i) {
			o << p.preLines[i] << endl;
		}
		
		/* print linearities */
		if ( p.l->any() ) {
			o << "linearity " << p.l->count();
			for (lrs::index_set_iter iter = lrs::begin(*p.l); 
					iter != lrs::end(*p.l); ++iter) {
				o << " " << *iter;
			}
			o << endl;
		}
		
		/* print constraint matrix */
		o << "begin" << endl;
		o << p.m->size() << " " << p.m->dim() << " rational" << endl;
		for (ind i = 0; i < p.m->size(); ++i) {
			for (ind j = 0; j < p.m->dim(); ++j) {
				o << " " << p.m->elem(i,j);
			}
			o << "\n";
		}
		o << "end" << endl;
		
		/* print post-matrix unparsed lines */
		for (uind i = 0; i < p.postLines.size(); ++i) {
			o << p.postLines[i] << endl;
		}
		
		/* print gram matrix */
		switch ( p.gs ) {
		case gram_omitted: /* do nothing */ break;
		case gram_inexact: o << "gram inexact" << endl; break;
		case gram_euclid: o << "gram Euclid" << endl; break;
		case gram_q: o << "gram Q" << endl; break;
		case gram_provided:
			o << "gram begin" << endl;
			for (uind i = 0; i < p.gm->dim(); ++i) {
				for (uind j = 0; j < p.gm->dim(); ++j) {
					o << " " << (*p.gm)(i, j);
				}
				o << endl;
			}
			o << "gram end" << endl;
			break;
		}
		
		/* print symmetry group */
		switch ( p.ss ) {
		case sym_omitted: /* do nothing */ break;
		case sym_auto: o << "symmetry auto" << endl; break;
		case sym_provided:
			o << "symmetry begin" << endl;
			permutation_list gmgs = small_gen_set(*p.g);
			typedef permutation_list::iterator perm_iter;
			for (perm_iter iter = gmgs.begin(); iter != gmgs.end(); ++iter) {
				o << in_str(**iter) << endl;
			}
			o << "symmetry end" << endl;
			break;
		}
		
		return o;
	}
	
	/** Gets a line that isn't a comment (begins with '*' or '#').
	 *  @param in		the input stream to read the line from
	 *  @param s		the string to store the line in
	 *  @param v		optional vector to push comment lines back to
	 *  @return the input stream
	 */
	std::istream& getContentLine(std::istream& in, string& s, 
								 std::vector<string>* v = 0) {
		std::getline(in, s);
		while ( in && (s.empty() || s[0] == '*' || s[0] == '#') ) {
			if ( v ) v->push_back(s);
			std::getline(in, s);
		}
		if ( ! in ) s = "";
		
		return in;
	}
	
	/** Checks if a string begins with a given prefix
	 *  @param s		the string to search
	 *  @param pre		the prefix to match
	 *  @return true if the prefix matches, false otherwise
	 */
	bool prefixMatch(std::string const& s, char const* pre) {
		for (unsigned int i = 0; pre[i]; ++i) {
			if ( 
				i >= s.length() 	/* s is shorter than prefix */
				|| s[i] != pre[i]	/* mismatch with prefix */
			) return false;
		}
		return true;				/* match with prefix in all locations */
	}
	
	parse_results_ptr parse(std::istream& in) {
		parse_results_ptr p = boost::make_shared<parse_results>();
		
		string s = ""; getContentLine(in, s, &p->preLines);
		/* temporary linearity vector */
		std::vector<ind> linV(0);
		
		bool firstLoop = true;
		/* parse options up to begin line */
		while ( ! prefixMatch(s, "begin") ) {
			
			if ( prefixMatch(s, "H-representation") ) {
				/* Set halfspace representation flag */
				p->rep = halfspace;
			} else if ( prefixMatch(s, "V-representation") ) {
				/* Set vertex representation flag */
				p->rep = vertex;
			} else if ( prefixMatch(s, "A-representation") ) {
				/* Set arrangement representation flag */
				p->rep = arrangement;
			} else if ( prefixMatch(s, "linearity") ) {
				/* parse linearities */
				std::istringstream read(s);
				ind k, t;
				read >> k; /* linearity count */
				linV.resize(k);
				/* read linearities */
				for (ind i = 0; i < k; ++i) { read >> linV[i]; }
			} else if (firstLoop) {
				p->name = s;
			} else {
				p->preLines.push_back(s);
			}
			
			/* get next line */
			firstLoop = false;
			getContentLine(in, s, &p->preLines);
		}
		
		p->m = parseMatrix(in);
		
		/* read linearities into index set */
		p->l = boost::make_shared<index_set>(p->m->size()+1);
		for (std::vector<ind>::iterator iter = linV.begin(); 
				iter != linV.end(); ++iter) p->l->set(*iter);
		
		/* continue parsing */
		while ( in ) {
			getContentLine(in, s, &p->postLines);
			
			if ( prefixMatch(s, "symmetry") ) {
				if ( prefixMatch(s, "symmetry auto") ) {
					p->ss = sym_auto;
				} else if ( prefixMatch(s, "symmetry begin") ) {
					p->ss = sym_provided;
					p->g = parsePermutationGroup(in, p->m->size());
				} else goto pushLine;
			} else if ( prefixMatch(s, "gram") ) {
				if ( prefixMatch(s, "gram auto") ) {
					p->gs = gram_euclid;
				} else if ( prefixMatch(s, "gram Euclid") ) {
					p->gs = gram_euclid;
				} else if ( prefixMatch(s, "gram inexact") ) {
					p->gs = gram_inexact;
				} else if ( prefixMatch(s, "gram Q") ) {
					p->gs = gram_q;
				} else if ( prefixMatch(s, "gram begin") ) {
					p->gs = gram_provided;
					p->gm = parseGram(in, p->m->size());
				} else goto pushLine;
			} else {
				if ( in ) 
					pushLine: p->postLines.push_back(s);
			}
		}
		
		return p;
	}
	
	matrix_ptr parseMatrix(std::istream& in) {
		/* get dimension line */
		string s = "";
		getContentLine(in, s);
		std::istringstream read(s);
		
		/* read dimensions */
		ind n, d;
		read >> n;
		read >> d;
		
		/* create new matrix and load data */
		matrix_ptr m = boost::make_shared<matrix>(n, d);
		for (ind i = 0; i < n; i++) {
			for (ind j = 0; j < d; j++) {
				in >> m->elem(i,j);
				m->elem(i,j).canonicalize();
			}
		}
		
		/* ignore up to end line */
		getContentLine(in, s);
		while ( s != string("end") ) getContentLine(in, s);
		
		return m;
	}
	
	permutation_group_ptr parsePermutationGroup(std::istream& in, ind n) {
		std::vector<permutation_ptr> generators;
		
		//read in generators
		string s = "";
		permutation_ptr p;
		
		std::getline(in, s);
		while ( s != string("symmetry end") ) {
			p.reset(new permutation(n, s));
			generators.push_back(p);
			
			std::getline(in, s);
		}
		
		return permlib::construct(n, generators.begin(), generators.end());
	}
	
	gram_matrix_ptr parseGram(std::istream& in, ind n) {
		/* create new matrix and load data */
		gram_matrix_ptr gm = boost::make_shared<gram_matrix>(n);
		for (ind i = 0; i < n; ++i) for (ind j = 0; j < n; ++j) {
			in >> gm->at(i,j);
			if ( gm->at(i,j) >= (int)gm->k() ) gm->k() = gm->at(i,j) + 1;
		}
		
		/* ignore up to end line */
		string s;
		getContentLine(in, s);
		while ( s != string("gram end") ) getContentLine(in, s);
		
		return gm;
	}
	
} /* namespace basil */
