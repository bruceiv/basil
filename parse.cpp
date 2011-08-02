#include <istream>
#include <ostream>
#include <vector>
#include <set>
#include <sstream>

#include <boost/make_shared.hpp>

#include <permlib/permlib_api.h>

#include "lrs/cobasis.hpp"

#include "basil.hpp"
#include "parse.hpp"

namespace basil {
	
	/** Converts a permutation to a string compatible with its input format. 
	 *  Derived from standard printing code in PermLib */
	string in_str(permutation const& p) {
		typedef permlib::dom_int dom_int;
		std::ostringstream o;
		
		std::set<dom_int> worked;
		for (dom_int x = 0; x < p.size(); ++x) {
			dom_int px, startX;
			startX = x;
			px = p/x;
			if (worked.count(x) || x == px) {
				continue;
			}

			if ( ! worked.empty() ) o << " ,";
			worked.insert(x);
			o << " " << (x+1);
			while (p/px != startX) {
				o << " " << (px+1);
				worked.insert(px);
				px = p/px;
			}
			worked.insert(px);
			o << " " << (px+1);
		}
		return o.str();
	}
	
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
		switch (p.gs) {
		case ommited: /* do nothing */ break;
		case inexact: o << "gram inexact" << endl; break;
		case exact: o << "gram auto" << endl; break;
		case provided:
			o << "gram begin" << endl;
			for (ind i = 0; i < p.gm->dim(); ++i) {
				for (ind j = 0; j < p.gm->dim(); ++j) {
					o << " " << (*p.gm)(i, j);
				}
				o << endl;
			}
			o << "gram end" << endl;
			break;
		}
		
		if ( p.g ) {
			o << "symmetry begin" << endl;
			/* gsgs is the strong generating set of *p.g, a list of pointers to 
			 * permutations */
			permutation_group::PERMlist& gsgs = p.g->S;
			typedef permutation_group::PERMlist::iterator perm_iter;
			for (perm_iter iter = gsgs.begin(); iter != gsgs.end(); ++iter) {
				o << in_str(**iter) << endl;
			}
			o << "symmetry end" << endl;
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
			
			if ( prefixMatch(s, "symmetry begin") ) {
				/* parse permutation group, starting here */
				p->g = parsePermutationGroup(in, p->m->size());
			} else if ( prefixMatch(s, "gram") ) {
				if ( prefixMatch(s, "gram auto") ) {
					p->gs = exact;
				} else if ( prefixMatch(s, "gram inexact") ) {
					p->gs = inexact;
				} else if ( prefixMatch(s, "gram begin") ) {
					p->gs = provided;
					p->gm = parseGram(in, p->m->size());
				}
			} else {
				if ( in ) p->postLines.push_back(s);
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
			in >> (*gm)(i,j);
		}
		
		/* ignore up to end line */
		string s;
		getContentLine(in, s);
		while ( s != string("gram end") ) getContentLine(in, s);
		
		return gm;
	}
	
} /* namespace basil */
