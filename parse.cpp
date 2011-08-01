#include <istream>
#include <ostream>
#include <vector>

#include <boost/make_shared.hpp>

#include <permlib/permlib_api.h>

#include "basil.hpp"
#include "parse.hpp"

namespace basil {
	
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
	
	parse_results parse(std::istream& in) {
		parse_results p;
		
		string s = ""; getContentLine(in, s, &p.preLines);
		/* temporary linearity vector */
		std::vector<ind> linV(0);
		
		/* parse options up to begin line */
		while ( ! prefixMatch(s, "begin") ) {
			
			if ( prefixMatch(s, "V-representation") ) {
				/* Set vertex representation flag */
				p.rep = vertex;
			} else if ( prefixMatch(s, "A-representation") ) {
				/* Set arrangement representation flag */
				p.rep = arrangement;
			} else if ( prefixMatch(s, "linearity") ) {
				/* parse linearities */
				std::istringstream read(s);
				ind k, t;
				read >> k; /* linearity count */
				linV.resize(k);
				/* read linearities */
				for (ind i = 0; i < k; ++i) { read >> linV[i]; }
			} else {
				p.preLines.push_back(s);
			}
			
			/* get next line */
			getContentLine(in, s, &p.preLines);
		}
		
		p.m = parseMatrix(in);
		
		/* read linearities into index set */
		p.l = boost::make_shared<index_set>(p.m->size()+1);
		for (std::vector<ind>::iterator iter = linV.begin(); 
				iter != linV.end(); ++iter) p.l->set(*iter);
		
		/* continue parsing */
		while ( in ) {
			getContentLine(in, s, &p.postLines);
			
			if ( prefixMatch(s, "symmetry begin") ) {
				/* parse permutation group, starting here */
				p.g = parsePermutationGroup(in, p.m->size());
			} else if ( prefixMatch(s, "gram") ) {
				if ( prefixMatch(s, "gram auto") ) {
					p.gs = exact;
				} else if ( prefixMatch(s, "gram inexact") ) {
					p.gs = inexact;
				} else if ( prefixMatch(s, "gram begin") ) {
					p.gs = provided;
					p.gm = parseGram(in, p.m->size());
				}
			} else {
				p.postLines.push_back(s);
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
