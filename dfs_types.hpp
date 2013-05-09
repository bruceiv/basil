#ifndef DFS_TYPES_HPP_
#define DFS_TYPES_HPP_

/** Data types used by the Basil DFS algorithm.
 *
 *  @author Aaron Moss
 */

/*  Copyright: Aaron Moss, 2012, moss.aaron@unb.ca  */

/*  This file is part of Basil.

    Basil is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Basil is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Basil.  If not, see <http://www.gnu.org/licenses/>.  */

#include <boost/unordered_map.hpp>

#include <gmpxx.h>

#include "basil.hpp"
#include "gram.hpp"

#include "lrs/cobasis.hpp"
#include "lrs/matrix.hpp"

namespace basil {

	////////////////////////////////////////////////////////////////////////
	// Typedefs for DFS
	////////////////////////////////////////////////////////////////////////

	typedef lrs::vector_mpq coordinates;

	typedef std::vector<index_set> index_set_list;

	/** Joint vertex-cobasis storage */
	struct vertex_data {

		/** Single-cobasis constructor. Initializes all fields as you would
		 *  think, where cobs is set up to be a set initially including
		 *  only cob.
		 */
		vertex_data(coordinates coords, index_set inc, index_set cob,
				mpz_class det, gram_matrix gram) : coords(coords),
				inc(inc), cobs(), det(det), gram(gram) {
			cobs.insert(cob);
		}

		/** Multiple-cobasis constructor. Initializes all fields to the
		 *  given values */
		vertex_data(coordinates coords, index_set inc,
				std::set<index_set> cobs, mpz_class det, gram_matrix gram)
				: coords(coords), inc(inc), cobs(cobs), det(det),
				gram(gram) { }

		/* Key data */
		/** Coordinates of the vertex */
		coordinates coords;
		/** Set of incident cobasis indices */
		index_set inc;
		/** Set of cobases for this vertex */
		std::set<index_set> cobs;

		/* Invariants */
		/** determinant */
		mpz_class det;
		/** gram matrix */
		gram_matrix gram;
	};
	typedef shared_ptr<vertex_data> vertex_data_ptr;
	typedef std::vector<vertex_data_ptr> vertex_data_list;

	/** map of vertex coordinates to a vertex data pointer */
	typedef
		boost::unordered_map<
			coordinates, vertex_data_ptr, lrs::vector_mpq_hash>
		coordinates_map;
	/** map of a cobasis to a vertex data pointer */
	typedef
		boost::unordered_map<index_set, vertex_data_ptr, lrs::index_set_hash>
		cobasis_map;
	/** map of a gram vector to its cobases, and the associated vertex
	 *  data */
	typedef
		boost::unordered_multimap<
			gram_matrix,
			std::pair<index_set, vertex_data_ptr>,
			gram_matrix_hash>
		cobasis_gram_map;
	/** map of a gram vector to its vertices */
	typedef
		boost::unordered_multimap<
			gram_matrix, vertex_data_ptr, gram_matrix_hash>
		vertex_gram_map;

} /* namespace basil */

#endif /* DFS_TYPES_HPP_ */
