/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_UTIL_VECTOR_HASH_H
#define WHATPROT_UTIL_VECTOR_HASH_H

// Standard C++ library headers:
#include <functional>
#include <vector>

// An std namespace injection is the accepted way of creating a new hash
// function.
namespace std {

// The hash formula used here is taken from the Boost library. Replicating it
// here avoids a Boost dependency. More information can be found in the
// following links:
//   -
//   https://www.boost.org/doc/libs/1_74_0/doc/html/hash/reference.html#boost.hash_combine
//   -
//   https://stackoverflow.com/questions/4948780/magic-number-in-boosthash-combine
template <typename K>
struct hash<vector<K>> {
public:
    size_t operator()(const vector<K>& v) const {
        hash<K> hash_k;
        size_t seed = 0;
        for (const K& k : v) {
            seed ^= hash_k(k) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

}  // namespace std

#endif  // WHATPROT_UTIL_VECTOR_HASH_H
