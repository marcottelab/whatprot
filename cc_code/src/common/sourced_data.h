// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_COMMON_SOURCED_DATA_H
#define FLUOROSEQ_COMMON_SOURCED_DATA_H

#include <type_traits>

namespace fluoroseq {

template<typename V, typename S>
class SourcedData {
public:
    SourcedData(V value, S source) : value(value), source(source) {}

    ~SourcedData() {
        if (std::is_pointer(V)) {
            delete value;
        }
    }

    V value;  // if D is a pointer type then value is owned.
    S source;  // if S is a pointer type then source is NOT owned.
};

template<typename S>
class SourceSet {
public:
    SourceSet(int num_sources,
              S* sources) : num_sources(num_sources), sources(sources) {}

    ~SourceSet() {
        delete[] sources;
    }

    // The list of sources (S* sources) is owned, but if S is a pointer type the
    // individual sources are NOT owned (for example sources[0] is not owned).
    S* sources;
    int num_sources;
};

template<typename S>
class SourceWithCount {
public:
    SourceWithCount(S source, int count) : source(source), count(count) {}

    S source;  // if S is a pointer type then source is NOT owned.
    int count;
};

template<typename S>
class SourceCountMap {
public:
    SourceCountMap(int num_sources, SourceWithCount<S>** sources_with_counts)
            : num_sources(num_sources),
              sources_with_counts(sources_with_counts) {}
    
    ~SourceCountMap() {
        for (int i = 0; i < num_sources; i++) {
            delete sources_with_counts[i];
        }
        delete[] sources_with_counts;
    }

    // The list (sources_with_counts) is owned. Its elements are owned (for
    // example, souces_with_counts[0] is owned). However, its elements do not
    // own the type S data member if S is a pointer type (for example,
    // sources_with_counts[0].source would NOT be owned if S is a pointer type).
    SourceWithCount<S>** sources_with_counts;
    int num_sources;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_SOURCED_DATA_H
