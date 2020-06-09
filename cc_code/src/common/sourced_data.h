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
        if (std::is_pointer<V>::value) {
            delete value;
        }
        if (std::is_pointer<S>::value) {
            delete source;
        }
    }

    V value;  // if D is a pointer type then value is owned.
    S source;  // if S is a pointer type then source is owned.
};

template<typename S>
class SourceSet {
public:
    SourceSet(int num_sources,
              S* sources) : num_sources(num_sources), sources(sources) {}

    ~SourceSet() {
        if (std::is_pointer<S>::value) {
            for (int i = 0; i < num_sources; i++) {
                delete sources[i];
            }
        }
        delete[] sources;
    }

    S* sources;  // owned.
    int num_sources;
};

template<typename S>
class SourceWithCount {
public:
    SourceWithCount(S source, int count) : source(source), count(count) {}

    ~SourceWithCount() {
        if (std::is_pointer<S>::value) {
            delete source;
        }
    }

    S source;  // if S is a pointer type then source is owned.
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

    SourceWithCount<S>** sources_with_counts;  // owned.
    int num_sources;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_SOURCED_DATA_H
