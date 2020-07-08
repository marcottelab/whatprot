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
class SourceCount {
public:
    SourceCount(S source, int count) : source(source), count(count) {}

    ~SourceCount() {
        if (std::is_pointer<S>::value) {
            delete source;
        }
    }

    S source;  // if S is a pointer type then source is owned.
    int count;  // A count of the number of sources.
};

template<typename S>
class SourceCountList {
public:
    SourceCountList(int num_sources, SourceCount<S>** sources)
            : num_sources(num_sources), sources(sources) {}
    
    ~SourceCountList() {
        for (int i = 0; i < num_sources; i++) {
            delete sources[i];
        }
        delete[] sources;
    }

    SourceCount<S>** sources;  // owned.
    int num_sources;
};

template<typename S>
class SourceCountHits {
public:
    SourceCountHits(S source, int count, int hits)
            : source(source), count(count), hits(hits) {}

    ~SourceCountHits() {
        if (std::is_pointer<S>::value) {
            delete source;
        }
    }

    S source;  // if S is a pointer type then source is owned.
    int count;  // A count of the number of sources.
    int hits;  // A count of the number of hits in a simulation.
};

template<typename S>
class SourceCountHitsList {
public:
    SourceCountHitsList(int num_sources, SourceCountHits<S>** sources)
            : num_sources(num_sources), sources(sources) {}
    
    ~SourceCountHitsList() {
        for (int i = 0; i < num_sources; i++) {
            delete sources[i];
        }
        delete[] sources;
    }

    SourceCountHits<S>** sources;
    int num_sources;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_SOURCED_DATA_H
