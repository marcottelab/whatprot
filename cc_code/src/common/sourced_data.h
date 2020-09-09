// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_COMMON_SOURCED_DATA_H
#define FLUOROSEQ_COMMON_SOURCED_DATA_H

#include "util/delete.h"

namespace fluoroseq {

template<typename V, typename S>
class SourcedData {
public:
    SourcedData(V value, S source) : value(value), source(source) {}

    ~SourcedData() {
        delete_if_pointer(value);
        delete_if_pointer(source);
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
        delete_array(num_sources, sources);
    }

    S* sources;  // owned.
    int num_sources;
};

template<typename S>
class SourceCount {
public:
    SourceCount(S source, int count) : source(source), count(count) {}
    
    SourceCount(const SourceCount& other) : source(other.source),
                                            count(other.count) {}

    ~SourceCount() {
        delete_if_pointer(source);
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
        delete_array(num_sources, sources);
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
        delete_if_pointer(source);
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

    SourceCountHitsList(const SourceCountHitsList& other)
            : num_sources(other.num_sources) {
        sources = new SourceCountHits<S>*[num_sources];
        for (int i = 0; i < num_sources; i++) {
            sources[i] = new SourceCountHits<S>(other.sources[i]->source,
                                                other.sources[i]->count,
                                                other.sources[i]->hits);
        }
    }

    SourceCountHitsList(SourceCountHitsList&& other) 
            : num_sources(other.num_sources), sources(other.sources) {
        other.sources = NULL;
    }
    
    ~SourceCountHitsList() {
        if (sources != NULL) {
            delete_array(num_sources, sources);
        }
    }

    SourceCountHits<S>** sources;
    int num_sources;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_SOURCED_DATA_H
