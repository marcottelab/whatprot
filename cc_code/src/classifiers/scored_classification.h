// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_CLASSIFIERS_SCORED_CLASSIFICATION_H
#define FLUOROSEQ_CLASSIFIERS_SCORED_CLASSIFICATION_H

#include "common/dye_seq.h"

namespace fluoroseq {

// ScoredClassification MUST remain compatible with offsetof in the <cstddef>
// header. This is NOT compatible with C++98, as C++98 would require
// ScoredClassification to be a POD (Plain Old Datatype), which is a much
// stricter requirement. This IS compatible with C++11. In C++11 this class
// is required to be a "standard layout" class. This gives the following
// requirements:
//     * no virtual functions and no virtual base classes.
//     * has the same access control (private, protected, public) for all its
//       non-static data members.
//     * either has no non-static data members in the most derived class and at
//       most one base class with non-static data members, or has no base
//       classes with non-static data members.
//     * its base class (if any) is itself also a standard-layout class. And,
//     * has no base classes of the same type as its first non-static data
//       member.
class ScoredClassification {
public:
    ScoredClassification(int id, double score, double total);
    ScoredClassification();
    double adjusted_score();

    double score;
    double total;
    int id;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_SCORED_CLASSIFICATION_H
