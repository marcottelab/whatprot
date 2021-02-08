/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TEST_UTIL_FAKEIT_H
#define WHATPROT_TEST_UTIL_FAKEIT_H

// Standard C++ library headers:
#include <string>
#include <type_traits>

// External headers:
#include "fakeit.hpp"

namespace whatprot {
namespace test_util {

template <typename T>
class CloseMatcherCreator : public fakeit::TypedMatcherCreator<T> {
public:
    CloseMatcherCreator(const T& expected, const T& tolerance)
            : expected(expected), tolerance(tolerance) {}

    virtual ~CloseMatcherCreator() = default;

    class Matcher : public fakeit::TypedMatcher<T> {
    public:
        Matcher(const T& expected, const T& tolerance)
                : expected(expected), tolerance(tolerance) {}

        virtual bool matches(const T& actual) const override {
            return (actual >= expected - tolerance)
                   && (actual <= expected + tolerance);
        }

        virtual std::string format() const override {
            return std::string("Approximately ")
                   + fakeit::TypeFormatter<T>::format(expected);
        }

        T expected;
        T tolerance;
    };

    virtual fakeit::TypedMatcher<T>* createMatcher() const override {
        return new Matcher(expected, tolerance);
    }

    T expected;
    T tolerance;
};

template <typename T>
CloseMatcherCreator<T> Close(const T& expected, const T& tolerance) {
    CloseMatcherCreator<T> close(expected, tolerance);
    return close;
}

// Type P should ALWAYS be T*. This is an ugly hack to get the matches()
// function to be recognized as an override of the matches() function in
// fakeit::TypedMatcher<T*>. The compiler doesn't seem to recognize that
// "const T*&" with "T = SomeType" is the same as "const T&" with
// "T = SomeType*". Although this is ugly, it's hidden from the user, who just
// calls the Ptr() function below, so it's probably OK.
template <typename T, typename P>
class PtrMatcherCreator : public fakeit::TypedMatcherCreator<P> {
public:
    PtrMatcherCreator(fakeit::TypedMatcherCreator<T>* matcher_creator)
            : matcher_creator(matcher_creator) {}

    virtual ~PtrMatcherCreator() {
        delete matcher_creator;
    }

    class Matcher : public fakeit::TypedMatcher<P> {
    public:
        Matcher(fakeit::TypedMatcher<T>* matcher) : matcher(matcher) {}

        virtual ~Matcher() {
            delete matcher;
        }

        virtual bool matches(const P& actual) const override {
            return matcher->matches(*actual);
        }

        virtual std::string format() const override {
            return std::string("Pointer ") + matcher->format();
        }

        typename fakeit::TypedMatcher<T>* matcher;
    };

    virtual fakeit::TypedMatcher<P>* createMatcher() const override {
        return new Matcher(matcher_creator->createMatcher());
    }

    fakeit::TypedMatcherCreator<T>* matcher_creator;
};

// Here we check whether we received a fakeit::TypedMatcherCreator as the
// argument. If we did, we enable this definition. The return type of the
// function when it's enabled is PtrMatcherCreator<T, T*>.
template <typename T, template <typename> class MC>
typename std::enable_if<
        std::is_base_of<fakeit::TypedMatcherCreator<T>, MC<T>>::value,
        PtrMatcherCreator<T, T*>>::type
Ptr(MC<T> matcher_creator) {
    // We need to copy into the heap so that we can access it from
    // PtrMatcherCreator and use polymorphism. This is also an additional reason
    // to template MC, the "Matcher Creator." Otherwise we wouldn't know what
    // type to construct.
    MC<T>* new_matcher_creator = new MC<T>(matcher_creator);
    PtrMatcherCreator<T, T*> ptr(new_matcher_creator);
    return ptr;
}

// If argument is given directly, instead of through a TypedMatcherCreator, as
// above, then we infer that an EqMatcherCreator for the given type is wanted.
template <typename T>
PtrMatcherCreator<T, T*> Ptr(const T& expected) {
    PtrMatcherCreator<T, T*> ptr(new fakeit::EqMatcherCreator<T>(expected));
    return ptr;
}

}  // namespace test_util
}  // namespace whatprot

#endif  // WHATPROT_TEST_UTIL_FAKEIT_H
