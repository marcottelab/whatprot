// Author: Matthew Beauregard Smith
#include <boost/test/unit_test.hpp>

#include "common/error_model.h"

#include <functional>

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
using std::function;
const double TOL = 0.000000001;
}

BOOST_AUTO_TEST_SUITE(common_suite);
BOOST_AUTO_TEST_SUITE(error_model_suite);

BOOST_AUTO_TEST_CASE(constructor_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = 1.0;
    double sigma = .16;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma);
    BOOST_TEST(em.p_edman_failure == p_edman_failure);
    BOOST_TEST(em.p_detach == p_detach);
    BOOST_TEST(em.p_bleach == p_bleach);
    BOOST_TEST(em.p_dud == p_dud);
    BOOST_TEST(em.distribution_type == dist_type);
    BOOST_TEST(em.mu == mu);
    BOOST_TEST(em.sigma == sigma);
}

BOOST_AUTO_TEST_CASE(pdf_state_zero_obs_zero_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = 1.0;
    double sigma = .16;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma);
    function<double (double, int)> pdf = em.pdf();
    double observed = 0.0;
    int state = 0;
    BOOST_TEST(pdf(observed, state) == 1.0);
}

BOOST_AUTO_TEST_CASE(pdf_state_zero_obs_one_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = 1.0;
    double sigma = .16;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma);
    function<double (double, int)> pdf = em.pdf();
    double observed = 1.0;
    int state = 0;
    BOOST_TEST(pdf(observed, state) == 0.0);
}

BOOST_AUTO_TEST_CASE(pdf_state_one_obs_zero_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = 1.0;
    double sigma = .16;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma);
    function<double (double, int)> pdf = em.pdf();
    double observed = 0.0;
    int state = 1;
    BOOST_TEST(pdf(observed, state) == 0.0);
}

BOOST_AUTO_TEST_CASE(pdf_state_one_obs_one_test, * tolerance(TOL)) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = 1.0;
    double sigma = .16;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma);
    function<double (double, int)> pdf = em.pdf();
    double observed = 1.0;
    int state = 1;
    // The test value was found using an online lognormal distribution pdf
    // calculator.
    BOOST_TEST(pdf(observed, state) == 2.4933892525089547);
}

BOOST_AUTO_TEST_SUITE_END();  // error_model_suite
BOOST_AUTO_TEST_SUITE_END();  // common_suite

}  // namespace fluoroseq
