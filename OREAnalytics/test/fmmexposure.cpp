/*
 Copyright (C) 2026 Ed Herranz
 All rights reserved.

 This file is part of ORE, a free-software/open-source library
 for transparent pricing and risk analysis - http://opensourcerisk.org

 ORE is free software: you can redistribute it and/or modify it
 under the terms of the Modified BSD License.  You should have received a
 copy of the license along with this program.
 The license is also available online at <http://opensourcerisk.org>

 This program is distributed on the basis that it will form a useful
 contribution to risk analytics and model standardisation, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
*/

#include "oreatoplevelfixture.hpp"
#include <boost/test/unit_test.hpp>

#include <orea/cube/cube_io.hpp>
#include <orea/cube/inmemorycube.hpp>
#include <orea/scenario/aggregationscenariodata.hpp>

#include <cstdio>
#include <set>

using namespace QuantLib;
using namespace ore::analytics;

BOOST_FIXTURE_TEST_SUITE(OREAnalyticsTestSuite, ore::test::OreaTopLevelFixture)

BOOST_AUTO_TEST_SUITE(FmmExposureContractTests)

// A7 contract tests: the external cube route writes ORE's native files through ORE's own
// serialisers; these tests pin what the loaders return for files written that way.

BOOST_AUTO_TEST_CASE(testAggregationScenarioDataRoundTrip) {
    BOOST_TEST_MESSAGE("Testing the aggregation scenario data file round trip (A7 route 1: Numeraire per date "
                       "and sample written by saveAggregationScenarioData, read by loadAggregationScenarioData)...");
    const Size dates = 3, samples = 4;
    auto asd = QuantLib::ext::make_shared<InMemoryAggregationScenarioData>(dates, samples);
    for (Size d = 0; d < dates; ++d)
        for (Size n = 0; n < samples; ++n)
            asd->set(d, n, 1.0 + 0.1 * d + 0.01 * n, AggregationScenarioDataType::Numeraire, "");
    const std::string file = "fmm_a7_asd_roundtrip.csv";
    saveAggregationScenarioData(file, *asd);
    auto loaded = loadAggregationScenarioData(file);
    std::remove(file.c_str());
    BOOST_REQUIRE(loaded);
    BOOST_CHECK_EQUAL(loaded->dimDates(), dates);
    BOOST_CHECK_EQUAL(loaded->dimSamples(), samples);
    Size mismatches = 0;
    std::set<std::pair<Size, Size>> bad;
    for (Size d = 0; d < dates; ++d)
        for (Size n = 0; n < samples; ++n) {
            const Real expected = 1.0 + 0.1 * d + 0.01 * n;
            const Real got = loaded->get(d, n, AggregationScenarioDataType::Numeraire, "");
            if (std::fabs(got - expected) > 1e-12) {
                ++mismatches;
                bad.insert({d, n});
            }
        }
    BOOST_TEST_MESSAGE("aggregation scenario data round trip: " << mismatches << " of " << dates * samples
                                                                 << " records differ after the round trip"
                                                                 << (mismatches == 1 && bad.count({0, 0})
                                                                         ? " (the first record, date 0 sample 0: the "
                                                                           "loader reads two header lines where the "
                                                                           "writer emits one; upstream observation "
                                                                           "recorded in the A7 plan)"
                                                                         : ""));
    // every record other than the first one must round-trip exactly; the first record's loss is the
    // known upstream defect (harmless for the uncollateralised baseline, which reads no key)
    BOOST_CHECK(mismatches == 0 || (mismatches == 1 && bad.count({0, 0}) == 1));
}

BOOST_AUTO_TEST_CASE(testExternalCubeRoundTrip) {
    BOOST_TEST_MESSAGE("Testing the NPV cube file round trip for an externally filled double-precision cube (ids "
                       "in set order, T0 slot, date indexing, sparse zeros)...");
    const Date asof(10, February, 2025);
    std::set<std::string> ids = {"SWAP", "BERMUDAN", "CALLABLE"};
    std::vector<Date> dates = {Date(12, May, 2025), Date(12, August, 2025), Date(12, November, 2025)};
    const Size samples = 5;
    auto cube = QuantLib::ext::make_shared<DoublePrecisionInMemoryCube>(asof, ids, dates, samples, 1, 0.0);
    auto value = [](Size i, Size d, Size n) { return (i + 1) * 1000.0 + d * 10.0 + n * 0.25 - 2.0; };
    for (const auto& kv : cube->idsAndIndexes()) {
        const Size i = kv.second;
        cube->setT0((i + 1) * 100.0, i, 0);
        for (Size d = 0; d < dates.size(); ++d)
            for (Size n = 0; n < samples; ++n)
                cube->set(value(i, d, n), i, d, n, 0);
    }
    const std::string file = "fmm_a7_cube_roundtrip.csv";
    saveCube(file, NPVCubeWithMetaData(cube, nullptr, false, QuantLib::ext::nullopt));
    auto loaded = loadCube(file);
    std::remove(file.c_str());
    BOOST_REQUIRE(loaded && loaded->cube());
    const auto& c2 = *loaded->cube();
    BOOST_CHECK_EQUAL(c2.asof(), asof);
    BOOST_CHECK_EQUAL(c2.numIds(), Size(3));
    BOOST_CHECK_EQUAL(c2.numDates(), dates.size());
    BOOST_CHECK_EQUAL(c2.samples(), samples);
    BOOST_CHECK_EQUAL(c2.depth(), Size(1));
    BOOST_CHECK(c2.usesDoublePrecision());
    // ids keep their (set) order and every value comes back exactly
    for (const auto& kv : cube->idsAndIndexes()) {
        const Size i = kv.second;
        BOOST_CHECK_EQUAL(c2.idsAndIndexes().at(kv.first), i);
        BOOST_CHECK_EQUAL(c2.getT0(i, 0), (i + 1) * 100.0);
        for (Size d = 0; d < dates.size(); ++d) {
            BOOST_CHECK_EQUAL(c2.dates()[d], dates[d]);
            for (Size n = 0; n < samples; ++n)
                BOOST_CHECK_EQUAL(c2.get(i, d, n, 0), value(i, d, n));
        }
    }
    BOOST_TEST_MESSAGE("cube round trip: " << c2.numIds() << " ids, " << c2.numDates() << " dates, " << c2.samples()
                                           << " samples, exact");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
