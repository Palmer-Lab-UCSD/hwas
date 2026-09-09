
#include <string>
#include <cstdint>
#include <cstdio>

#include <gtest/gtest.h>

#include <grm.h>


TEST(TestGrm, Inititalize) {
    grm::ugrm_t<float> float_grm_null = grm::UnnormalizedGrm<float>::init(0);
    EXPECT_TRUE(float_grm_null == nullptr);

    // verify default values
    grm::ugrm_t<float> float_grm = grm::UnnormalizedGrm<float>::init(1);
    EXPECT_EQ(float_grm->nsamples(), 1);
    EXPECT_EQ(float_grm->cap(), 1);
    EXPECT_FALSE(float_grm->is_null());

    grm::ugrm_t<int> int_grm = grm::UnnormalizedGrm<int>::init(1);
    EXPECT_EQ(int_grm->nsamples(), 1);
    EXPECT_EQ(int_grm->cap(), 1);
    EXPECT_FALSE(int_grm->is_null());
}


TEST(TestGrm, FloatConstructorValidInput) {
    constexpr uint32_t nsamps_true = 4;
    constexpr uint32_t cap_true = nsamps_true * (nsamps_true + 1) / 2;

    float float_default {};

    grm::ugrm_t<float> g = grm::UnnormalizedGrm<float>::init(nsamps_true);

    EXPECT_EQ(g->nsamples(), nsamps_true);
    EXPECT_EQ(g->cap(), cap_true);
    EXPECT_FALSE(g->is_null());

    // data should be zero-initialized
    float val = 42.0;
    grm::Status status = grm::Status::ErrInternal;

    for (uint32_t i = 0; i < g->nsamples(); i++) {
        for (uint32_t j = 0; j < g->nsamples(); j++) {
            status = g->get(i, j, &val);
            EXPECT_EQ(status, grm::Status::Success);

            EXPECT_FLOAT_EQ(val, float_default);
            val = 42.0;
        }
    }

}

TEST(TestGrm, IntConstructorValidInput) {
    constexpr uint32_t nsamps_true = 4;
    constexpr uint32_t cap_true = nsamps_true * (nsamps_true + 1) / 2;
    int int_default {};

    grm::ugrm_t<int> g = grm::UnnormalizedGrm<int>::init(nsamps_true);

    EXPECT_EQ(g->nsamples(), nsamps_true);
    EXPECT_EQ(g->cap(), cap_true);
    EXPECT_FALSE(g->is_null());

    // data should be zero-initialized
    int val = 42;
    grm::Status status = grm::Status::ErrInternal;

    for (uint32_t i = 0; i < g->nsamples(); i++) {
        for (uint32_t j = 0; j < g->nsamples(); j++) {
            status = g->get(i, j, &val);
            EXPECT_EQ(status, grm::Status::Success);

            EXPECT_EQ(val, int_default);
            val = 42;
        }
    }
}

// This test aims to also test indexing, as if indexing is wrong then
// I should not recover the correct values
TEST(TestGrm, GetSetValues) {
    // Verify the manual example from grm.h comments for n=3
    grm::ugrm_t<float> g = grm::UnnormalizedGrm<float>::init(3);

    // Upper triangle and diagonal
    EXPECT_EQ(g->set(0, 0, static_cast<float>(1)), grm::Status::Success);
    EXPECT_EQ(g->set(0, 1, static_cast<float>(2)), grm::Status::Success);
    EXPECT_EQ(g->set(0, 2, static_cast<float>(3)), grm::Status::Success);
    EXPECT_EQ(g->set(1, 1, static_cast<float>(4)), grm::Status::Success);
    EXPECT_EQ(g->set(1, 2, static_cast<float>(5)), grm::Status::Success);
    EXPECT_EQ(g->set(2, 2, static_cast<float>(6)), grm::Status::Success);

    float val = 0;
    EXPECT_EQ(g->get(0, 0, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(1));
    EXPECT_EQ(g->get(0, 1, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(2));
    EXPECT_EQ(g->get(0, 2, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(3));
    EXPECT_EQ(g->get(1, 1, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(4));
    EXPECT_EQ(g->get(1, 2, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(5));
    EXPECT_EQ(g->get(2, 2, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(6));

    // test get for lower triangle
    EXPECT_EQ(g->get(1, 0, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(2));
    EXPECT_EQ(g->get(2, 0, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(3));
    EXPECT_EQ(g->get(2, 1, &val), grm::Status::Success);
    EXPECT_FLOAT_EQ(val, static_cast<float>(5));
}


TEST(TestGrm, MidxToArrBoundsCheck) {
    grm::ugrm_t<float> g = grm::UnnormalizedGrm<float>::init(3);
    float val {};
    EXPECT_EQ(g->get(3, 0,  &val), grm::Status::ErrIndexOutofBounds);
    EXPECT_EQ(g->get(1, -1, &val), grm::Status::ErrIndexOutofBounds);
    EXPECT_EQ(g->get(3, 3,  &val), grm::Status::ErrIndexOutofBounds);
    EXPECT_EQ(g->get(1, 43, &val), grm::Status::ErrIndexOutofBounds);


    EXPECT_EQ(g->set(3, 0,  val), grm::Status::ErrIndexOutofBounds);
    EXPECT_EQ(g->set(1, -1, val), grm::Status::ErrIndexOutofBounds);
    EXPECT_EQ(g->set(3, 3,  val), grm::Status::ErrIndexOutofBounds);
    EXPECT_EQ(g->set(1, 43, val), grm::Status::ErrIndexOutofBounds);
}


TEST(TestGrm, SingleSampleMatrix) {
    grm::ugrm_t<float> g = grm::UnnormalizedGrm<float>::init(1);

    EXPECT_EQ(g->nsamples(), static_cast<uint32_t>(1));
    EXPECT_EQ(g->cap(), static_cast<uint32_t>(1));
    EXPECT_FALSE(g->is_null());

    constexpr int ntrue_vals = 4;
    float true_vals[ntrue_vals] = { 2.5, 0, 42, static_cast<float>(34354233535123) };
    float test_val {};
    grm::Status status = grm::Status::ErrInternal;

    for (int i = 0; i < ntrue_vals; i++) {
        status = g->set(0, 0, true_vals[i]);
        EXPECT_EQ(status, grm::Status::Success);

        status = g->get(0, 0, &test_val);
        EXPECT_EQ(status, grm::Status::Success);
        EXPECT_FLOAT_EQ(true_vals[i], test_val);
    }
}

