
#include <string>
#include <cstdint>
#include <cstdio>

#include <gtest/gtest.h>

#include <grm.h>


////////////////////////////////////////////////////////////////////
// COORDINATES TESTS
////////////////////////////////////////////////////////////////////

TEST(TestGrm, DefaultConstructor) {
    // verify default values
    grm::Grm<float> float_grm {};
    EXPECT_EQ(float_grm.nsamps, static_cast<uint32_t>(0));
    EXPECT_EQ(float_grm.capacity, static_cast<uint32_t>(0));
    EXPECT_TRUE(float_grm.data == nullptr);

    grm::Grm<int> int_grm {};
    EXPECT_EQ(int_grm.nsamps, static_cast<uint32_t>(0));
    EXPECT_EQ(int_grm.capacity, static_cast<uint32_t>(0));
    EXPECT_TRUE(int_grm.data == nullptr);
}


TEST(TestGrm, ConstructorValidInput) {
    constexpr uint32_t nsamps_true = 4;
    constexpr uint32_t cap_true = nsamps_true * (nsamps_true + 1) / 2;

    float float_default {};

    grm::Grm<float> g { nsamps_true };

    EXPECT_EQ(g.nsamps, nsamps_true);
    EXPECT_EQ(g.capacity, cap_true);
    EXPECT_NE(g.data, nullptr);

    // data should be zero-initialized
    for (uint64_t i = 0; i < g.capacity; i++)
        EXPECT_FLOAT_EQ(g.data[i], float_default);

    int int_default {};
    grm::Grm<int> gint { nsamps_true };

    EXPECT_EQ(gint.nsamps, nsamps_true);
    EXPECT_EQ(gint.capacity, cap_true);
    EXPECT_NE(gint.data, nullptr);

    // data should be zero-initialized
    for (uint64_t i = 0; i < gint.capacity; i++)
        EXPECT_FLOAT_EQ(gint.data[i], int_default);

}


TEST(TestGrm, ConstructorZero) {
    constexpr uint32_t nsamps_true = 0;
    constexpr uint32_t cap_true = 0;

    grm::Grm<float> gfloat { nsamps_true };

    EXPECT_EQ(gfloat.nsamps, nsamps_true);
    EXPECT_EQ(gfloat.data, nullptr);
    EXPECT_EQ(gfloat.capacity, cap_true);


    grm::Grm<int> gint { nsamps_true };

    EXPECT_EQ(gint.nsamps, nsamps_true);
    EXPECT_EQ(gint.data, nullptr);
    EXPECT_EQ(gint.capacity, cap_true);
}


TEST(TestGrm, MidxToArr) {
    // Verify the manual example from grm.h comments for n=3
    grm::Grm<float> g { 3 };
    uint32_t idx = 0;

    // Upper triangle and diagonal
    EXPECT_EQ(g.midx_to_arr(0, 0, &idx), 0);
    EXPECT_EQ(idx, static_cast<uint32_t>(0));

    EXPECT_EQ(g.midx_to_arr(0, 1, &idx), 0);
    EXPECT_EQ(idx, static_cast<uint32_t>(1));

    EXPECT_EQ(g.midx_to_arr(0, 2, &idx), 0);
    EXPECT_EQ(idx,static_cast<uint32_t>(2));

    EXPECT_EQ(g.midx_to_arr(1, 1, &idx), 0);
    EXPECT_EQ(idx,static_cast<uint32_t>(3));

    EXPECT_EQ(g.midx_to_arr(1, 2, &idx), 0);
    EXPECT_EQ(idx,static_cast<uint32_t>(4));

    EXPECT_EQ(g.midx_to_arr(2, 2, &idx), 0);
    EXPECT_EQ(idx,static_cast<uint32_t>(5));

    // Lower triangle should map to same idx by symmetry
    EXPECT_EQ(g.midx_to_arr(1, 0, &idx), 0);
    EXPECT_EQ(idx, static_cast<uint32_t>(1));

    EXPECT_EQ(g.midx_to_arr(2, 0, &idx), 0);
    EXPECT_EQ(idx, static_cast<uint32_t>(2)); 

    EXPECT_EQ(g.midx_to_arr(2, 1, &idx), 0);
    EXPECT_EQ(idx, static_cast<uint32_t>(4));
}


// TEST(TestGrm, MidxToArrBoundsCheck) {
//     grm::Grm g { 3 };
//     uint64_t idx = 0;
// 
//     EXPECT_EQ(g.midx_to_arr(3, 0, &idx), grm::ERROR_IDX_ARR_BOUNDS);
//     EXPECT_EQ(g.midx_to_arr(0, 3, &idx), grm::ERROR_IDX_ARR_BOUNDS);
//     EXPECT_EQ(g.midx_to_arr(3, 3, &idx), grm::ERROR_IDX_ARR_BOUNDS);
// }
// 
// 
// TEST(TestGrm, OperatorParensSetAndGet) {
//     grm::Grm g { 3 };
// 
//     // Set via operator()
//     g(0, 0) = 1.0f;
//     g(0, 1) = 2.0f;
//     g(0, 2) = 3.0f;
//     g(1, 1) = 4.0f;
//     g(1, 2) = 5.0f;
//     g(2, 2) = 6.0f;
// 
//     // Read back via operator() const
//     const grm::Grm& cg = g;
//     EXPECT_FLOAT_EQ(cg(0, 0), 1.0f);
//     EXPECT_FLOAT_EQ(cg(0, 1), 2.0f);
//     EXPECT_FLOAT_EQ(cg(0, 2), 3.0f);
//     EXPECT_FLOAT_EQ(cg(1, 1), 4.0f);
//     EXPECT_FLOAT_EQ(cg(1, 2), 5.0f);
//     EXPECT_FLOAT_EQ(cg(2, 2), 6.0f);
// }
// 
// 
// TEST(TestGrm, OperatorParensSymmetry) {
//     grm::Grm g { 3 };
// 
//     g(0, 1) = 7.5f;
//     g(2, 0) = 3.3f;
// 
//     const grm::Grm& cg = g;
// 
//     // (i,j) should equal (j,i) due to symmetry
//     EXPECT_FLOAT_EQ(cg(0, 1), cg(1, 0));
//     EXPECT_FLOAT_EQ(cg(0, 2), cg(2, 0));
//     EXPECT_FLOAT_EQ(cg(0, 1), 7.5f);
//     EXPECT_FLOAT_EQ(cg(0, 2), 3.3f);
// }
// 
// 
// TEST(TestGrm, SetAndGet) {
//     grm::Grm g { 3 };
// 
//     EXPECT_EQ(g.set(0, 0, 1.1f), grm::SUCCESS);
//     EXPECT_EQ(g.set(1, 2, 2.2f), grm::SUCCESS);
// 
//     float val = 0.0f;
//     EXPECT_EQ(g.get(0, 0, &val), grm::SUCCESS);
//     EXPECT_FLOAT_EQ(val, 1.1f);
// 
//     EXPECT_EQ(g.get(1, 2, &val), grm::SUCCESS);
//     EXPECT_FLOAT_EQ(val, 2.2f);
// 
//     // symmetry
//     EXPECT_EQ(g.get(2, 1, &val), grm::SUCCESS);
//     EXPECT_FLOAT_EQ(val, 2.2f);
// }
// 
// 
// TEST(TestGrm, SetGetBoundsCheck) {
//     grm::Grm g { 3 };
// 
//     EXPECT_EQ(g.set(3, 0, 1.0f), grm::ERROR_IDX_ARR_BOUNDS);
//     EXPECT_EQ(g.set(0, 3, 1.0f), grm::ERROR_IDX_ARR_BOUNDS);
// 
//     float val = 0.0f;
//     EXPECT_EQ(g.get(3, 0, &val), grm::ERROR_IDX_ARR_BOUNDS);
//     EXPECT_EQ(g.get(0, 3, &val), grm::ERROR_IDX_ARR_BOUNDS);
// }
// 
// 
// TEST(TestGrm, MoveConstructor) {
//     grm::Grm src { 3 };
//     src(0, 0) = 1.0f;
//     src(1, 2) = 5.0f;
// 
//     grm::Grm dst { std::move(src) };
// 
//     EXPECT_EQ(dst.n_samples, static_cast<uint64_t>(3));
//     EXPECT_NE(dst.data, nullptr);
//     EXPECT_FLOAT_EQ(dst(0, 0), 1.0f);
//     EXPECT_FLOAT_EQ(dst(1, 2), 5.0f);
// 
//     EXPECT_EQ(src.n_samples, static_cast<uint64_t>(0));
//     EXPECT_EQ(src.data, nullptr);
// }
// 
// 
// TEST(TestGrm, MoveAssignment) {
//     grm::Grm src { 2 };
//     src(0, 0) = 1.0f;
//     src(0, 1) = 2.0f;
//     src(1, 1) = 3.0f;
// 
//     grm::Grm dst {};
//     dst = std::move(src);
// 
//     EXPECT_EQ(dst.n_samples, static_cast<uint64_t>(2));
//     EXPECT_FLOAT_EQ(dst(0, 0), 1.0f);
//     EXPECT_FLOAT_EQ(dst(0, 1), 2.0f);
//     EXPECT_FLOAT_EQ(dst(1, 1), 3.0f);
// 
//     EXPECT_EQ(src.n_samples, static_cast<uint64_t>(0));
//     EXPECT_EQ(src.data, nullptr);
// }
// 
// 
// TEST(TestGrm, SingleSampleMatrix) {
//     grm::Grm g { 1 };
// 
//     EXPECT_EQ(g.n_samples, static_cast<uint64_t>(1));
//     EXPECT_EQ(g.size(), static_cast<uint64_t>(1));
//     EXPECT_NE(g.data, nullptr);
// 
//     g(0, 0) = 2.5f;
//     EXPECT_FLOAT_EQ(g(0, 0), 2.5f);
// 
//     float val = 0.0f;
//     EXPECT_EQ(g.get(0, 0, &val), grm::SUCCESS);
//     EXPECT_FLOAT_EQ(val, 2.5f);
// 
//     EXPECT_EQ(g.set(0, 0, 3.0f), grm::SUCCESS);
//     EXPECT_EQ(g.get(0, 0, &val), grm::SUCCESS);
//     EXPECT_FLOAT_EQ(val, 3.0f);
// }
// 
// 
// // Set values via the lower triangle (i > j), confirm they are
// // accessible via the upper triangle (i < j)
// TEST(TestGrm, SetViaLowerTriangleGetViaUpper) {
//     grm::Grm g { 4 };
// 
//     // set off-diagonal values using lower triangle indices
//     g.set(1, 0, 1.1f);
//     g.set(2, 0, 2.2f);
//     g.set(2, 1, 3.3f);
//     g.set(3, 0, 4.4f);
//     g.set(3, 1, 5.5f);
//     g.set(3, 2, 6.6f);
// 
//     // get via upper triangle
//     float val = 0.0f;
//     g.get(0, 1, &val); EXPECT_FLOAT_EQ(val, 1.1f);
//     g.get(0, 2, &val); EXPECT_FLOAT_EQ(val, 2.2f);
//     g.get(1, 2, &val); EXPECT_FLOAT_EQ(val, 3.3f);
//     g.get(0, 3, &val); EXPECT_FLOAT_EQ(val, 4.4f);
//     g.get(1, 3, &val); EXPECT_FLOAT_EQ(val, 5.5f);
//     g.get(2, 3, &val); EXPECT_FLOAT_EQ(val, 6.6f);
// }
// 
// 
// // Validate indexing for a 5x5 matrix
// TEST(TestGrm, LargerMatrixIndexing) {
//     uint64_t n = 5;
//     grm::Grm g { n };
// 
//     // fill every element with a unique value
//     float counter = 1.0f;
//     for (uint64_t i = 0; i < n; i++)
//         for (uint64_t j = i; j < n; j++)
//             g.set(i, j, counter++);
// 
//     // verify all values
//     counter = 1.0f;
//     float val = 0.0f;
//     for (uint64_t i = 0; i < n; i++) {
//         for (uint64_t j = i; j < n; j++) {
//             g.get(i, j, &val);
//             EXPECT_FLOAT_EQ(val, counter++);
//         }
//     }
// }
// 
// 
// // Verify that move assignment replaces a non-empty Grm
// TEST(TestGrm, MoveAssignmentReplacesExisting) {
//     grm::Grm dst { 2 };
//     dst(0, 0) = 99.0f;
// 
//     grm::Grm src { 4 };
//     src(0, 0) = 1.0f;
//     src(3, 3) = 42.0f;
// 
//     dst = std::move(src);
// 
//     EXPECT_EQ(dst.n_samples, static_cast<uint64_t>(4));
//     EXPECT_FLOAT_EQ(dst(0, 0), 1.0f);
//     EXPECT_FLOAT_EQ(dst(3, 3), 42.0f);
// }
// 
// 
// ////////////////////////////////////////////////////////////////////
// // MACRO TEST
// ////////////////////////////////////////////////////////////////////
// 
// TEST(TestMatrixMacro, ManualValidation) {
//     // Validate MATRIX_IDX_TO_ARRAY against the worked example
//     // in grm.h for a 3x3 matrix
//     uint64_t n = 3;
//     EXPECT_EQ(MATRIX_IDX_TO_ARRAY(0, 0, n), static_cast<uint64_t>(0));
//     EXPECT_EQ(MATRIX_IDX_TO_ARRAY(0, 1, n), static_cast<uint64_t>(1));
//     EXPECT_EQ(MATRIX_IDX_TO_ARRAY(0, 2, n), static_cast<uint64_t>(2));
//     EXPECT_EQ(MATRIX_IDX_TO_ARRAY(1, 1, n), static_cast<uint64_t>(3));
//     EXPECT_EQ(MATRIX_IDX_TO_ARRAY(1, 2, n), static_cast<uint64_t>(4));
//     EXPECT_EQ(MATRIX_IDX_TO_ARRAY(2, 2, n), static_cast<uint64_t>(5));
// }
// 
// 
// TEST(TestMatrixMacro, LargerMatrix) {
//     // 4x4 matrix: upper triangle has 10 elements (0..9)
//     uint64_t n = 4;
//     uint64_t expected = 0;
//     for (uint64_t i = 0; i < n; i++)
//         for (uint64_t j = i; j < n; j++)
//             EXPECT_EQ(MATRIX_IDX_TO_ARRAY(i, j, n), expected++);
// }
// 
// 
// ////////////////////////////////////////////////////////////////////
// // FULL GRM FILE WRITE/READ TESTS
// ////////////////////////////////////////////////////////////////////
// 
// // Helper to build a complete Hdr + Grm for file I/O tests
// static void build_test_data(grm::Hdr* hdr, grm::Grm* g) {
//     hdr->grm_type = grm::EHC;
// 
//     char contig[] = "chr1";
//     *hdr->coords = grm::Coordinates { contig, 2 };
//     hdr->coords->pos[0] = static_cast<uint64_t>(50);
//     hdr->coords->pos[1] = static_cast<uint64_t>(150);
// 
//     *hdr->samples = grm::Samples { 3 };
//     hdr->samples->names[0] = "s1";
//     hdr->samples->names[1] = "s2";
//     hdr->samples->names[2] = "s3";
// 
//     *g = grm::Grm { 3 };
//     g->set(0, 0, 1.0f);
//     g->set(0, 1, 0.5f);
//     g->set(0, 2, 0.2f);
//     g->set(1, 1, 1.0f);
//     g->set(1, 2, 0.3f);
//     g->set(2, 2, 1.0f);
// }
// 
// 
// TEST(TestGrmFile, WriteReadRoundTrip) {
//     grm::Hdr hdr_w {};
//     grm::Grm grm_w {};
//     build_test_data(&hdr_w, &grm_w);
// 
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
// 
//     ASSERT_EQ(grm::write(&fio, &hdr_w, &grm_w), grm::SUCCESS);
// 
//     rewind(fio.fid);
// 
//     grm::Hdr hdr_r {};
//     grm::Grm grm_r {};
//     ASSERT_EQ(grm::read(&fio, &hdr_r, &grm_r), grm::SUCCESS);
// 
//     // verify header
//     EXPECT_EQ(hdr_r.grm_type, grm::EHC);
//     EXPECT_EQ(hdr_r.coords->contig, "chr1");
//     EXPECT_EQ(hdr_r.coords->len, static_cast<uint64_t>(2));
//     EXPECT_EQ(hdr_r.samples->len, static_cast<uint64_t>(3));
//     EXPECT_EQ(hdr_r.samples->names[0], "s1");
//     EXPECT_EQ(hdr_r.samples->names[1], "s2");
//     EXPECT_EQ(hdr_r.samples->names[2], "s3");
// 
//     // verify grm data
//     EXPECT_EQ(grm_r.n_samples, static_cast<uint64_t>(3));
//     float val = 0.0f;
//     grm_r.get(0, 0, &val); EXPECT_FLOAT_EQ(val, 1.0f);
//     grm_r.get(0, 1, &val); EXPECT_FLOAT_EQ(val, 0.5f);
//     grm_r.get(0, 2, &val); EXPECT_FLOAT_EQ(val, 0.2f);
//     grm_r.get(1, 1, &val); EXPECT_FLOAT_EQ(val, 1.0f);
//     grm_r.get(1, 2, &val); EXPECT_FLOAT_EQ(val, 0.3f);
//     grm_r.get(2, 2, &val); EXPECT_FLOAT_EQ(val, 1.0f);
// }
// 
// 
// TEST(TestGrmFile, ReadBadMagicNumber) {
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
// 
//     // write a wrong magic number
//     uint32_t bad_magic = 0x12345678;
//     fwrite(&bad_magic, sizeof(bad_magic), 1, fio.fid);
//     rewind(fio.fid);
// 
//     grm::Hdr hdr {};
//     grm::Grm g {};
//     EXPECT_EQ(grm::read(&fio, &hdr, &g), grm::ERROR_NOT_A_GRM_FILE);
// }
// 
// 
// TEST(TestGrmFile, WriteNullArgs) {
//     grm::Hdr hdr {};
//     grm::Grm g { 2 };
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
// 
//     EXPECT_EQ(grm::write(nullptr, &hdr, &g), grm::ERROR_NULLPTR_ARG);
//     EXPECT_EQ(grm::write(&fio, static_cast<const grm::Hdr*>(nullptr), &g),
//               grm::ERROR_NULLPTR_ARG);
//     EXPECT_EQ(grm::write(&fio, &hdr, static_cast<const grm::Grm*>(nullptr)),
//               grm::ERROR_NULLPTR_ARG);
// }
// 
// 
// // The read function for grm::read(fio, hdr, grm) should return
// // ERROR_NULLPTR_ARG for null pointer arguments (consistent with write).
// TEST(TestGrmFile, ReadNullArgs) {
//     grm::Hdr hdr {};
//     grm::Grm g {};
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
// 
//     EXPECT_EQ(grm::read(nullptr, &hdr, &g), grm::ERROR_NULLPTR_ARG);
//     EXPECT_EQ(grm::read(&fio, static_cast<grm::Hdr*>(nullptr), &g),
//               grm::ERROR_NULLPTR_ARG);
//     EXPECT_EQ(grm::read(&fio, &hdr, static_cast<grm::Grm*>(nullptr)),
//               grm::ERROR_NULLPTR_ARG);
// }
// 
// 
// // Reading from an empty file should fail
// TEST(TestGrmFile, ReadEmptyFile) {
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
// 
//     grm::Hdr hdr {};
//     grm::Grm g {};
//     EXPECT_EQ(grm::read(&fio, &hdr, &g), grm::ERROR_ON_READ);
// }
// 
// 
// TEST(TestGrmFile, WriteNullFid) {
//     grm::Hdr hdr {};
//     grm::Grm g { 1 };
//     io::FileIO fio { nullptr };
// 
//     EXPECT_EQ(grm::write(&fio, &hdr, &g), grm::ERROR_NULLPTR_ARG);
// }
// 
// 
// TEST(TestGrmFile, ReadNullFid) {
//     grm::Hdr hdr {};
//     grm::Grm g {};
//     io::FileIO fio { nullptr };
// 
//     EXPECT_EQ(grm::read(&fio, &hdr, &g), grm::ERROR_NULLPTR_ARG);
// }
// 
// 
// // Verify the magic number (FILE_TYPE_SPEC) is the first 4 bytes written
// TEST(TestGrmFile, MagicNumberWrittenFirst) {
//     grm::Hdr hdr {};
//     grm::Grm g {};
//     build_test_data(&hdr, &g);
// 
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
//     ASSERT_EQ(grm::write(&fio, &hdr, &g), grm::SUCCESS);
// 
//     rewind(fio.fid);
// 
//     uint32_t magic = 0;
//     size_t nread = fread(&magic, sizeof(magic), 1, fio.fid);
//     ASSERT_EQ(nread, static_cast<size_t>(1));
//     EXPECT_EQ(magic, grm::FILE_TYPE_SPEC);
// }
// 
// 
// // Single sample GRM file round-trip
// TEST(TestGrmFile, WriteReadSingleSample) {
//     grm::Hdr hdr {};
//     hdr.grm_type = grm::DS;
// 
//     char contig[] = "chr1";
//     *hdr.coords = grm::Coordinates { contig, 1 };
//     hdr.coords->pos[0] = 500;
// 
//     *hdr.samples = grm::Samples { 1 };
//     hdr.samples->names[0] = "lone_sample";
// 
//     grm::Grm g { 1 };
//     g.set(0, 0, 0.75f);
// 
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
//     ASSERT_EQ(grm::write(&fio, &hdr, &g), grm::SUCCESS);
// 
//     rewind(fio.fid);
// 
//     grm::Hdr hdr_r {};
//     grm::Grm g_r {};
//     ASSERT_EQ(grm::read(&fio, &hdr_r, &g_r), grm::SUCCESS);
// 
//     EXPECT_EQ(hdr_r.grm_type, grm::DS);
//     EXPECT_EQ(hdr_r.samples->len, static_cast<uint64_t>(1));
//     EXPECT_EQ(hdr_r.samples->names[0], "lone_sample");
// 
//     EXPECT_EQ(g_r.n_samples, static_cast<uint64_t>(1));
//     float val = 0.0f;
//     g_r.get(0, 0, &val);
//     EXPECT_FLOAT_EQ(val, 0.75f);
// }
// 
// 
// // Read from a file that has the correct magic number but is
// // truncated after it (no header data follows)
// TEST(TestGrmFile, ReadTruncatedAfterMagic) {
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
// 
//     fwrite(&grm::FILE_TYPE_SPEC, sizeof(grm::FILE_TYPE_SPEC), 1, fio.fid);
//     rewind(fio.fid);
// 
//     grm::Hdr hdr {};
//     grm::Grm g {};
//     EXPECT_EQ(grm::read(&fio, &hdr, &g), grm::ERROR_ON_READ);
// }
// 
// 
// // Write then read two independent GRM files sequentially to the same
// // file, verifying both are recovered correctly
// TEST(TestGrmFile, WriteReadTwoSequential) {
//     grm::Hdr hdr1 {};
//     hdr1.grm_type = grm::EHC;
//     char c1[] = "chr1";
//     *hdr1.coords = grm::Coordinates { c1, 1 };
//     hdr1.coords->pos[0] = 10;
//     *hdr1.samples = grm::Samples { 2 };
//     hdr1.samples->names[0] = "a";
//     hdr1.samples->names[1] = "b";
// 
//     grm::Grm g1 { 2 };
//     g1.set(0, 0, 1.0f);
//     g1.set(0, 1, 0.5f);
//     g1.set(1, 1, 1.0f);
// 
//     grm::Hdr hdr2 {};
//     hdr2.grm_type = grm::EAC;
//     char c2[] = "chr2";
//     *hdr2.coords = grm::Coordinates { c2, 1 };
//     hdr2.coords->pos[0] = 20;
//     *hdr2.samples = grm::Samples { 2 };
//     hdr2.samples->names[0] = "x";
//     hdr2.samples->names[1] = "y";
// 
//     grm::Grm g2 { 2 };
//     g2.set(0, 0, 2.0f);
//     g2.set(0, 1, 0.8f);
//     g2.set(1, 1, 2.0f);
// 
//     io::FileIO fio { tmpfile() };
//     ASSERT_NE(fio.fid, nullptr);
// 
//     ASSERT_EQ(grm::write(&fio, &hdr1, &g1), grm::SUCCESS);
//     ASSERT_EQ(grm::write(&fio, &hdr2, &g2), grm::SUCCESS);
// 
//     rewind(fio.fid);
// 
//     // read first
//     grm::Hdr r1 {};
//     grm::Grm rg1 {};
//     ASSERT_EQ(grm::read(&fio, &r1, &rg1), grm::SUCCESS);
//     EXPECT_EQ(r1.grm_type, grm::EHC);
//     EXPECT_EQ(r1.samples->names[0], "a");
//     float val = 0.0f;
//     rg1.get(0, 1, &val);
//     EXPECT_FLOAT_EQ(val, 0.5f);
// 
//     // read second
//     grm::Hdr r2 {};
//     grm::Grm rg2 {};
//     ASSERT_EQ(grm::read(&fio, &r2, &rg2), grm::SUCCESS);
//     EXPECT_EQ(r2.grm_type, grm::EAC);
//     EXPECT_EQ(r2.samples->names[0], "x");
//     rg2.get(0, 1, &val);
//     EXPECT_FLOAT_EQ(val, 0.8f);
// }
