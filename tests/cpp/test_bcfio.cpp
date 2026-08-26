
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <string>

#include <gtest/gtest.h>

namespace htslib {
extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
}
}

#include <bcfio.h>

const char VCF_NAME[] { "inst/exdata/geno_test_data.vcf" };
const char VCFGZ_NAME[] { "inst/exdata/geno_test_data.vcf.gz" };
const char BCF_NAME[] { "inst/exdata/geno_test_data.bcf" };

uint8_t K_FOUNDERS = 8;
uint8_t N_SAMPS = 11;

const char POS_FILE_VALID[] { "inst/exdata/pos.include" };
const char POS_FILE_DUPLICATE[] { "inst/exdata/pos_duplicates.include" };
const char POS_FILE_INVALID[] { "inst/exdata/pos_invalid.include" };


// TODO:
//      [] sample subset from file unit tests
//      [] update commented unit tests


///////////////////////////////////////////////////////////////////////////
// Test bcfio::BcfHeader
///////////////////////////////////////////////////////////////////////////

// TEST(TestBcfHeader, BcfHdrFmtGt) {
//     htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
//     bcfio::BcfHeader hdr { fid };
// 
//     EXPECT_FALSE(hdr.isnull());
//     
//     bcfio::BcfHdrAttr attr {};
// 
//     int status = hdr.get_format_attr("GT", &attr);
//     EXPECT_EQ(status, 0);
//     EXPECT_EQ(attr.number, static_cast<uint8_t>(1));
//     EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
//     EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_STR));
// 
//     if (fid) htslib::hts_close(fid);
// }
// 
// 
// TEST(TestBcfHeader, BcfHdrFmtGp) {
//     htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
//     bcfio::BcfHeader hdr { fid };
// 
//     EXPECT_FALSE(hdr.isnull());
//     
//     bcfio::BcfHdrAttr attr {};
// 
//     int status = hdr.get_format_attr("GP", &attr);
//     EXPECT_EQ(status, 0);
//     EXPECT_EQ(attr.number, static_cast<uint8_t>(3));
//     EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
//     EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));
// 
//     if (fid) htslib::hts_close(fid);
// }
// 
// TEST(TestBcfHeader, BcfHdrFmtDs) {
//     htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
//     bcfio::BcfHeader hdr { fid };
// 
//     EXPECT_FALSE(hdr.isnull());
//     
//     bcfio::BcfHdrAttr attr {};
// 
//     int status = hdr.get_format_attr("DS", &attr);
//     EXPECT_EQ(status, 0);
//     EXPECT_EQ(attr.number, static_cast<uint8_t>(1));
//     EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
//     EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));
// 
//     if (fid) htslib::hts_close(fid);
// }
// 
// 
// TEST(TestBcfHeader, BcfHdrFmtErr) {
//     htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
//     bcfio::BcfHeader hdr { fid };
// 
//     EXPECT_FALSE(hdr.isnull());
//     
//     bcfio::BcfHdrAttr attr {};
// 
//     int status = hdr.get_format_attr("DOESNOTEXIST", &attr);
//     EXPECT_NE(status, 0);
// 
//     if (fid) htslib::hts_close(fid);
// }
// 
// 
// TEST(TestBcfHeader, BcfHdrFilter) {
//     htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
//     bcfio::BcfHeader hdr { fid };
// 
//     EXPECT_FALSE(hdr.isnull());
//     
//     bcfio::BcfHdrAttr attr {};
// 
//     int status = hdr.get_filter_attr("PASS", &attr);
//     EXPECT_EQ(status, 0);
// 
//     status = hdr.get_filter_attr("PASSING", &attr);
//     EXPECT_NE(status, 0);
// 
//     if (fid) htslib::hts_close(fid);
// }
// 
// 
// TEST(TestBcfHeader, BcfHdrInfoEaf) {
//     htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
//     bcfio::BcfHeader hdr { fid };
// 
//     EXPECT_FALSE(hdr.isnull());
//     
//     bcfio::BcfHdrAttr attr {};
// 
//     int status = hdr.get_info_attr("EAF", &attr);
//     EXPECT_EQ(status, 0);
//     EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));
//     EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_VAR));
// 
//     if (fid) htslib::hts_close(fid);
// }
// 
// 
// TEST(TestBcfHeader, BcfHdrInfoErc) {
//     htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
//     bcfio::BcfHeader hdr { fid };
// 
//     EXPECT_FALSE(hdr.isnull());
//     
//     bcfio::BcfHdrAttr attr {};
// 
//     int status = hdr.get_info_attr("ERC", &attr);
//     EXPECT_EQ(status, 0);
//     EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));
//     EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_VAR));
// 
//     if (fid) htslib::hts_close(fid);
// }
// 
// 
// TEST(TestBcfHeader, BcfHdrInfoErr) {
//     htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
//     bcfio::BcfHeader hdr { fid };
// 
//     EXPECT_FALSE(hdr.isnull());
//     
//     bcfio::BcfHdrAttr attr {};
// 
//     int status = hdr.get_info_attr("NOTAINFOMEMBER", &attr);
//     EXPECT_NE(status, 0);
// 
//     if (fid) htslib::hts_close(fid);
// }
// 
// 
TEST(TestHFileRead, Initialise) {
    struct Results {
        bool is_nullptr;
        bool is_bcf;
    };

    constexpr int num_files = 4;
    bcfio::hfile_conn_t hfs[num_files] = {
        nullptr,
        bcfio::hread(VCF_NAME),
        bcfio::hread(VCFGZ_NAME),
        bcfio::hread(BCF_NAME)
    };

    Results results[num_files] = {
        { true, false },
        { false, true },
        { false, true },
        { false, true }
    };

    for (int i = 0; i < num_files; i++) {
        if (results[i].is_nullptr)
            EXPECT_TRUE(hfs[i] == nullptr);
        else
            EXPECT_FALSE(hfs[i] == nullptr);

        if (results[i].is_bcf)
            EXPECT_TRUE(hfs[i]->is_bcf());
        else if (!results[i].is_nullptr && !results[i].is_bcf)
            EXPECT_FALSE(hfs[i]->is_bcf());
    }
}

TEST(TestBcfRecord, Constructor) {
    bcfio::brec_t<float> bfloat = bcfio::BcfRecord<float>::init();
    ASSERT_TRUE(bfloat != nullptr);
    EXPECT_EQ(bfloat->data_cap, 0);
    EXPECT_EQ(bfloat->ncol, 0);
    EXPECT_EQ(bfloat->nrow, 0);
    EXPECT_TRUE(bfloat->data == nullptr);

    bcfio::brec_t<int32_t> bint32 = bcfio::BcfRecord<int32_t>::init();
    ASSERT_TRUE(bint32 != nullptr);
    EXPECT_EQ(bint32->data_cap, 0);
    EXPECT_EQ(bint32->ncol, 0);
    EXPECT_EQ(bint32->nrow, 0);
    EXPECT_TRUE(bint32->data == nullptr);

    bcfio::brec_t<int> bint = bcfio::BcfRecord<int>::init();
    ASSERT_TRUE(bint != nullptr);
    EXPECT_EQ(bint->data_cap, 0);
    EXPECT_EQ(bint->ncol, 0);
    EXPECT_EQ(bint->nrow, 0);
    EXPECT_TRUE(bint->data == nullptr);
}


TEST(TestBcfRecord, GetGenomicCoords) {
    const char VCF_NAME[] { "inst/exdata/geno_test_data.vcf" };
    const char VCFGZ_NAME[] { "inst/exdata/geno_test_data.vcf.gz" };
    const char BCF_NAME[] { "inst/exdata/geno_test_data.bcf" };

    constexpr int nfiles = 3;

    bcfio::bid_t bids[nfiles] = {
        bcfio::bread(VCF_NAME),
        bcfio::bread(VCFGZ_NAME),
        bcfio::bread(BCF_NAME)
    };

    const char true_chrom_name[] = "chr12";
    constexpr int num_pos = 8;
    int64_t true_pos[num_pos] = {
        788, 1321, 1335, 1661, 1714, 2088, 2090, 2631
    };

    bcfio::brec_t<float> brec = bcfio::BcfRecord<float>::init();
    ASSERT_TRUE(brec != nullptr);

    const char* chrom_out = nullptr;
    int64_t pos_out = -1;
    bcfio::Status status = bcfio::Status::ErrInternal;
    bcfio::Bcf* bid = nullptr;

    EXPECT_TRUE(bcfio::chrom<float>(nullptr, brec.get()) == nullptr);
    EXPECT_TRUE(bcfio::chrom<float>(nullptr, nullptr) == nullptr);

    for (int i = 0; i < nfiles; i++) {

        bid = bids[i].get();
        EXPECT_TRUE(bcfio::chrom<float>(bid, nullptr) == nullptr);

        for (int j = 0; bcfio::next_record<float>(bid, brec.get(), "DS") == bcfio::Status::Success; j++) {
            status = bcfio::pos<float>(brec.get(), &pos_out);
            EXPECT_EQ(status, bcfio::Status::Success);
            EXPECT_EQ(pos_out, true_pos[j]);

            status = bcfio::pos<float>(nullptr, &pos_out);
            EXPECT_EQ(status, bcfio::Status::ErrInvalidInput);

            status = bcfio::pos<float>(nullptr, nullptr);
            EXPECT_EQ(status, bcfio::Status::ErrInvalidInput);

            status = bcfio::pos<float>(brec.get(), nullptr);
            EXPECT_EQ(status, bcfio::Status::ErrInvalidInput);

            chrom_out = bcfio::chrom<float>(bid, brec.get());
            EXPECT_STREQ(chrom_out, true_chrom_name);
        }
    }
}

TEST(BcfRecord, ChromFailures) {
    // default constructors have nullptr vals
    constexpr int num_bids = 3;
    bcfio::bid_t bids[num_bids] = {
        nullptr,
        bcfio::bread(BCF_NAME), 
        bcfio::bread(BCF_NAME)
    };
    for (int i = 1; i < num_bids; i++)
        ASSERT_TRUE(bids[i] != nullptr);
    bids[1]->close();

    constexpr int num_brec = 3;
    bcfio::brec_t<float> brecs[num_brec] = {
        nullptr,
        bcfio::BcfRecord<float>::init(),
        bcfio::BcfRecord<float>::init()
    };
    for (int i = 1; i < num_brec; i++)
        ASSERT_TRUE(brecs[i] != nullptr);

    bcfio::Status status = bcfio::next_record<float>(bids[1].get(), 
            brecs[1].get(), 
            "DS");
    ASSERT_EQ(status, bcfio::Status::ErrBcfNotOpen);
    htslib::bcf_destroy(brecs[1]->rec);
    brecs[1]->rec = nullptr;

    status = bcfio::next_record<float>(bids[2].get(), 
            brecs[2].get(), 
            "DS");
    ASSERT_EQ(status, bcfio::Status::Success);

    const char* chrom_return_vals[num_bids * num_brec] = {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        "chr12"
    };

    bcfio::Bcf* bid = nullptr;
    bcfio::BcfRecord<float>* brec = nullptr;
    int i = 0;
    for (int ibid = 0; ibid < num_bids; ibid++) {

        if (bids[ibid] != nullptr)
            bid = bids[ibid].get();
        else
            bid = nullptr;

        for (int ibrec = 0; ibrec < num_brec; ibrec++) {
            if (brecs[ibrec] != nullptr)
                brec = brecs[ibrec].get();
            else
                brec = nullptr;

            if (chrom_return_vals[i] == nullptr)
                EXPECT_TRUE(bcfio::chrom(bid, brec) == nullptr);
            else 
                EXPECT_STREQ(bcfio::chrom(bid, brec),
                        chrom_return_vals[i]);

            i++;
        }
    }
}


////////////////////////////////////////////////////////////////////
// Test bcfio::Bcf
////////////////////////////////////////////////////////////////////


TEST(TestBcf, DefaultConstructor) {
    bcfio::bid_t bid = bcfio::bread(BCF_NAME);
    bid->close(); 
    EXPECT_FALSE(bcfio::is_open(bid.get()));
}

TEST(TestBcf, Constructor) {
    bcfio::bid_t bid = bcfio::bread(VCF_NAME);

    EXPECT_TRUE(bcfio::is_open(bid.get()));

    uint32_t n = 0;
    bcfio::Status status = bcfio::num_samples(bid.get(), &n);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(n, N_SAMPS);
}

TEST(TestBcf, OpenFailure) {
    bcfio::bid_t bid = bcfio::bread("");
    EXPECT_EQ(bid, nullptr);
}

TEST(TestBcf, OpenVCFSuccess) {
    bcfio::bid_t bid = bcfio::bread(VCF_NAME);
    ASSERT_NE(bid, nullptr);
    EXPECT_TRUE(bcfio::is_open(bid.get()));
}

TEST(TestBcf, OpenVCFGZSuccess) {
    bcfio::bid_t bid = bcfio::bread(VCFGZ_NAME);
    ASSERT_NE(bid, nullptr);
    EXPECT_TRUE(bcfio::is_open(bid.get()));
}

TEST(TestBcf, OpenBCFSuccess) {
    bcfio::bid_t bid = bcfio::bread(BCF_NAME);
    ASSERT_NE(bid, nullptr);
    EXPECT_TRUE(bcfio::is_open(bid.get()));
}

TEST(TestBcf, Closing) {
    constexpr uint32_t nfiles = 3;
    bcfio::bid_t bids[nfiles] = {
        bcfio::bread(VCF_NAME),
        bcfio::bread(VCFGZ_NAME),
        bcfio::bread(BCF_NAME)
    };

    bcfio::Bcf* bid = nullptr;
    for (int i = 0; i < nfiles; i++) {
        bid = bids[i].get();
        EXPECT_TRUE(bid->fid != nullptr);
        EXPECT_TRUE(bid->hdr != nullptr);

        bid->close();
        EXPECT_TRUE(bid->fid == nullptr);
        EXPECT_TRUE(bid->hdr == nullptr);

        // confirm that nothing happens upon repeated close
        bid->close();
        EXPECT_TRUE(bid->fid == nullptr);
        EXPECT_TRUE(bid->hdr == nullptr);
    }
}
 

TEST(TestBcf, Kfmt) {
    bcfio::bid_t bid = bcfio::bread(BCF_NAME);
    
    // DS is alt allele dosage, which is more clearly defined as 
    // the expected count of alt alleles under the trained HMM
    uint16_t k = 0;
    bcfio::Status status = bcfio::Status::ErrInternal;

    status = bcfio::k_fmt(bid.get(), "DS", &k);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(k, 1);

    status = bcfio::k_fmt(bid.get(), "HD", &k);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(k, K_FOUNDERS);

    // TODO: what happens if I submit "GT", it exists but is a string
    //      not float
    // error detection
    EXPECT_EQ(bcfio::k_fmt(bid.get(), "WRONG_ID", &k),
            bcfio::Status::ErrInvalidId);
    EXPECT_EQ(bcfio::k_fmt(bid.get(), "", &k),
            bcfio::Status::ErrInvalidId);
    EXPECT_EQ(bcfio::k_fmt(bid.get(), nullptr, &k),
            bcfio::Status::ErrInvalidInput);
    EXPECT_EQ(bcfio::k_fmt(bid.get(), "HD", nullptr),
            bcfio::Status::ErrInvalidInput);
}


TEST(TestBcf, SubsetSamplesEasy) {
    bcfio::bid_t bid = bcfio::bread(BCF_NAME);

    constexpr uint32_t nsamps = 3;
    const char* sample_subset[nsamps] = { "S01", "S07", "S08" };
    const char sample_list[] = "S01,S07,S08";

    bcfio::Status status = bcfio::Status::ErrInternal;
    status = bcfio::subset_samples(bid.get(), sample_list);
    ASSERT_EQ(status, bcfio::Status::Success);

    uint32_t n = 0;
    status = bcfio::num_samples(bid.get(), &n);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(n, nsamps);

    for (int i = 0; i < nsamps; i++)
        EXPECT_STREQ(sample_subset[i], bid->hdr->samples[i]);
}


TEST(TestBcf, SubsetSamplesWrongSampName) {
    bcfio::bid_t bid = bcfio::bread(BCF_NAME);

    constexpr uint32_t nsamps = 3;
    const char* sample_subset[nsamps] = { "S01", "S72", "S08" };

    constexpr uint32_t nsamps_correct = 2;
    const uint32_t correct_samp_idx[nsamps_correct] = { 0, 2 };
    const char sample_list[] = "S01,S08";

    bcfio::Status status = bcfio::subset_samples(bid.get(), sample_list);
    ASSERT_EQ(status, bcfio::Status::Success);

    uint32_t n = 0;
    status = bcfio::num_samples(bid.get(), &n);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(n, nsamps_correct);

    uint32_t idx;
    for (uint32_t i = 0; i < nsamps_correct; i++) {
        idx = correct_samp_idx[i];
        EXPECT_STREQ(sample_subset[idx], bid->hdr->samples[i]);
    }
}


TEST(TestBcf, SubsetSamplesNullptr) {
    bcfio::bid_t bid = bcfio::bread(BCF_NAME);

    bcfio::Status status = bcfio::subset_samples(bid.get(), nullptr);
    ASSERT_EQ(status, bcfio::Status::Success);

    uint32_t n = 0;
    status = bcfio::num_samples(bid.get(), &n);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(n, 0);
}


struct SampleStringData {
    SampleStringData(uint32_t count)
        : nsamps(count),
        samp_names(count > 0 ? new char*[count] : nullptr),
        samp_list(nullptr) {};

    ~SampleStringData() {
        if (samp_names) {
            for (uint32_t i = 0; i < nsamps; i++) {
                if (samp_names[i])
                    delete[] samp_names[i];
            }

            delete[] samp_names;
            samp_names = nullptr;
        }

        if (samp_list) {
            delete[] samp_list;
            samp_list = nullptr;
        }
    }

    uint32_t nsamps;
    char** samp_names;
    char* samp_list;
};


std::unique_ptr<SampleStringData> 
NewSampleStringData(uint32_t count, ...) {
    if (count == 0)
        return std::make_unique<SampleStringData>(0);

    std::unique_ptr<SampleStringData> data = 
        std::make_unique<SampleStringData>(count);

    va_list args;
    va_start(args, count);

    // copy words from variadic arguments to the array of C-style
    // character strings
    char* s = nullptr;
    size_t len = 0;
    size_t total_len = 0;
    uint32_t i = 0;
    for (i = 0; i < count; i++) {
        s = va_arg(args, char*);
        if (s == nullptr) break;

        // strlen does not count null character, hence + 1
        len = std::strlen(s) + 1;
        total_len += len - 1;

        data->samp_names[i] = new char[len];
        std::memset(data->samp_names[i], '\0', len);

        // strncpy(dst, src, n)
        std::strncpy(data->samp_names[i], s, len);
    }

    va_end(args);

    if (i != count)
        return nullptr;

    // make comma delimited list of sample names, simply replace
    // the null chars with comma's except for last word
    count = data->nsamps;
    data->samp_list = new char[total_len + count];
    std::memset(data->samp_list, '\0', total_len + count);

    size_t j = 0;
    for (i = 0; i < count-1; i++) {
        len = std::strlen(data->samp_names[i]) + 1;
        std::strncpy(data->samp_list + j, data->samp_names[i], len);

        data->samp_list[j + len - 1] = ',';
        j += len;
    }

    len = std::strlen(data->samp_names[i]) + 1;
    std::strncpy(data->samp_list + j, data->samp_names[i], len);
    
    return data;
}

TEST(TestBcf, SubsetSamplesSubsequentSets) {
    bcfio::bid_t bid = bcfio::bread(BCF_NAME);

    constexpr int ntests = 3;
    std::unique_ptr<SampleStringData> data[ntests] = {
        NewSampleStringData(3, "S01", "S05", "S08"),
        NewSampleStringData(2, "S01", "S08"),
        NewSampleStringData(0)
    };

    // NewSampleStringData(4, "S02", "S03", "S04", "S07")
    bcfio::Status status = bcfio::Status::ErrInternal;
    uint32_t n = 0;
    for (int i = 0; i < ntests; i++) {
        status = bcfio::subset_samples(bid.get(), 
                data[i]->samp_list);
        ASSERT_EQ(status, bcfio::Status::Success);

        status = bcfio::num_samples(bid.get(), &n);
        EXPECT_EQ(status, bcfio::Status::Success);
        EXPECT_EQ(n, data[i]->nsamps);

        for (int j = 0; j < data[i]->nsamps; j++)
            EXPECT_STREQ(data[i]->samp_names[j], 
                    bid->hdr->samples[j]);
    }
}

TEST(TestBcf, SubsetSamplesDiffSubsequentSets) {

    const char* bcf_names[3] = { VCF_NAME, VCFGZ_NAME, BCF_NAME };
    bcfio::bid_t bid;

    uint32_t n = 0;
    bcfio::Status status = bcfio::Status::ErrInternal;
    for (int bcf_i = 0; bcf_i < 3; bcf_i++) {
        bid = bcfio::bread(bcf_names[bcf_i]);
    
        constexpr int ntests = 2;
        std::unique_ptr<SampleStringData> data[ntests] = {
            NewSampleStringData(3, "S01", "S05", "S08"),
            NewSampleStringData(4, "S04", "S05", "S06", "S08")
        };
    
        status = bcfio::subset_samples(bid.get(), 
                    data[0]->samp_list);
        ASSERT_EQ(status, bcfio::Status::Success);

        status = bcfio::num_samples(bid.get(), &n);
        EXPECT_EQ(status, bcfio::Status::Success);
        EXPECT_EQ(n, data[0]->nsamps);
        
        status = bcfio::subset_samples(bid.get(), 
                    data[1]->samp_list);
        EXPECT_EQ(status, bcfio::Status::WarnSampleSetMismatch);

        bid->close();
    }
}


TEST(TestBcf, SamplesExclusion) {
    bcfio::bid_t bid = bcfio::bread(BCF_NAME);

    constexpr int ntests = 6;
    std::unique_ptr<SampleStringData> data[ntests] = {
        NewSampleStringData(5, "S01", "S02", "S03", "S05", "S08"),
        NewSampleStringData(1, "^S03"),
        NewSampleStringData(4, "S01", "S02", "S05", "S08"),
        NewSampleStringData(4, "S01", "^S02", "S05", "S08"),
        NewSampleStringData(2, "^S01", "S05"),
        NewSampleStringData(1, "S08")
    };

    bcfio::Status status = bcfio::subset_samples(bid.get(), 
                data[0]->samp_list);
    ASSERT_EQ(status, bcfio::Status::Success);

    uint32_t n = 0;
    status = bcfio::num_samples(bid.get(), &n);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(n, data[0]->nsamps);
    
    // Test simple sample exclusion
    status = bcfio::subset_samples(bid.get(), 
                data[1]->samp_list);
    ASSERT_EQ(status, bcfio::Status::Success);

    status = bcfio::num_samples(bid.get(), &n);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(n, data[2]->nsamps);

    for (int j = 0; j < data[2]->nsamps; j++)
        EXPECT_STREQ(data[2]->samp_names[j], bid->hdr->samples[j]);


    // Can't mix and match exlcusion inclusion cases
    status = bcfio::subset_samples(bid.get(), 
                data[3]->samp_list);
    ASSERT_EQ(status, bcfio::Status::WarnSampleSetMismatch);


    // Multiple sample exclusion
    status = bcfio::subset_samples(bid.get(), 
                data[4]->samp_list);
    ASSERT_EQ(status, bcfio::Status::Success);

    status = bcfio::num_samples(bid.get(), &n);
    EXPECT_EQ(status, bcfio::Status::Success);
    EXPECT_EQ(n, data[5]->nsamps);

    for (uint32_t j = 0; j < n; j++)
        EXPECT_STREQ(data[5]->samp_names[j], bid->hdr->samples[j]);
}


TEST(TestBcfInfo, GetFilename) {
    constexpr int num_files = 4;
    bcfio::bid_t bids[num_files] = {
        nullptr,
        bcfio::bread(VCF_NAME),
        bcfio::bread(VCFGZ_NAME),
        bcfio::bread(BCF_NAME)
    };
    const char* filenames[num_files] = { 
        nullptr,
        VCF_NAME, 
        VCFGZ_NAME, 
        BCF_NAME
    };

    for (int i = 0; i < num_files; i++) {
        if (bids[i] == nullptr) {
            EXPECT_TRUE(filenames[i] == nullptr);
            continue;
        }

        EXPECT_STREQ(bcfio::get_filename(bids[i].get()), filenames[i]);
    }
}


TEST(TestBcfInfo, NumPositions) {
    constexpr int num_files = 5;
    bcfio::bid_t bids[num_files] = {
        nullptr,
        bcfio::bread(VCF_NAME),
        bcfio::bread(VCFGZ_NAME),
        bcfio::bread(BCF_NAME),
        bcfio::bread(BCF_NAME)
    };
    bids[4]->close();

    int64_t npos = 8;
    int64_t ntest = 0;
    bcfio::Status statuses[num_files] = { 
        bcfio::Status::ErrBcfNotOpen, 
        bcfio::Status::Success, 
        bcfio::Status::Success, 
        bcfio::Status::Success,
        bcfio::Status::ErrBcfNotOpen
    };

    bcfio::Status status = bcfio::Status::ErrInternal;

    for (int i = 0; i < num_files; i++) {
        status = bcfio::num_pos(bids[i].get(), &ntest);

        EXPECT_EQ(status, statuses[i]);

        if (statuses[i] != bcfio::Status::Success)
            continue;

        EXPECT_EQ(ntest, npos);
     }

}

TEST(BcfUtils, IsBcf) {
    EXPECT_TRUE(bcfio::is_bcf(BCF_NAME));
    EXPECT_TRUE(bcfio::is_bcf(VCF_NAME));
}


TEST(GenomicCoord, Comparison) {
    struct GenomicCoordTestData {
        bcfio::GenomicCoord lhs;
        bcfio::GenomicCoord rhs;
        int truth_val;
    };

    // the lhs compared to the rhs encoding:
    //  -1: less than,
    //  0: equality,
    //  1: greater than,
    
    constexpr int ntests = 5;
    GenomicCoordTestData gr[ntests] = {
        { {"chr12", 2325}, {"chr12", 5422},      -1},
        { {"chr12", 2325}, {"chr12", 2325},      0 },
        { {"chr13", 2325}, {"chr12", 2325},      1 },
        { {"chr13", 2325}, {"chr12", 42325},     1 },
        { {"chr12", 2325}, {"chr13", 25},        -1}
    };

    for (int i = 0; i < ntests; i++) {
        switch (gr[i].truth_val) {
        case -1:
            EXPECT_TRUE(gr[i].lhs < gr[i].rhs);
            break;
        case 0:
            EXPECT_TRUE(gr[i].lhs == gr[i].rhs);
            break;
        case 1:
            EXPECT_TRUE(gr[i].lhs > gr[i].rhs);
            break;
        default:
            printf("Couldn't interpret the test.\n");
            ASSERT_EQ(1,2);
        }
    }
}


TEST(PositionsFile, FactoryFunction) {
    struct ValidPosTestData {
        bcfio::pos_file_t pfid;
        bool is_null;
    };

    constexpr int num_files = 3;
    ValidPosTestData pdata[num_files] = {
        { bcfio::PositionsFile::read(POS_FILE_VALID), false },
        { bcfio::PositionsFile::read(POS_FILE_DUPLICATE), false },
        { nullptr, true }
    };

    for (int i = 0; i < num_files; i++) {
        if (pdata[i].is_null) {
            EXPECT_TRUE(pdata[i].pfid == nullptr);
            EXPECT_FALSE(bcfio::is_open(pdata[i].pfid.get()));
        } else {
            EXPECT_FALSE(pdata[i].pfid == nullptr);
            EXPECT_TRUE(bcfio::is_open(pdata[i].pfid.get()));
        }
    }
}


TEST(PositionsFile, GetLineValidRecords) {
    struct GetLineTestData {
        const char* tline;
        bcfio::Status status;
    };

    bcfio::pos_file_t pfid = bcfio::PositionsFile::read(POS_FILE_VALID);
    ASSERT_TRUE(pfid != nullptr);

    constexpr int num_pos = 3;
    GetLineTestData true_vals[num_pos] {
        { "chr12:1321", bcfio::Status::Success },
        { "chr12:1714", bcfio::Status::Success },
        { "chr12:2631", bcfio::Status::Success }
    };

    bcfio::Status status = bcfio::Status::ErrInternal;
    for (int i = 0; i < num_pos; i++) {
        status = pfid->getline();
        EXPECT_EQ(status, true_vals[i].status);
        EXPECT_STREQ(pfid->buf(), true_vals[i].tline);
    }

    // repeated calls when EOF reach returns EndOfFile and empty buffer
    status = pfid->getline();
    EXPECT_EQ(status, bcfio::Status::EndOfFile);
    EXPECT_STREQ("", pfid->buf());

}

TEST(PositionsFile, NextRecordValidFile) {
    struct GetRecordTestData {
        bcfio::GenomicCoord gc;
        bcfio::Status status;
    };

    bcfio::pos_file_t pfid = bcfio::PositionsFile::read(POS_FILE_VALID);
    ASSERT_TRUE(pfid != nullptr);

    constexpr int num_pos = 3;
    GetRecordTestData true_vals[num_pos] {
        { { "chr12", 1321 }, bcfio::Status::Success },
        { { "chr12", 1714 }, bcfio::Status::Success },
        { { "chr12", 2631 }, bcfio::Status::Success }
    };

    bcfio::GenomicCoord null_gc {};
    null_gc.ctg = std::string("");
    null_gc.pos = 0;

    bcfio::GenomicCoord wrong_gc {};
    wrong_gc.ctg = std::string("wrong_ctg");
    wrong_gc.pos = 35223523;

    bcfio::GenomicCoord gc {};
    bcfio::Status status = bcfio::Status::ErrInternal;
    for (int i = 0; i < num_pos; i++) {
        status = pfid->next_record(&gc);
        EXPECT_EQ(status, true_vals[i].status);
        EXPECT_EQ(gc, true_vals[i].gc);
        EXPECT_NE(gc, wrong_gc);
        EXPECT_NE(gc, null_gc);
    }

    // repeated calls when EOF reach returns EndOfFile and empty buffer
    status = pfid->next_record(&gc);
    EXPECT_EQ(status, bcfio::Status::EndOfFile);
    EXPECT_EQ(gc, null_gc);
    EXPECT_NE(gc, wrong_gc);
}

TEST(PositionsFile, GetlineInvalidFile) {
    struct GetlineResultsData {
        const char* tline;
        bcfio::Status status;
    };

    constexpr int num_lines = 6;
    GetlineResultsData truth_vals[num_lines] = {
        { "", bcfio::Status::WarnEmptyLine },
        { "chr12:1321", bcfio::Status::Success },
        { "chr12:1714", bcfio::Status::Success },
        { "", bcfio::Status::WarnEmptyLine },
        { "", bcfio::Status::ErrParsePositionsFileCoordStrTooLong },
        { "chr12:2631", bcfio::Status::Success }
    };

    bcfio::pos_file_t pfid = bcfio::PositionsFile::read(POS_FILE_INVALID);

    bcfio::Status status = bcfio::Status::ErrInternal;
    for (int i = 0; i < num_lines; i++) {
        status = pfid->getline();
        EXPECT_EQ(status, truth_vals[i].status);
        EXPECT_STREQ(pfid->buf(), truth_vals[i].tline);
    }

    status = pfid->getline();
    EXPECT_EQ(status, bcfio::Status::EndOfFile);
    EXPECT_STREQ("", pfid->buf());
}

// TODO: test bcfio::PositionsFile::next_record for violation in 
// coordinate string model
//
// TODO: test bcfio::set_pos_from_file
//
//
//
//
// TEST(TestReadBcf, LoadRecord) {
// 
//     bcfio::ReadBcf bcf = bcfio::open(VCF_NAME, "r");
//     bcfio::BcfFloatRecord rec {};
// 
//     bcf.next_record(&rec, "HD");
// 
//     int32_t k_founders = bcf.k_fmt("HD");
//     EXPECT_FALSE(k_founders <= 0);
//     EXPECT_EQ(rec.size(), bcf.n_samples() * static_cast<uint32_t>(k_founders));
//     EXPECT_EQ(bcf.n_samples(), rec.nrows());
//     EXPECT_EQ(static_cast<uint64_t>(k_founders), rec.ncols());
//     EXPECT_TRUE(rec.is_snp());
// 
//     // EXPECT_EQ(record.chrom(), "chr12");
//     // EXPECT_EQ(record.pos(), 788);
//     // EXPECT_EQ(record.id(), ".");
//     // EXPECT_EQ(record.ref(), 'A');
//     // EXPECT_EQ(record.alt(), 'G');
//     // EXkECT_EQ(record.qual(), ".");
//     // EXPECT_EQ(record.filter(), "PASS");
//     // EXPECT_EQ(record.info(), "EAF=0.00228;INFO_SCORE=1;HWE=1;ERC=0.01949;EAC=7.94153;PAF=0.00245;REF_PANEL=0");
//     // EXPECT_EQ(record.format(), "GT:GP:DS:HD");
// }
// 
// 
// 
// struct Buff {
//     Buff(const uint32_t size_in)
//         : size(size_in), array(new char[size]) { reset(); };
//     ~Buff() { if (array) delete[] array; };
// 
//     void reset() {
//         std::memset(array, '\0', size);
//     }
// 
//     uint32_t size;
//     char* array;
// };
// 
// struct FloatArray {
//     FloatArray(const uint32_t size_in)
//         : size(size_in), array(new float[size]) { reset(); };
//     ~FloatArray() { if (array) delete[] array; };
// 
//     void reset() {
//         for (uint32_t i = 0; i < size; i++)
//             array[i] = static_cast<float>(0);
//     }
// 
//     uint32_t size;
//     float *array;
// };
// 
// 
// int get_sample_truth_vals(FILE* fid,
//         FloatArray* data, 
//         Buff* buff) {
// 
//     size_t buff_idx = 0;
//     size_t data_idx = 0;
//     int c;
//     while ((c = fgetc(fid)) != EOF) {
// 
//         if (c == ',' || c == '\n') {
//             if (buff_idx >= buff->size-1) 
//                 return -1;
//             buff->array[buff_idx] = '\0';
// 
//             if (data_idx >= data->size)
//                 return -1;
// 
//             data->array[data_idx++] = atof(buff->array);
//             buff_idx = 0;
//             buff->reset();
// 
//             if (c == '\n') break;
// 
//             continue;
//         }
// 
//         buff->array[buff_idx++] = c;
//     }
//     
//     return 0;
// }
// 
// 
// 
// TEST(TestReadBcf, HDRecordValue) {
// 
//     bcfio::ReadBcf bcf = bcfio::open(VCF_NAME, "r");
//     bcfio::BcfFloatRecord rec {};
// 
//     Buff buff_fname { 100 };
//     Buff buff_data { 1000 };
//     FloatArray data { static_cast<uint32_t>(bcf.k_fmt("HD")) };
// 
//     FILE* fid;
//     // iterate positions
//     size_t pos = 1;
//     int status;
//     while (bcf.next_record(&rec, "HD") == 0) {
//         
//         snprintf(buff_fname.array, 
//                 buff_fname.size, 
//                 "tests/hd_%02zu.csv", pos++);
// 
//         fid = fopen(buff_fname.array, "r");
// 
//         EXPECT_EQ(rec.nrows(), bcf.n_samples());
//         EXPECT_EQ(rec.ncols(), static_cast<uint64_t>(bcf.k_fmt("HD")));
// 
//         // loop over samples
//         for (uint64_t i = 0; i < rec.nrows(); i++) {
// 
//             status = get_sample_truth_vals(fid, &data, &buff_data);
//             if (status != 0)
//                 printf("\n\nERROR\n\n");
// 
//            // loop over haplotypes
//             for (uint64_t j = 0; j < rec.ncols(); j++)
//                 EXPECT_EQ(data.array[j], rec.get(i, j).value());
//         }
//         fclose(fid);
//     }
// }
