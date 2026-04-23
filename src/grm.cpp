// Palmer Lab at UCSD
//
//
// ACKNOWLEDGMENT
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Opus 4.6, the AI assistant from Anthropic
// with minor recommendations incorporated.
//
//

#include <grm.hpp>


////////////////////////////////////////////////////////////////////
// COORDINATES CLASS
////////////////////////////////////////////////////////////////////

grm::Coordinates::Coordinates(Coordinates&& other)
    : contig(std::move(other.contig)),
    len(other.len), 
    pos(std::move(other.pos)) {

    other.len = 0;
    other.contig = "";
}

grm::Coordinates& grm::Coordinates::operator=(Coordinates&& other) {
    if (this == &other)
       return *this; 

    len = other.len;
    contig = std::move(other.contig);
    pos = std::move(other.pos);

    other.len = 0;
    other.contig = "";

    return *this;
}

// remember that Coordinates* should be uninstantiated
grm::STATUS grm::write(io::FileIO* fio, const Coordinates* coords) {

    if (!fio || !fio->fid || !coords)
        return grm::ERROR_NULLPTR_ARG;

    size_t nwritten = 0;

    // write contig name to file
    uint64_t nchar = coords->contig.size();
    nwritten = fwrite(&nchar, sizeof(nchar), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    nwritten = fwrite(coords->contig.c_str(), 
            sizeof(char), 
            nchar,
            fio->fid);
    if (static_cast<uint64_t>(nwritten) != nchar)
        return grm::ERROR_ON_WRITE;

    // write positions
    uint64_t npos = coords->len;
    nwritten = fwrite(&npos, sizeof(npos), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    nwritten = fwrite(coords->pos.get(), 
            sizeof(uint64_t),
            npos,
            fio->fid);
    if (static_cast<uint64_t>(nwritten) != npos)
        return grm::ERROR_ON_WRITE;

    return grm::SUCCESS;
}


grm::STATUS grm::read(io::FileIO* fio, Coordinates* coords) {

    if (!fio || !fio->fid || !coords)
        return grm::ERROR_NULLPTR_ARG;

    // I create a temporary Coordinates class, because I don't want
    // the input coords instance to partial update upon an error
    grm::Coordinates tmpc {};

    // will store the number of bytes read at each step
    size_t nread = 0;

    // read contig name
    uint64_t size_contig_name = 0;
    nread = fread(&size_contig_name, sizeof(size_contig_name), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;

    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(size_contig_name + 1);
    std::memset(buffer.get(), '\0', size_contig_name + 1);

    nread = fread(buffer.get(), sizeof(char), size_contig_name, fio->fid);
    if (static_cast<uint64_t>(nread) != size_contig_name)
        return grm::ERROR_ON_READ;

    tmpc.contig = std::string(buffer.get(), size_contig_name);

    // read in positions
    uint64_t npos = 0;
    nread = fread(&npos, sizeof(npos), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;

    tmpc.len = npos;

    tmpc.pos = std::make_unique<uint64_t[]>(npos);
    nread = fread(tmpc.pos.get(), sizeof(uint64_t), npos, fio->fid);    
    if (static_cast<uint64_t>(nread) != npos)
        return grm::ERROR_ON_READ;

    *coords = std::move(tmpc);

    return grm::SUCCESS;
}

////////////////////////////////////////////////////////////////////
// SAMPLES CLASS
////////////////////////////////////////////////////////////////////

grm::Samples::Samples(grm::Samples&& other) 
    : len(other.len), names(std::move(other.names)) {
        other.len = 0;
}


grm::Samples& grm::Samples::operator=(grm::Samples&& other) {
    if (this == &other)
        return *this;
    
    len = other.len;
    names = std::move(other.names);
    
    other.len = 0;
    
    return *this;
}


grm::STATUS grm::write(io::FileIO* fio, const grm::Samples* samples) {

    if (!fio || !fio->fid || !samples)
        return grm::ERROR_NULLPTR_ARG;

    size_t nwritten = 0;
    uint64_t nsamps = samples->len;
    nwritten = fwrite(&nsamps, sizeof(nsamps), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    // When it comes time to read the data, I need to make a character
    // buffer to temporarily place the read string.  To make this
    // buffer, I need to know the length of string with the greatest
    // number of characters.  Here I find that number and store in
    // the binary file.
    
    uint64_t nchar_max = 0;
    uint64_t tmp = 0;
    for (uint64_t n = 0; n < nsamps; n++)
        if ((tmp = samples->names[n].size()) > nchar_max) nchar_max = tmp;

    if (nchar_max == 0)
        return grm::ERROR_INVALID_ARG;

    nwritten = fwrite(&nchar_max, sizeof(nchar_max), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    // Write each string to file;
    uint64_t nchar = 0;
    for (uint64_t n = 0; n < nsamps; n++) {
        nchar = samples->names[n].size(); 

        nwritten = fwrite(&nchar, sizeof(nchar), 1, fio->fid);
        if (nwritten != 1)
            return grm::ERROR_ON_WRITE;

        nwritten = fwrite(samples->names[n].c_str(),
                sizeof(char), 
                nchar,
                fio->fid);
        if (nwritten != nchar)
            return grm::ERROR_ON_WRITE;
    }

    return grm::SUCCESS;
}

grm::STATUS grm::read(io::FileIO* fio, grm::Samples* samples) {

    if (!fio || !fio->fid || !samples)
        return grm::ERROR_NULLPTR_ARG;

    // I create a temporary Sample class, because I don't want the 
    // input samples instance to partial update upon an error
    grm::Samples tmp_samps {};

    size_t nread;
    uint64_t n_samples = 0;

    nread = fread(&n_samples, sizeof(n_samples), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;
    
    tmp_samps.len = n_samples;
    tmp_samps.names = std::make_unique<std::string[]>(n_samples); 

    // Get the number of characters of the longest string
    size_t nchar_max = 0;
    nread = fread(&nchar_max, sizeof(size_t), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;


    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(nchar_max+1);
    std::memset(buffer.get(), '\0', nchar_max + 1);

    // read sample names
    uint64_t nchar = 0;
    for (uint64_t n = 0; n < n_samples; n++) {
        nread = fread(&nchar, sizeof(nchar), 1, fio->fid);
        if (nread != 1)
            return grm::ERROR_ON_READ;
        
        nread = fread(buffer.get(), sizeof(char), nchar, fio->fid);
        if (static_cast<uint64_t>(nread) != nchar)
            return grm::ERROR_ON_READ;

        tmp_samps.names[n] = std::string(buffer.get(), nchar);

        std::memset(buffer.get(), '\0', nchar);
        nchar = 0;
    }

    *samples = std::move(tmp_samps);

    return grm::SUCCESS;
}

////////////////////////////////////////////////////////////////////
// HDR CLASS
////////////////////////////////////////////////////////////////////

grm::Hdr::Hdr() 
    : prog_version(constants::PROG_VERSION),
    file_version(grm::FILE_VERSION),
    grm_type(UNSPECIFIED), 
    coords(std::make_unique<Coordinates>()),
    samples(std::make_unique<Samples>()) {};


grm::Hdr::Hdr(Hdr&& other)
    : prog_version(other.prog_version),
    file_version(other.file_version),
    grm_type(other.grm_type),
    coords(nullptr), samples(nullptr) {

    other.prog_version = constants::PROG_VERSION;
    other.file_version = grm::FILE_VERSION;
    other.grm_type = grm::UNSPECIFIED;

    coords = std::move(other.coords);
    samples = std::move(other.samples);
}

grm::Hdr& grm::Hdr::operator=(Hdr&& other) {
    if (this == &other)
        return *this;

    prog_version = other.prog_version;
    other.prog_version = constants::PROG_VERSION;

    file_version = other.file_version;
    other.file_version = FILE_VERSION;

    grm_type = other.grm_type;
    other.grm_type = grm::UNSPECIFIED;

    coords = std::move(other.coords);
    samples = std::move(other.samples);
    return *this;
}


grm::STATUS grm::write(io::FileIO* fio, const Hdr* hdr) {

    if (!fio || !fio->fid || !hdr)
        return grm::ERROR_NULLPTR_ARG;

    size_t nwritten = 0;
    uint32_t tmp_version = hdr->prog_version.pack();
    nwritten = fwrite(&tmp_version, sizeof(tmp_version), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    tmp_version = hdr->file_version.pack();
    nwritten = fwrite(&tmp_version, sizeof(tmp_version), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    nwritten = fwrite(&hdr->grm_type, sizeof(GrmType), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    grm::STATUS status = grm::FAILED;
    if ((status = write(fio, hdr->coords.get())) != grm::SUCCESS)
        return status;

    if ((status = write(fio, hdr->samples.get())) != grm::SUCCESS)
        return status;

    return grm::SUCCESS;
}

grm::STATUS grm::read(io::FileIO* fio, Hdr* hdr) {

    if (!fio || !fio->fid || !hdr)
        return grm::ERROR_NULLPTR_ARG;

    Hdr tmp_hdr {};

    size_t nread = 0;
    uint32_t version = 0;
    nread = fread(&version, sizeof(version), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;
    tmp_hdr.prog_version = utils::Version::unpack(version);

    nread = fread(&version, sizeof(version), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;
    tmp_hdr.file_version = utils::Version::unpack(version);

    nread = fread(&tmp_hdr.grm_type, sizeof(grm::GrmType), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;

    grm::STATUS status = grm::FAILED;
    status = read(fio, tmp_hdr.coords.get());
    if (status != grm::SUCCESS)
        return status;

    status = read(fio, tmp_hdr.samples.get());
    if (status != grm::SUCCESS)
        return status;

    *hdr = std::move(tmp_hdr);
    return grm::SUCCESS;
}


////////////////////////////////////////////////////////////////////
// GRM CLASS
////////////////////////////////////////////////////////////////////
//
// Recall that the GRM is a symmetric matrix, therefore we only need
// to store the upper triagonal and diagonal element values.  
// Consequently, the size of the array storing the data is n*(n +1)/2.
//

grm::Grm::Grm(): n_samples(0), data(nullptr) {};

grm::Grm::Grm(uint64_t n_samps)
    : n_samples(n_samps),
        data(size() != 0 ? std::make_unique<float[]>(size()) : nullptr) {

    if (data)
       std::memset(data.get(), 0, size() * sizeof(float));
}

grm::Grm::Grm(grm::Grm&& other)
    : n_samples(other.n_samples), data(std::move(other.data)) {
    other.n_samples = 0;    
};

grm::Grm& grm::Grm::operator=(grm::Grm&& other) {
    if (this == &other)
        return *this;

    n_samples = other.n_samples;
    other.n_samples = 0;

    data = std::move(other.data);

    return *this;
}

// The number of upper diagonal + diagonal elements of the GRM
uint64_t grm::Grm::size() const { 
    return n_samples * (n_samples + 1) / 2; 
};


grm::STATUS grm::Grm::midx_to_arr(const uint64_t i, const uint64_t j,
        uint64_t* idx) const {

    if (i >= n_samples || j >= n_samples)
        return grm::ERROR_IDX_ARR_BOUNDS;

    // remember that by symmetry, the matrix is equal to its transpose
    if (i <= j)
        *idx = MATRIX_IDX_TO_ARRAY(i, j, n_samples);
    else 
        *idx = MATRIX_IDX_TO_ARRAY(j, i, n_samples);

    return grm::SUCCESS;
}


float grm::Grm::operator()(const uint64_t i, const uint64_t j) const {
    if (i > j)
        return data[MATRIX_IDX_TO_ARRAY(j, i, n_samples)];
    return data[MATRIX_IDX_TO_ARRAY(i, j, n_samples)];
}


float& grm::Grm::operator()(const uint64_t i, const uint64_t j) {
    if (i > j)
        return data[MATRIX_IDX_TO_ARRAY(j, i, n_samples)];
    return data[MATRIX_IDX_TO_ARRAY(i, j, n_samples)];
}


grm::STATUS grm::Grm::get(const uint64_t i, const uint64_t j, float *val) const {
    uint64_t idx = 0;
    grm::STATUS status = grm::FAILED;
    if ((status = midx_to_arr(i, j, &idx)) != grm::SUCCESS) 
        return status;

    *val = data[idx];

    return status;
}


grm::STATUS grm::Grm::set(const uint64_t i, const uint64_t j, const float val) {
    uint64_t idx = 0;
    grm::STATUS status = grm::FAILED;
    if ((status = midx_to_arr(i, j, &idx)) != grm::SUCCESS) 
        return status;

    data[idx] = val;

    return status;
}


grm::STATUS grm::write(io::FileIO *fio, 
        const grm::Hdr* hdr, const grm::Grm* grmatrix) {
    
    if (!fio || !fio->fid || !hdr || !grmatrix)
        return grm::ERROR_NULLPTR_ARG;

    size_t nwritten = 0;

    // write file type specification
    nwritten = fwrite(&FILE_TYPE_SPEC, 
            sizeof(FILE_TYPE_SPEC), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    grm::STATUS status = grm::FAILED;

    // srite meta data stored in header;
    if ((status = grm::write(fio, hdr)) != grm::SUCCESS)
        return status;

    uint64_t ndata = static_cast<uint64_t>(grmatrix->size());
    nwritten = fwrite(&ndata, sizeof(ndata), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    nwritten = fwrite(grmatrix->data.get(), 
            sizeof(grmatrix->data[0]),
            ndata,
            fio->fid);
    if (static_cast<uint64_t>(nwritten) != ndata)
        return grm::ERROR_ON_WRITE;

    return grm::SUCCESS;
}


grm::STATUS grm::read(io::FileIO* fio,
        grm::Hdr* hdr, grm::Grm* grmatrix) {
    
    if (!fio || !fio->fid || !hdr || !grmatrix)
        return grm::ERROR_NULLPTR_ARG;

    size_t nread = 0;
    uint32_t ftype = 0;

    nread = fread(&ftype, sizeof(ftype), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;

    if (ftype != FILE_TYPE_SPEC)
        return grm::ERROR_NOT_A_GRM_FILE;

    grm::STATUS status = grm::read(fio, hdr);
    if (status != grm::SUCCESS)
        return status;

    uint64_t n_samples = hdr->samples->len;
    grm::Grm tmp_grm { n_samples };
    uint64_t ndata = 0;
    nread = fread(&ndata, sizeof(ndata), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;

    if (ndata != n_samples*(n_samples + 1) / 2)
        return grm::ERROR_ON_READ;

    nread = fread(tmp_grm.data.get(), 
            sizeof(tmp_grm.data[0]), ndata, fio->fid);

    if (nread != ndata)
        return grm::ERROR_ON_READ;

    *grmatrix = std::move(tmp_grm);

    return grm::SUCCESS;
}
