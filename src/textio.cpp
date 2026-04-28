#include <hwas_types.h>


bool textio::Buffer<T>::isopen() {
    return buf_size_ == 0 || !buf_ || buf_idx_ >= buf_size_;
}

int textio::Buffer<T>::append(T c) {
    if (buf_idx_ >= buf_size_ - 1)
        return -1;
    
    buf_[buf_idx_++] = c;
    return 0;
}

void textio::Buffer<T>::reset() noexcept {
    std::memset(buf_, '\0', sizeof(T) * buf_size_);
    buf_idx_ = 0;
}


//@return   0       success
//          -1      error: argument nullptr or file not open
//          -2      error: end of file not reached
//          -3      error: failed to return to file handle origin
int textio::count_lines(textio::FileIO* fio, uint64_t* count) {
    if (!fio || !textio::isopen(fio))
        return -1;

    uint64_t counter = 0;

    int prev_char = '\0';
    int c = '\0';
    FILE* fid = fio->fid_;

    // don't count empty lines
    while ((c = getchar(fid)) != EOF) {
        if (c == '\n' && prev_char != '\n') counter++;
        prev_char = c;
    }

    // EOF is serving as a end of line character
    // that is why we need to increment
    if (c == EOF && prev_char != textio::END_OF_LINE)
        counter++;

    // return file handle to the beginning of file for parsing
    if (std::fseek(fid, 0, SEEK_SET) != 0) return -3;
    if (!std::feof(EOF)) return -2;

    *count = counter;
    return 0;
}


// @return  0       success
//          1       EOF 
//          -1`     error: input argument invalid
//          -2      error: buffer overload
//          -3      error: file io error
int textio::get_line(textio::FileIO* fio, textio::Buffer<char>* cbufptr) {

    if (!fio || !textio::isopen(fio) || !cbuf || cbuf->size == 0)
        return -1;

    FILE* fid = fio->fid_;

    int c = '\0';
    while ((c = std::fgetc(fid)) != EOF) {
        if (c == textio::END_OF_LINE)
            break;

        if (cbuf->append(c) < 0)
            return -2;
    }
    
    if (cbuf->append(textio::END_OF_LINE) < 0)
        return -2;

    if (std::ferror(fid))
        return -3

    if (std::feof(fid))
        return textio::END_OF_FILE;

    return 0;
}


// @return  0           success
//          -1          error: invalid args
//          -2          error: end valid line termination
int textio::parse_line(textio::Buffer<char>* cbuf, 
        textio::Buffer<uint64_t>* ibuf,
        const char* delimiter) {

    if (!delimiter || !(cbuf->isopen() || !(ibuf->isopen()))
        return -1;

    char c = '\0';
    ibuf.reset();
    ibuf.append(0)
    for (uint64_t i = 0;
            i < cbuf->buf_size_ && (c = cbuf->buf_[i]) != '\0';
            i++) {

        if (c == delimiter) {
            cbuf->buf_[i] = '\0';
            ibuf.append(i+1);
        }

    }

    if (c != '\0' && i >= cbuf->buf_size_)
        return -2;

    return 0;
}
