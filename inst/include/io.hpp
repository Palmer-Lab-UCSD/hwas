#ifndef HEADER_IO_H
#define HEADER_IO_H

#include <cstdio>
#include <string>
#include <memory>

namespace io {

enum STATUS { 
    SUCCESS,
    FERROR, 
    FEOF, 
    INVALID_ARG_ERROR, 
    FSEEK_ERROR, 
    FEOF_ERROR,
    END_OF_BUF_ERROR
};


struct FileIO {
    FileIO(FILE* fid): fid(fid) {}; 
    ~FileIO() { if (fid) { fclose(fid); fid = nullptr; } };

    FileIO(const FileIO&) = delete;
    FileIO& operator=(const FileIO&) = delete;
    
    FileIO(FileIO&& other) noexcept : fid(other.fid) {
        other.fid = nullptr;
    }

    FileIO& operator=(FileIO&& other) {
        // protect against self assignment
        if (this == &other) return *this;

        if (fid) fclose(fid);

        fid = other.fid;
        other.fid = nullptr;

        return *this;
    }

    FILE *fid;
};


// @title: open a file 
// @description:
// @param filename: name and path of file to open
// @param mode: a mode in the set of those in the C library function fopen
// @return a pointer to opened file
inline FileIO open(const char *filename, const char *mode) {
    if (!mode  || !filename)
        return nullptr;

    FILE *fid = fopen(filename, mode);
    if (!fid)
        return nullptr;

    // compiler implements elision
    return FileIO(fid);
}


// @title: file object
// @description: This class manages the lifetime of a C-style file stream
//      by RAII.  To contruct an instance of the class use the "open" function
//      below.
// @param fid: an opened C-style file stream
// int bseek(FileIO *fio);



// @title: File statistics
// @description: This object is returned by any function meant to calculate
//      file character statistics.
// struct FileStats {
//     size_t nchar = 0;
//     size_t nwords = 0;
//     size_t nlines = 0;
//     size_t nblanklines = 0;
// }


// @title: word count
// @description: Similar to the UNIX/Linux wc command line program, wc
//      calculates the number of characters, words, lines, etc. that
//      the specified file contains.
// @param tio: an instance of TextIO
// @param fs: the structure that the file statistics will be stored
// @return a STATUS code that specifies whether the function was successful
//      or failed.
// STATUS wc(TextIO *tio, FileStats *fs);




// STATUS getline(TextIO *tio, Array<char> linebuf);
}

#endif
