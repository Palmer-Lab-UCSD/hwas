#ifndef STATUS_H
#define STATUS_H

enum struct Status : int {
    WarnEmptyLine                           = 5,
    WarnSampleSetMismatch                   = 4,
    Success                                 = 0,
    EndOfFile                               = -1,
    ErrNotImplemented                       = -3,
    ErrHtslib                               = -4,
    ErrBcfNotOpen                           = -5,
    ErrBcfRecordInvalid                     = -6,
    ErrInternal                             = -7,
    ErrInvalidInput                         = -8,
    ErrParseBcf                             = -9,
    ErrInvalidId                            = -10,
    ErrBcfOpenFailure                       = -11,
    ErrDuplicatePositions                   = -12,
    ErrParsePositionsFileInvalidCoord       = -16,
    ErrParsePositionsFileCoordStrTooLong    = -17,
    ErrCouldNotReadFile                     = -18,
    ErrCouldNotInsertCoordInPosSet          = -19,
    ErrParseUnrecoverable                   = -20,
    ErrNotSymmetricMatrix                   = -21,
    ErrNotSquareMatrix                      = -22,
    ErrInvalidMatrix                        = -23,
    ErrHeritabilityOutOfRange               = -24,
    ErrDimensionsNotEqual                   = -25,
    ErrIndexOutofBounds                     = -26,
};

const char* status_msg(Status status);

#endif
