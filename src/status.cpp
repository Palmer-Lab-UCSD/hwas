#include <status.h>

const char* status_msg(Status status) {
    switch (status) {
    case Status::WarnEmptyLine:
        return "Warning: Line is empty";
    case Status::WarnSampleSetMismatch:
        return "Warning: Sample list contains names not in bcf";
    case Status::Success:
        return "Success";
    case Status::EndOfFile:
        return "Reached end of file.";
    case Status::ErrNotImplemented:
        return "Not implemented.";
    case Status::ErrHtslib:
        return "Likely a problem with htslib interface, please"
            " contact the maintainers.";
    case Status::ErrBcfNotOpen:
        return "Bcf file not open for reading";
    case Status::ErrBcfRecordInvalid:
        return "Likely invalid Bcf Record.";
    case Status::ErrInternal:
        return "Internal error, please contact maintainers.";
    case Status::ErrInvalidInput:
        return "Invalid input value";
    case Status::ErrParseBcf:
        return "Error parsing Bcf file, please check whether"
            " the file is correctly formatted.  If formatted"
            " correctly please contact maintainers.";
    case Status::ErrInvalidId:
        return "Invalid id for the bcf query";
    case Status::ErrBcfOpenFailure:
        return "Failed trying to open file, please check that"
            " the specified file is a valid vcf, vcf.gz, or bcf"
            " formatted file.";
    case Status::ErrDuplicatePositions:
        return "Positions file has dupliate positions.";
    case Status::ErrParsePositionsFileInvalidCoord:
        return "Invalid contig:pos detected in positions file.";
    case Status::ErrParsePositionsFileCoordStrTooLong:
        return "contig:pos string too long";
    case Status::ErrCouldNotReadFile:
        return "Could not open file for reading.";
    case Status::ErrCouldNotInsertCoordInPosSet:
        return "Could not add coordinate to position set, may be duplicate.";
    case Status::ErrParseUnrecoverable:
        return "File egregiously violages expected contents, exit.";
    case Status::ErrNotSymmetricMatrix:
        return "Matrix is not symmetric as required";
    case Status::ErrNotSquareMatrix:
        return "Matrix is not square as required";
    case Status::ErrInvalidMatrix:
        return "Invalid matrix";
    case Status::ErrHeritabilityOutOfRange:
        return "Heritability must be on interval (0,1).";
    case Status::ErrDimensionsNotEqual:
        return "Dimension(s) of data objects are not equal.";
    case Status::ErrIndexOutofBounds:
        return "Index out of bounds";
    default:
        break;
    }

    return "Unexpected status, please contact maintainers.";
}
