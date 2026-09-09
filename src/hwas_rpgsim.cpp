#include <Rpgsim.h>


const char* rpgsim::status_msg(rpgsim::Status status) {
    switch (status) {
    case rpgsim::Status::Success:
        return "Success";
    case rpgsim::Status::ErrNotSymmetricMatrix:
        return "Not a symmetric matrix";
    case rpgsim::Status::ErrHeritabilityOutOfRange:
        return "Heritability must be on interval (0,1).";
    default:
        break;
    }

    return "Unexpected status, please contact maintainers.";
}
