// Palmer Lab at UCSD 2026
//
#include <grm.h>


const char* grm::status_msg(grm::Status status) {
    switch (status) {
    case grm::Status::Success:
        return "Success";
    case grm::Status::ErrDimensionsNotEqual:
        return "Dimension(s) of data do not match grm.";
    case grm::Status::ErrIndexOutofBounds:
        return "Index out of bounds";
    default:
        break;
    }

    return "Unexpected status, please contact maintainers.";
}

