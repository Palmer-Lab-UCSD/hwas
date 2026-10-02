
#include <Rmatutils.h>


bool matutils::is_symm_(const Rcpp::NumericMatrix& A) {
    double tol = 1e-10;

    int nrow = A.nrow();
    if (nrow != A.ncol())
        return false;

    double upper_bound = 0;
    double lower_bound = 0;

    for (int i = 0; i < nrow; i++) {
        for (int j = i+1; j < nrow; j++) {

            lower_bound = A(j, i) - tol;
            upper_bound = A(j, i) + tol;

            if (A(i, j) < lower_bound || A(i, j) > upper_bound)
                return false;
        }
    }

    return true;
}

Status matutils::matrix_sum(const Rcpp::NumericMatrix& A, 
        double& val) {

    double s = 0;
    for (int i = 0; i < A.nrow(); i++) {
        for (int j = 0; j < A.ncol(); j++)
            s += A(i, j);
    }

    val = s;

    return Status::Success;
}


Status matutils::matrix_trace(const Rcpp::NumericMatrix& A, 
        double& val) {

    if (A.nrow() != A.ncol())
        return Status::ErrNotSquareMatrix;

    double mtrace = 0;
    for (int i = 0; i < A.nrow(); i++)
        mtrace += A(i, i);

    val = mtrace;
    return Status::Success;
}
