//
// Created by tyler on 6/25/17.
//

#ifndef BLAZE_ITERATIVE_PRECONDITIONBICGSTAB_HPP
#define BLAZE_ITERATIVE_PRECONDITIONBICGSTAB_HPP

#include "PreconditionBiCGSTABTag.hpp"

BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

namespace detail {

template<typename MatrixType, typename T>
void decomposition(std::string type, const MatrixType &A, MatrixType &K1, MatrixType &K2){
    
    
    if (type.compare("Cholesky")){
        
        MatrixType L;
        llh( A, L );  // LLH decomposition of a row-major matrix
        K1 = L;
        K2 = ctrans(L); 
    }
    
    if (type.compare("QR")){
        
        qr( A, K1, K2 );  //QR decomposition of a row-major matrix
        
    }
    
    if (type.compare("RQ")){
        
        rq( A, K1, K2 ); //RQ decomposition of a row-major matrix
    }
    
    if (type.compare("LU")){
        
    MatrixType P;
    lu( A, K1, K2, P ); 
    K2 *= P;
        
    }
    
     
}    
    
    
    
/**
 *  Implementation of the Preconditioned BiCGSTAB method, following the
 *  preconditioned version on Wikipedia using a LU Decomposition.
 */
template<typename MatrixType, typename T>
void solve_impl(
        DynamicVector<T> &x,
        const MatrixType &A,
        const DynamicVector<T> &b,
        PreconditionBiCGSTABTag &tag)
{

    // Decomposition A = K1 * K2
    MatrixType K1;
    MatrixType K2;
    decomposition<MatrixType,T>("LU",A,K1,K2);
    
    // Compute inverse
    auto Kinv = inv(K1*K2);
    auto K1inv = inv(K1 ); 
    
    DynamicVector<T> r = b - A * x;
    DynamicVector<T> p(r);
    DynamicVector<T> v(r);
    DynamicVector<T> y(r);
    DynamicVector<T> r0(r);
    DynamicVector<T> s(p.size());
    DynamicVector<T> t(p.size());
    DynamicVector<T> z(p.size());
    DynamicVector<T> error(r);

    auto absolute_residual_0 = trans(r) * r;
    auto absolute_residual = absolute_residual_0;

    auto rho_prev = T(1);
    auto w = T(1);
    auto alpha = T(1);


    std::size_t iteration{0};
    while (true) {

        auto rho = trans(r0) * r;
        auto beta = (rho * alpha) / (rho_prev * w);

        p = r + beta * (p - w * v);
        //v = A * p;
        y = Kinv * p;
        v = A * y;
        
        alpha = rho / (trans(r0) * v);
        
        auto h = x + alpha * y;
        error = b - A*h;
        absolute_residual = trans(error)*error;
        auto relative_residual = absolute_residual/absolute_residual_0;
         if(tag.do_log()) {
            tag.log_residual(relative_residual);
        }

        if(tag.terminateIteration(iteration,absolute_residual,relative_residual)) {
            x = h;
            break;
        }
        
        s = r - alpha * v;
        z = Kinv * s;
        t = A*z;

        // sometimes, t will be zero, so trans(t)*t is zero.
        // This happens if the solution is exactly correct,
        // So best to set w=0, and loop will terminate below.
        auto t_dot_t = blaze::dot(K1inv*t,K1inv*t);
        if(t_dot_t == 0)
            w = 0;
        else
            w = blaze::dot(K1inv*t,K1inv*s)/t_dot_t; // replace with 12


        x += alpha*y + w*z;

        error = b - A*x;
        absolute_residual = trans(error)*error;
        relative_residual = absolute_residual/absolute_residual_0;
        if(tag.do_log()) {
            tag.log_residual(relative_residual);
        }

        if(tag.terminateIteration(iteration,absolute_residual,relative_residual)) {
            break;
        }

        r = s - w*t;
        rho_prev = rho;

        ++iteration;
    }
}

} //end namespace detail

ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_BICGSTAB_HPP
