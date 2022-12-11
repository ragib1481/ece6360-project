//
// Created by ragib1481 on 12/6/22.
//

#ifndef PROJECT_MYMATH_CUH
#define PROJECT_MYMATH_CUH

#include <math.h>
#include "helper.cuh"

#define el 0.5772156649015329

__host__ __device__
int factorial(int n) {
    int val = 1;
    for (int i = 2; i <= n; i++) {
        val *= i;
    }
    return val;
}

namespace redefined {
    template <typename P>
    __host__ __device__ int msta1(P x,int mp)
    {
        P a0,f0,f1,f;
        int i,n0,n1,nn;

        a0 = fabs(x);
        n0 = (int)(1.1*a0)+1;
        f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-mp;
        n1 = n0+5;
        f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-mp;
        for (i=0;i<20;i++) {
            nn = (int)(n1-(n1-n0)/(1.0-f0/f1));
            f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-mp;
            if (std::abs(nn-n1) < 1) break;
            n0 = n1;
            f0 = f1;
            n1 = nn;
            f1 = f;
        }
        return nn;
    }

    template <typename P>
    __host__ __device__ int msta2(P x,int n,int mp)
    {
        P a0,ejn,hmp,f0,f1,f,obj;
        int i,n0,n1,nn;

        a0 = fabs(x);
        hmp = 0.5*mp;
        ejn = 0.5*log10(6.28*n)-n*log10(1.36*a0/n);
        if (ejn <= hmp) {
            obj = mp;
            n0 = (int)(1.1*a0);
            if (n0 < 1) n0 = 1;
        }
        else {
            obj = hmp+ejn;
            n0 = n;
        }
        f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-obj;
        n1 = n0+5;
        f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-obj;
        for (i=0;i<20;i++) {
            nn = (int)(n1-(n1-n0)/(1.0-f0/f1));
            f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-obj;
            if (std::abs(nn-n1) < 1) break;
            n0 = n1;
            f0 = f1;
            n1 = nn;
            f1 = f;
        }
        return nn+10;
    }

    template <typename P >
    __host__ __device__ P gamma(P x)
    {
        const P EPS = std::numeric_limits<P>::epsilon();
        const P FPMIN_MAG = std::numeric_limits<P>::min();
        const P FPMIN = std::numeric_limits<P>::lowest();
        const P FPMAX = std::numeric_limits<P>::max();

        int i,k,m;
        P ga,gr,r,z;

        static P g[] = {
                1.0,
                0.5772156649015329,
                -0.6558780715202538,
                -0.420026350340952e-1,
                0.1665386113822915,
                -0.421977345555443e-1,
                -0.9621971527877e-2,
                0.7218943246663e-2,
                -0.11651675918591e-2,
                -0.2152416741149e-3,
                0.1280502823882e-3,
                -0.201348547807e-4,
                -0.12504934821e-5,
                0.1133027232e-5,
                -0.2056338417e-6,
                0.6116095e-8,
                0.50020075e-8,
                -0.11812746e-8,
                0.1043427e-9,
                0.77823e-11,
                -0.36968e-11,
                0.51e-12,
                -0.206e-13,
                -0.54e-14,
                0.14e-14};

        if (x > 171.0) return FPMAX;    // This value is an overflow flag.
        if (x == (int)x) {
            if (x > 0.0) {
                ga = 1.0;               // use factorial
                for (i=2;i<x;i++) {
                    ga *= i;
                }
            }
            else
                ga = FPMAX;
        }
        else {
            if (fabs(x) > 1.0) {
                z = fabs(x);
                m = (int)z;
                r = 1.0;
                for (k=1;k<=m;k++) {
                    r *= (z-k);
                }
                z -= m;
            }
            else
                z = x;
            gr = g[24];
            for (k=23;k>=0;k--) {
                gr = gr*z+g[k];
            }
            ga = 1.0/(gr*z);
            if (fabs(x) > 1.0) {
                ga *= r;
                if (x < 0.0) {
                    ga = -M_PI/(x*ga*sin(M_PI*x));
                }
            }
        }
        return ga;
    }

    template <typename P>
    __host__ __device__ int bessjyv(P v,P x,P &vm,P *jv,P *yv, P *djv,P *dyv){
        P v0, vl, vg;
        P vv,a,a0,r,x2,bjv0,bjv1,bjvl,f,f0,f1,f2;
        P r0,r1,ck,cs,cs0,cs1,sk,qx,px,byv0,byv1,rp,xk,rq;
        P b,ec,w0,w1,bju0,bju1,pv0,pv1,byvk;
        int j,k,l,m,n,kz;

        const P EPS = std::numeric_limits<P>::epsilon();
        const P FPMIN_MAG = std::numeric_limits<P>::min();
        const P FPMIN = std::numeric_limits<P>::lowest();
        const P FPMAX = std::numeric_limits<P>::max();

        x2 = x*x;
        n = (int)v;
        v0 = v-n;
        if ((x < 0.0) || (v < 0.0)) return 1;
        if (x < EPS) {
            for (k=0;k<=n;k++) {
                jv[k] = 0.0;
                yv[k] = FPMIN;
                djv[k] = 0.0;
                dyv[k] = FPMAX;
                if (v0 == 0.0) {
                    jv[0] = 1.0;
                    djv[1] = 0.5;
                }
                else djv[0] = FPMAX;
            }
            vm = v;
            return 0;
        }
        if (x <= 12.0) {
            for (l=0;l<2;l++) {
                vl = v0 + l;
                bjvl = 1.0;
                r = 1.0;
                for (k=1;k<=40;k++) {
                    r *= -0.25*x2/(k*(k+vl));
                    bjvl += r;
                    if (fabs(r) < fabs(bjvl)*EPS) break;
                }
                vg = 1.0 + vl;
                a = pow(0.5*x,vl)/gamma(vg);
                if (l == 0) bjv0 = bjvl*a;
                else bjv1 = bjvl*a;
            }
        }
        else {
            if (x >= 50.0) kz = 8;
            else if (x >= 35.0) kz = 10;
            else kz = 11;
            for (j=0;j<2;j++) {
                vv = 4.0*(j+v0)*(j+v0);
                px = 1.0;
                rp = 1.0;
                for (k=1;k<=kz;k++) {
                    rp *= (-0.78125e-2)*(vv-pow(4.0*k-3.0,2.0))*
                          (vv-pow(4.0*k-1.0,2.0))/(k*(2.0*k-1.0)*x2);
                    px += rp;
                }
                qx = 1.0;
                rq = 1.0;
                for (k=1;k<=kz;k++) {
                    rq *= (-0.78125e-2)*(vv-pow(4.0*k-1.0,2.0))*
                          (vv-pow(4.0*k+1.0,2.0))/(k*(2.0*k+1.0)*x2);
                    qx += rq;
                }
                qx *= 0.125*(vv-1.0)/x;
                xk = x-(0.5*(j+v0)+0.25)*M_PI;
                a0 = sqrt(M_2_PI/x);
                ck = cos(xk);
                sk = sin(xk);

                if (j == 0) {
                    bjv0 = a0*(px*ck-qx*sk);
                    byv0 = a0*(px*sk+qx*ck);
                }
                else if (j == 1) {
                    bjv1 = a0*(px*ck-qx*sk);
                    byv1 = a0*(px*sk+qx*ck);
                }
            }
        }
        jv[0] = bjv0;
        jv[1] = bjv1;
        djv[0] = v0*jv[0]/x-jv[1];
        djv[1] = -(1.0+v0)*jv[1]/x+jv[0];
        if ((n >= 2) && (n <= (int)(0.9*x))) {
            f0 = bjv0;
            f1 = bjv1;
            for (k=2;k<=n;k++) {
                f = 2.0*(k+v0-1.0)*f1/x-f0;
                jv[k] = f;
                f0 = f1;
                f1 = f;
            }
        }
        else if (n >= 2) {
            m = msta1(x,200);
            if (m < n) n = m;
            else m = msta2(x,n,15);
            f2 = 0.0;
            f1 = FPMIN_MAG;
            for (k=m;k>=0;k--) {
                f = 2.0*(v0+k+1.0)*f1/x-f2;
                if (k <= n) jv[k] = f;
                f2 = f1;
                f1 = f;
            }
            if (fabs(bjv0) > fabs(bjv1)) cs = bjv0/f;
            else cs = bjv1/f2;
            for (k=0;k<=n;k++) {
                jv[k] *= cs;
            }
        }
        for (k=2;k<=n;k++) {
            djv[k] = -(k+v0)*jv[k]/x+jv[k-1];
        }
        if (x <= 12.0) {
            if (v0 != 0.0) {
                for (l=0;l<2;l++) {
                    vl = v0 +l;
                    bjvl = 1.0;
                    r = 1.0;
                    for (k=1;k<=40;k++) {
                        r *= -0.25*x2/(k*(k-vl));
                        bjvl += r;
                        if (fabs(r) < fabs(bjvl)*1e-15) break;
                    }
                    vg = 1.0-vl;
                    b = pow(2.0/x,vl)/gamma(vg);
                    if (l == 0) bju0 = bjvl*b;
                    else bju1 = bjvl*b;
                }
                pv0 = M_PI*v0;
                pv1 = M_PI*(1.0+v0);
                byv0 = (bjv0*cos(pv0)-bju0)/sin(pv0);
                byv1 = (bjv1*cos(pv1)-bju1)/sin(pv1);
            }
            else {
                ec = log(0.5*x)+el;
                cs0 = 0.0;
                w0 = 0.0;
                r0 = 1.0;
                for (k=1;k<=30;k++) {
                    w0 += 1.0/k;
                    r0 *= -0.25*x2/(k*k);
                    cs0 += r0*w0;
                }
                byv0 = M_2_PI*(ec*bjv0-cs0);
                cs1 = 1.0;
                w1 = 0.0;
                r1 = 1.0;
                for (k=1;k<=30;k++) {
                    w1 += 1.0/k;
                    r1 *= -0.25*x2/(k*(k+1));
                    cs1 += r1*(2.0*w1+1.0/(k+1.0));
                }
                byv1 = M_2_PI*(ec*bjv1-1.0/x-0.25*x*cs1);
            }
        }
        yv[0] = byv0;
        yv[1] = byv1;
        for (k=2;k<=n;k++) {
            byvk = 2.0*(v0+k-1.0)*byv1/x-byv0;
            yv[k] = byvk;
            byv0 = byv1;
            byv1 = byvk;
        }
        dyv[0] = v0*yv[0]/x-yv[1];
        for (k=1;k<=n;k++) {
            dyv[k] = -(k+v0)*yv[k]/x+yv[k-1];
        }
        vm = n + v0;
        return 0;
    }

    template <typename P>
    __host__ __device__ int sphBessjyv(const int v, P z, P &vm, P* cjv, P* cyv, P* cjvp, P* cyvp){

        //first, compute the bessel functions of fractional order
        bessjyv<P>(v + (P)0.5, z, vm, cjv, cyv, cjvp, cyvp);

        if(z == 0){													//handle degenerate case of z = 0
            memset(cjv, 0, sizeof(P) * (v+1));
            cjv[0] = 1;
        }

        //iterate through each and scale
        for(int n = 0; n<=v; n++){

            if(z != 0){												//handle degenerate case of z = 0
                cjv[n] = cjv[n] * sqrt(M_PI/(z * 2.0));
                cyv[n] = cyv[n] * sqrt(M_PI/(z * 2.0));
            }

            cjvp[n] = -1.0 / (z * 2.0) * cjv[n] + cjvp[n] * sqrt(M_PI / (z * 2.0));
            cyvp[n] = -1.0 / (z * 2.0) * cyv[n] + cyvp[n] * sqrt(M_PI / (z * 2.0));
        }

        return 0;

    }

    template <typename T, typename P>
    __host__ __device__ int spHankel1(const int v, P x, T* hv, T* dhv) {
        P* jv  = new P[v+1];
        P* yv  = new P[v+1];
        P* djv = new P[v+1];
        P* dyv = new P[v+1];
        P vm;

        sphBessjyv<P>(v, x, vm, jv, yv, djv, dyv);
        for (int i = 0; i <= v; i++) {
            hv[i]  = T(jv[i], yv[i]) ;
            dhv[i] = T(djv[i], dyv[i]) ;
        }
        delete[] jv;
        delete[] djv;
        delete[] yv;
        delete[] dyv;

        return 0;
    }

    template <typename T>
    __host__ __device__ T legendre(const int n, const T& x) {
        T pn;
        T pn_2 = x;
        T pn_1 = (static_cast<T>(3.0) * x * x - static_cast<T>(1.0)) * static_cast<T>(0.5);

        if (n == 0) return static_cast<T>(1.0);
        if (n == 1) return pn_2;
        if (n == 2) return pn_1;
        if (n == -1) return static_cast<T>(1.0);

        for (int l = 3; l <= n; l++) {
            pn = ((static_cast<T>(2.0) * static_cast<T>(l) - static_cast<T>(1.0)) * x * pn_1 - (static_cast<T>(l) - static_cast<T>(1.0)) * pn_2) / static_cast<T>(l);
            pn_2 = pn_1;
            pn_1 = pn;
        }

        return pn;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spBesselJComplex(int n, thrust::complex<T> x) {
        float temp;
        thrust::complex<float> val(0, 0);
        for (int s = 0; s < 20; s++) {
            temp = pow<T>(-1.0, (T)s) / ((T)factorial(s) * tgamma((T)(s + n) + 1.5));
            val += pow<T>(x/thrust::complex<T>(2.0f, 0.0f), 2.0*s+n+0.5) * thrust::complex<T>(temp, 0);
        }

        return val * sqrt(thrust::complex<T>(M_PI / 2.0, 0)/ x);
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spBesselYComplex(int n, thrust::complex<T> x) {
        thrust::complex<T> mult(1, 0);
        if ((n % 2) == 0)  mult.real(-1.0);
        return mult * spBesselJComplex<T>(-n-1, x);
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spHankel1Complex(int n, thrust::complex<T> x) {
        return spBesselJComplex<T>(n, x) + thrust::complex<T>(0.0, 1.0) * spBesselYComplex<T>(n, x);
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spBesselJPComplex(int n, thrust::complex<T> x) {
        thrust::complex<T> val = thrust::complex<T>(n, 0) * spBesselJComplex<T>(n-1, x);
        val -= thrust::complex<T>(n+1, 0) * spBesselJComplex<T>(n+1, x);
        val /= thrust::complex<T>(2*n+1, 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spBesselYPComplex(int n, thrust::complex<T> x) {
        thrust::complex<T> val = thrust::complex<T>(n, 0) * spBesselYComplex<T>(n-1, x);
        val -= thrust::complex<T>(n+1, 0) * spBesselYComplex<T>(n+1, x);
        val /= thrust::complex<T>(2*n+1, 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spHankel1PComplex(int n, thrust::complex<T> x) {
        thrust::complex<T> val = thrust::complex<T>(n, 0) * spHankel1Complex<T>(n-1, x);
        val -= thrust::complex<T>(n+1, 0) * spHankel1Complex<T>(n+1, x);
        val /= thrust::complex<T>(2*n+1, 0);
        return val;
    }

}


#endif //PROJECT_MYMATH_CUH
