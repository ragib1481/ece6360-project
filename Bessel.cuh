//
// Created by ragib1481 on 12/12/22.
//

#ifndef PROJECT_BESSEL_CUH
#define PROJECT_BESSEL_CUH

#define _USE_MATH_DEFINES
#include <math.h>
#define eps 1e-15
#define el 0.5772156649015329

namespace stim{
    double PI = M_PI;
    // thrust::complex<double> cii(0.0,1.0);
    // thrust::complex<double> cone(1.0,0.0);
    // thrust::complex<double> czero(0.0,0.0);

    template< typename P >
    __host__ __device__
    P gamma(P x)
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

    template<typename P>
    __host__ __device__
    int bessjy01a(P x,P &j0,P &j1,P &y0,P &y1,
                  P &j0p,P &j1p,P &y0p,P &y1p)
    {
        const P EPS = std::numeric_limits<P>::epsilon();
        const P FPMIN_MAG = std::numeric_limits<P>::min();
        const P FPMIN = std::numeric_limits<P>::lowest();
        const P FPMAX = std::numeric_limits<P>::max();

        P x2,r,ec,w0,w1,r0,r1,cs0,cs1;
        P cu,p0,q0,p1,q1,t1,t2;
        int k,kz;
        static P a[] = {
                -7.03125e-2,
                0.112152099609375,
                -0.5725014209747314,
                6.074042001273483,
                -1.100171402692467e2,
                3.038090510922384e3,
                -1.188384262567832e5,
                6.252951493434797e6,
                -4.259392165047669e8,
                3.646840080706556e10,
                -3.833534661393944e12,
                4.854014686852901e14,
                -7.286857349377656e16,
                1.279721941975975e19};
        static P b[] = {
                7.32421875e-2,
                -0.2271080017089844,
                1.727727502584457,
                -2.438052969955606e1,
                5.513358961220206e2,
                -1.825775547429318e4,
                8.328593040162893e5,
                -5.006958953198893e7,
                3.836255180230433e9,
                -3.649010818849833e11,
                4.218971570284096e13,
                -5.827244631566907e15,
                9.476288099260110e17,
                -1.792162323051699e20};
        static P a1[] = {
                0.1171875,
                -0.1441955566406250,
                0.6765925884246826,
                -6.883914268109947,
                1.215978918765359e2,
                -3.302272294480852e3,
                1.276412726461746e5,
                -6.656367718817688e6,
                4.502786003050393e8,
                -3.833857520742790e10,
                4.011838599133198e12,
                -5.060568503314727e14,
                7.572616461117958e16,
                -1.326257285320556e19};
        static P b1[] = {
                -0.1025390625,
                0.2775764465332031,
                -1.993531733751297,
                2.724882731126854e1,
                -6.038440767050702e2,
                1.971837591223663e4,
                -8.902978767070678e5,
                5.310411010968522e7,
                -4.043620325107754e9,
                3.827011346598605e11,
                -4.406481417852278e13,
                6.065091351222699e15,
                -9.833883876590679e17,
                1.855045211579828e20};

        if (x < 0.0) return 1;
        if (x == 0.0) {
            j0 = 1.0;
            j1 = 0.0;
            y0 = -FPMIN;
            y1 = -FPMIN;
            j0p = 0.0;
            j1p = 0.5;
            y0p = FPMAX;
            y1p = FPMAX;
            return 0;
        }
        x2 = x*x;
        if (x <= 12.0) {
            j0 = 1.0;
            r = 1.0;
            for (k=1;k<=30;k++) {
                r *= -0.25*x2/(k*k);
                j0 += r;
                if (fabs(r) < fabs(j0)*1e-15) break;
            }
            j1 = 1.0;
            r = 1.0;
            for (k=1;k<=30;k++) {
                r *= -0.25*x2/(k*(k+1));
                j1 += r;
                if (fabs(r) < fabs(j1)*1e-15) break;
            }
            j1 *= 0.5*x;
            ec = log(0.5*x)+el;
            cs0 = 0.0;
            w0 = 0.0;
            r0 = 1.0;
            for (k=1;k<=30;k++) {
                w0 += 1.0/k;
                r0 *= -0.25*x2/(k*k);
                r = r0 * w0;
                cs0 += r;
                if (fabs(r) < fabs(cs0)*1e-15) break;
            }
            y0 = M_2_PI*(ec*j0-cs0);
            cs1 = 1.0;
            w1 = 0.0;
            r1 = 1.0;
            for (k=1;k<=30;k++) {
                w1 += 1.0/k;
                r1 *= -0.25*x2/(k*(k+1));
                r = r1*(2.0*w1+1.0/(k+1));
                cs1 += r;
                if (fabs(r) < fabs(cs1)*1e-15) break;
            }
            y1 = M_2_PI * (ec*j1-1.0/x-0.25*x*cs1);
        }
        else {
            if (x >= 50.0) kz = 8;          // Can be changed to 10
            else if (x >= 35.0) kz = 10;    //  "       "        12
            else kz = 12;                   //  "       "        14
            t1 = x-M_PI_4;
            p0 = 1.0;
            q0 = -0.125/x;
            for (k=0;k<kz;k++) {
                p0 += a[k]*pow(x,-2*k-2);
                q0 += b[k]*pow(x,-2*k-3);
            }
            cu = sqrt(M_2_PI/x);
            j0 = cu*(p0*cos(t1)-q0*sin(t1));
            y0 = cu*(p0*sin(t1)+q0*cos(t1));
            t2 = x-0.75*M_PI;
            p1 = 1.0;
            q1 = 0.375/x;
            for (k=0;k<kz;k++) {
                p1 += a1[k]*pow(x,-2*k-2);
                q1 += b1[k]*pow(x,-2*k-3);
            }
            j1 = cu*(p1*cos(t2)-q1*sin(t2));
            y1 = cu*(p1*sin(t2)+q1*cos(t2));
        }
        j0p = -j1;
        j1p = j0-j1/x;
        y0p = -y1;
        y1p = y0-y1/x;
        return 0;
    }
//
//  INPUT:
//      double x    -- argument of Bessel function
//
//  OUTPUT:
//      double j0   -- Bessel function of 1st kind, 0th order
//      double j1   -- Bessel function of 1st kind, 1st order
//      double y0   -- Bessel function of 2nd kind, 0th order
//      double y1   -- Bessel function of 2nd kind, 1st order
//      double j0p  -- derivative of Bessel function of 1st kind, 0th order
//      double j1p  -- derivative of Bessel function of 1st kind, 1st order
//      double y0p  -- derivative of Bessel function of 2nd kind, 0th order
//      double y1p  -- derivative of Bessel function of 2nd kind, 1st order
//
//  RETURN:
//      int error code: 0 = OK, 1 = error
//
//  This algorithm computes the functions using polynomial approximations.
//
    template<typename P>
    __host__ __device__
    int bessjy01b(P x,P &j0,P &j1,P &y0,P &y1,
                  P &j0p,P &j1p,P &y0p,P &y1p)
    {
        P t,t2,dtmp,a0,p0,q0,p1,q1,ta0,ta1;
        if (x < 0.0) return 1;
        if (x == 0.0) {
            j0 = 1.0;
            j1 = 0.0;
            y0 = -1e308;
            y1 = -1e308;
            j0p = 0.0;
            j1p = 0.5;
            y0p = 1e308;
            y1p = 1e308;
            return 0;
        }
        if(x <= 4.0) {
            t = x/4.0;
            t2 = t*t;
            j0 = ((((((-0.5014415e-3*t2+0.76771853e-2)*t2-0.0709253492)*t2+
                     0.4443584263)*t2-1.7777560599)*t2+3.9999973021)*t2
                  -3.9999998721)*t2+1.0;
            j1 = t*(((((((-0.1289769e-3*t2+0.22069155e-2)*t2-0.0236616773)*t2+
                        0.1777582922)*t2-0.8888839649)*t2+2.6666660544)*t2-
                     3.999999971)*t2+1.9999999998);
            dtmp = (((((((-0.567433e-4*t2+0.859977e-3)*t2-0.94855882e-2)*t2+
                        0.0772975809)*t2-0.4261737419)*t2+1.4216421221)*t2-
                     2.3498519931)*t2+1.0766115157)*t2+0.3674669052;
            y0 = M_2_PI*log(0.5*x)*j0+dtmp;
            dtmp = (((((((0.6535773e-3*t2-0.0108175626)*t2+0.107657607)*t2-
                        0.7268945577)*t2+3.1261399273)*t2-7.3980241381)*t2+
                     6.8529236342)*t2+0.3932562018)*t2-0.6366197726;
            y1 = M_2_PI*log(0.5*x)*j1+dtmp/x;
        }
        else {
            t = 4.0/x;
            t2 = t*t;
            a0 = sqrt(M_2_PI/x);
            p0 = ((((-0.9285e-5*t2+0.43506e-4)*t2-0.122226e-3)*t2+
                   0.434725e-3)*t2-0.4394275e-2)*t2+0.999999997;
            q0 = t*(((((0.8099e-5*t2-0.35614e-4)*t2+0.85844e-4)*t2-
                      0.218024e-3)*t2+0.1144106e-2)*t2-0.031249995);
            ta0 = x-M_PI_4;
            j0 = a0*(p0*cos(ta0)-q0*sin(ta0));
            y0 = a0*(p0*sin(ta0)+q0*cos(ta0));
            p1 = ((((0.10632e-4*t2-0.50363e-4)*t2+0.145575e-3)*t2
                   -0.559487e-3)*t2+0.7323931e-2)*t2+1.000000004;
            q1 = t*(((((-0.9173e-5*t2+0.40658e-4)*t2-0.99941e-4)*t2
                      +0.266891e-3)*t2-0.1601836e-2)*t2+0.093749994);
            ta1 = x-0.75*M_PI;
            j1 = a0*(p1*cos(ta1)-q1*sin(ta1));
            y1 = a0*(p1*sin(ta1)+q1*cos(ta1));
        }
        j0p = -j1;
        j1p = j0-j1/x;
        y0p = -y1;
        y1p = y0-y1/x;
        return 0;
    }

    template<typename P>
    __host__ __device__
    int msta1(P x,int mp)
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
            if (abs(nn-n1) < 1) break;
            n0 = n1;
            f0 = f1;
            n1 = nn;
            f1 = f;
        }
        return nn;
    }

    template<typename P>
    __host__ __device__
    int msta2(P x,int n,int mp)
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
            if (abs(nn-n1) < 1) break;
            n0 = n1;
            f0 = f1;
            n1 = nn;
            f1 = f;
        }
        return nn+10;
    }
//
//  INPUT:
//  double x    -- argument of Bessel function of 1st and 2nd kind.
//  int n       -- order
//
//  OUPUT:
//
//  int nm      -- highest order actually computed (nm <= n)
//  double jn[] -- Bessel function of 1st kind, orders from 0 to nm
//  double yn[] -- Bessel function of 2nd kind, orders from 0 to nm
//  double j'n[]-- derivative of Bessel function of 1st kind,
//                      orders from 0 to nm
//  double y'n[]-- derivative of Bessel function of 2nd kind,
//                      orders from 0 to nm
//
//  Computes Bessel functions of all order up to 'n' using recurrence
//  relations. If 'nm' < 'n' only 'nm' orders are returned.
//
    template<typename P>
    __host__ __device__
    int bessjyna(int n,P x,int &nm,P *jn,P *yn,
                 P *jnp,P *ynp)
    {
        P bj0,bj1,f,f0,f1,f2,cs;
        int i,k,m,ecode;

        nm = n;
        if ((x < 0.0) || (n < 0)) return 1;
        if (x < 1e-15) {
            for (i=0;i<=n;i++) {
                jn[i] = 0.0;
                yn[i] = -1e308;
                jnp[i] = 0.0;
                ynp[i] = 1e308;
            }
            jn[0] = 1.0;
            jnp[1] = 0.5;
            return 0;
        }
        ecode = bessjy01a(x,jn[0],jn[1],yn[0],yn[1],jnp[0],jnp[1],ynp[0],ynp[1]);
        if (n < 2) return 0;
        bj0 = jn[0];
        bj1 = jn[1];
        if (n < (int)0.9*x) {
            for (k=2;k<=n;k++) {
                jn[k] = 2.0*(k-1.0)*bj1/x-bj0;
                bj0 = bj1;
                bj1 = jn[k];
            }
        }
        else {
            m = msta1(x,200);
            if (m < n) nm = m;
            else m = msta2(x,n,15);
            f2 = 0.0;
            f1 = 1.0e-100;
            for (k=m;k>=0;k--) {
                f = 2.0*(k+1.0)/x*f1-f2;
                if (k <= nm) jn[k] = f;
                f2 = f1;
                f1 = f;
            }
            if (fabs(bj0) > fabs(bj1)) cs = bj0/f;
            else cs = bj1/f2;
            for (k=0;k<=nm;k++) {
                jn[k] *= cs;
            }
        }
        for (k=2;k<=nm;k++) {
            jnp[k] = jn[k-1]-k*jn[k]/x;
        }
        f0 = yn[0];
        f1 = yn[1];
        for (k=2;k<=nm;k++) {
            f = 2.0*(k-1.0)*f1/x-f0;
            yn[k] = f;
            f0 = f1;
            f1 = f;
        }
        for (k=2;k<=nm;k++) {
            ynp[k] = yn[k-1]-k*yn[k]/x;
        }
        return 0;
    }
//
//  Same input and output conventions as above. Different recurrence
//  relations used for 'x' < 300.
//
    template<typename P>
    __host__ __device__
    int bessjynb(int n,P x,int &nm,P *jn,P *yn,
                 P *jnp,P *ynp)
    {
        P t1,t2,f,f1,f2,bj0,bj1,bjk,by0,by1,cu,s0,su,sv;
        P ec,bs,byk,p0,p1,q0,q1;
        static P a[] = {
                -0.7031250000000000e-1,
                0.1121520996093750,
                -0.5725014209747314,
                6.074042001273483};
        static P b[] = {
                0.7324218750000000e-1,
                -0.2271080017089844,
                1.727727502584457,
                -2.438052969955606e1};
        static P a1[] = {
                0.1171875,
                -0.1441955566406250,
                0.6765925884246826,
                -6.883914268109947};
        static P b1[] = {
                -0.1025390625,
                0.2775764465332031,
                -1.993531733751297,
                2.724882731126854e1};

        int i,k,m;
        nm = n;
        if ((x < 0.0) || (n < 0)) return 1;
        if (x < 1e-15) {
            for (i=0;i<=n;i++) {
                jn[i] = 0.0;
                yn[i] = -1e308;
                jnp[i] = 0.0;
                ynp[i] = 1e308;
            }
            jn[0] = 1.0;
            jnp[1] = 0.5;
            return 0;
        }
        if (x <= 300.0 || n > (int)(0.9*x)) {
            if (n == 0) nm = 1;
            m = msta1(x,200);
            if (m < nm) nm = m;
            else m = msta2(x,nm,15);
            bs = 0.0;
            su = 0.0;
            sv = 0.0;
            f2 = 0.0;
            f1 = 1.0e-100;
            for (k = m;k>=0;k--) {
                f = 2.0*(k+1.0)/x*f1 - f2;
                if (k <= nm) jn[k] = f;
                if ((k == 2*(int)(k/2)) && (k != 0)) {
                    bs += 2.0*f;
//                su += pow(-1,k>>1)*f/(double)k;
                    su += (-1)*((k & 2)-1)*f/(P)k;
                }
                else if (k > 1) {
//                sv += pow(-1,k>>1)*k*f/(k*k-1.0);
                    sv += (-1)*((k & 2)-1)*(P)k*f/(k*k-1.0);
                }
                f2 = f1;
                f1 = f;
            }
            s0 = bs+f;
            for (k=0;k<=nm;k++) {
                jn[k] /= s0;
            }
            ec = log(0.5*x) +0.5772156649015329;
            by0 = M_2_PI*(ec*jn[0]-4.0*su/s0);
            yn[0] = by0;
            by1 = M_2_PI*((ec-1.0)*jn[1]-jn[0]/x-4.0*sv/s0);
            yn[1] = by1;
        }
        else {
            t1 = x-M_PI_4;
            p0 = 1.0;
            q0 = -0.125/x;
            for (k=0;k<4;k++) {
                p0 += a[k]*pow(x,-2*k-2);
                q0 += b[k]*pow(x,-2*k-3);
            }
            cu = sqrt(M_2_PI/x);
            bj0 = cu*(p0*cos(t1)-q0*sin(t1));
            by0 = cu*(p0*sin(t1)+q0*cos(t1));
            jn[0] = bj0;
            yn[0] = by0;
            t2 = x-0.75*M_PI;
            p1 = 1.0;
            q1 = 0.375/x;
            for (k=0;k<4;k++) {
                p1 += a1[k]*pow(x,-2*k-2);
                q1 += b1[k]*pow(x,-2*k-3);
            }
            bj1 = cu*(p1*cos(t2)-q1*sin(t2));
            by1 = cu*(p1*sin(t2)+q1*cos(t2));
            jn[1] = bj1;
            yn[1] = by1;
            for (k=2;k<=nm;k++) {
                bjk = 2.0*(k-1.0)*bj1/x-bj0;
                jn[k] = bjk;
                bj0 = bj1;
                bj1 = bjk;
            }
        }
        jnp[0] = -jn[1];
        for (k=1;k<=nm;k++) {
            jnp[k] = jn[k-1]-k*jn[k]/x;
        }
        for (k=2;k<=nm;k++) {
            byk = 2.0*(k-1.0)*by1/x-by0;
            yn[k] = byk;
            by0 = by1;
            by1 = byk;
        }
        ynp[0] = -yn[1];
        for (k=1;k<=nm;k++) {
            ynp[k] = yn[k-1]-k*yn[k]/x;
        }
        return 0;

    }

//  The following routine computes Bessel Jv(x) and Yv(x) for
//  arbitrary positive order (v). For negative order, use:
//
//      J-v(x) = Jv(x)cos(v pi) - Yv(x)sin(v pi)
//      Y-v(x) = Jv(x)sin(v pi) + Yv(x)cos(v pi)
//
    template<typename P>
    __host__ __device__
    int bessjyv(P v,P x,P &vm,P *jv,P *yv,
                P *djv,P *dyv)
    {
        P v0,vl,vg,vv,a,a0,r,x2,bjv0,bjv1,bjvl,f,f0,f1,f2;
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

    template<typename P>
    __host__ __device__
    int bessjyv_sph(int v, P z, P &vm, P* cjv,
                    P* cyv, P* cjvp, P* cyvp){

        //first, compute the bessel functions of fractional order
        bessjyv<P>(v + (P)0.5, z, vm, cjv, cyv, cjvp, cyvp);

        if(z == 0){													//handle degenerate case of z = 0
            memset(cjv, 0, sizeof(P) * (v+1));
            cjv[0] = 1;
        }

        //iterate through each and scale
        for(int n = 0; n<=v; n++){

            if(z != 0){												//handle degenerate case of z = 0
                cjv[n] = cjv[n] * sqrt(stim::PI/(z * 2.0));
                cyv[n] = cyv[n] * sqrt(stim::PI/(z * 2.0));
            }

            cjvp[n] = -1.0 / (z * 2.0) * cjv[n] + cjvp[n] * sqrt(stim::PI / (z * 2.0));
            cyvp[n] = -1.0 / (z * 2.0) * cyv[n] + cyvp[n] * sqrt(stim::PI / (z * 2.0));
        }

        return 0;

    }

    template<typename P>
    __host__ __device__
    int cbessjy01(thrust::complex<P> z, thrust::complex<P> &cj0, thrust::complex<P> &cj1,
                  thrust::complex<P> &cy0, thrust::complex<P> &cy1, thrust::complex<P> &cj0p,
                  thrust::complex<P> &cj1p, thrust::complex<P> &cy0p, thrust::complex<P> &cy1p)
    {
        thrust::complex<P> z1,z2,cr,cp,cs,cp0,cq0,cp1,cq1,ct1,ct2,cu;
        thrust::complex<double> cii(0.0,1.0);
        thrust::complex<double> cone(1.0,0.0);
        thrust::complex<double> czero(0.0,0.0);
        P a0,w0,w1;
        int k,kz;

        static P a[] = {
                -7.03125e-2,
                0.112152099609375,
                -0.5725014209747314,
                6.074042001273483,
                -1.100171402692467e2,
                3.038090510922384e3,
                -1.188384262567832e5,
                6.252951493434797e6,
                -4.259392165047669e8,
                3.646840080706556e10,
                -3.833534661393944e12,
                4.854014686852901e14,
                -7.286857349377656e16,
                1.279721941975975e19};
        static P b[] = {
                7.32421875e-2,
                -0.2271080017089844,
                1.727727502584457,
                -2.438052969955606e1,
                5.513358961220206e2,
                -1.825775547429318e4,
                8.328593040162893e5,
                -5.006958953198893e7,
                3.836255180230433e9,
                -3.649010818849833e11,
                4.218971570284096e13,
                -5.827244631566907e15,
                9.476288099260110e17,
                -1.792162323051699e20};
        static P a1[] = {
                0.1171875,
                -0.1441955566406250,
                0.6765925884246826,
                -6.883914268109947,
                1.215978918765359e2,
                -3.302272294480852e3,
                1.276412726461746e5,
                -6.656367718817688e6,
                4.502786003050393e8,
                -3.833857520742790e10,
                4.011838599133198e12,
                -5.060568503314727e14,
                7.572616461117958e16,
                -1.326257285320556e19};
        static P b1[] = {
                -0.1025390625,
                0.2775764465332031,
                -1.993531733751297,
                2.724882731126854e1,
                -6.038440767050702e2,
                1.971837591223663e4,
                -8.902978767070678e5,
                5.310411010968522e7,
                -4.043620325107754e9,
                3.827011346598605e11,
                -4.406481417852278e13,
                6.065091351222699e15,
                -9.833883876590679e17,
                1.855045211579828e20};

        a0 = abs(z);
        z2 = z*z;
        z1 = z;
        if (a0 == 0.0) {
            cj0 = cone;
            cj1 = czero;
            cy0 = thrust::complex<P>(-1e308,0);
            cy1 = thrust::complex<P>(-1e308,0);
            cj0p = czero;
            cj1p = thrust::complex<P>(0.5,0.0);
            cy0p = thrust::complex<P>(1e308,0);
            cy1p = thrust::complex<P>(1e308,0);
            return 0;
        }
        if (real(z) < 0.0) z1 = -z;
        if (a0 <= 12.0) {
            cj0 = cone;
            cr = cone;
            for (k=1;k<=40;k++) {
                cr *= -0.25*z2/(P)(k*k);
                cj0 += cr;
                if (abs(cr) < abs(cj0)*eps) break;
            }
            cj1 = cone;
            cr = cone;
            for (k=1;k<=40;k++) {
                cr *= -0.25*z2/(k*(k+1.0));
                cj1 += cr;
                if (abs(cr) < abs(cj1)*eps) break;
            }
            cj1 *= 0.5*z1;
            w0 = 0.0;
            cr = cone;
            cs = czero;
            for (k=1;k<=40;k++) {
                w0 += 1.0/k;
                cr *= -0.25*z2/(P)(k*k);
                cp = cr*w0;
                cs += cp;
                if (abs(cp) < abs(cs)*eps) break;
            }
            cy0 = M_2_PI*((log(0.5*z1)+el)*cj0-cs);
            w1 = 0.0;
            cr = cone;
            cs = cone;
            for (k=1;k<=40;k++) {
                w1 += 1.0/k;
                cr *= -0.25*z2/(k*(k+1.0));
                cp = cr*(2.0*w1+1.0/(k+1.0));
                cs += cp;
                if (abs(cp) < abs(cs)*eps) break;
            }
            cy1 = M_2_PI*((log(0.5*z1)+el)*cj1-1.0/z1-0.25*z1*cs);
        }
        else {
            if (a0 >= 50.0) kz = 8;         // can be changed to 10
            else if (a0 >= 35.0) kz = 10;   //   "      "     "  12
            else kz = 12;                   //   "      "     "  14
            ct1 = z1 - M_PI_4;
            cp0 = cone;
            for (k=0;k<kz;k++) {
                cp0 += a[k]*pow(z1,-2.0*k-2.0);
            }
            cq0 = -0.125/z1;
            for (k=0;k<kz;k++) {
                cq0 += b[k]*pow(z1,-2.0*k-3.0);
            }
            cu = sqrt(M_2_PI/z1);
            cj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1));
            cy0 = cu*(cp0*sin(ct1)+cq0*cos(ct1));
            ct2 = z1 - 0.75*M_PI;
            cp1 = cone;
            for (k=0;k<kz;k++) {
                cp1 += a1[k]*pow(z1,-2.0*k-2.0);
            }
            cq1 = 0.375/z1;
            for (k=0;k<kz;k++) {
                cq1 += b1[k]*pow(z1,-2.0*k-3.0);
            }
            cj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2));
            cy1 = cu*(cp1*sin(ct2)+cq1*cos(ct2));
        }
        if (real(z) < 0.0) {
            if (imag(z) < 0.0) {
                cy0 -= 2.0*cii*cj0;
                cy1 = -(cy1-2.0*cii*cj1);
            }
            else if (imag(z) > 0.0) {
                cy0 += 2.0*cii*cj0;
                cy1 = -(cy1+2.0*cii*cj1);
            }
            cj1 = -cj1;
        }
        cj0p = -cj1;
        cj1p = cj0-cj1/z;
        cy0p = -cy1;
        cy1p = cy0-cy1/z;
        return 0;
    }

    template<typename P>
    __host__ __device__
    int cbessjyna(int n, thrust::complex<P> z, int &nm, thrust::complex<P> *cj,
                  thrust::complex<P> *cy, thrust::complex<P> *cjp, thrust::complex<P> *cyp)
    {
        thrust::complex<P> cbj0,cbj1,cby0,cby1,cj0,cjk,cj1,cf,cf1,cf2;
        thrust::complex<P> cs,cg0,cg1,cyk,cyl1,cyl2,cylk,cp11,cp12,cp21,cp22;
        thrust::complex<P> ch0,ch1,ch2;
        thrust::complex<double> cii(0.0,1.0);
        thrust::complex<double> cone(1.0,0.0);
        thrust::complex<double> czero(0.0,0.0);
        P a0,yak,ya1,ya0,wa;
        int m,k,lb,lb0;

        if (n < 0) return 1;
        a0 = abs(z);
        nm = n;
        if (a0 < 1.0e-100) {
            for (k=0;k<=n;k++) {
                cj[k] = czero;
                cy[k] = thrust::complex<P> (-1e308,0);
                cjp[k] = czero;
                cyp[k] = thrust::complex<P>(1e308,0);
            }
            cj[0] = cone;
            cjp[1] = thrust::complex<P>(0.5,0.0);
            return 0;
        }
        cbessjy01(z,cj[0],cj[1],cy[0],cy[1],cjp[0],cjp[1],cyp[0],cyp[1]);
        cbj0 = cj[0];
        cbj1 = cj[1];
        cby0 = cy[0];
        cby1 = cy[1];
        if (n <= 1) return 0;
        if (n < (int)0.25*a0) {
            cj0 = cbj0;
            cj1 = cbj1;
            for (k=2;k<=n;k++) {
                cjk = 2.0*(k-1.0)*cj1/z-cj0;
                cj[k] = cjk;
                cj0 = cj1;
                cj1 = cjk;
            }
        }
        else {
            m = msta1(a0,200);
            if (m < n) nm = m;
            else m = msta2(a0,n,15);
            cf2 = czero;
            cf1 = thrust::complex<P> (1.0e-20,0.0);
            for (k=m;k>=0;k--) {
                cf = 2.0*(k+1.0)*cf1/z-cf2;
                if (k <=nm) cj[k] = cf;
                cf2 = cf1;
                cf1 = cf;
            }
            if (abs(cbj0) > abs(cbj1)) cs = cbj0/cf;
            else cs = cbj1/cf2;
            for (k=0;k<=nm;k++) {
                cj[k] *= cs;
            }
        }
        for (k=2;k<=nm;k++) {
            cjp[k] = cj[k-1]-(P)k*cj[k]/z;
        }
        ya0 = abs(cby0);
        lb = 0;
        cg0 = cby0;
        cg1 = cby1;
        for (k=2;k<=nm;k++) {
            cyk = 2.0*(k-1.0)*cg1/z-cg0;
            yak = abs(cyk);
            ya1 = abs(cg0);
            if ((yak < ya0) && (yak < ya1)) lb = k;
            cy[k] = cyk;
            cg0 = cg1;
            cg1 = cyk;
        }
        lb0 = 0;
        if ((lb > 4) && (imag(z) != 0.0)) {
            while (lb != lb0) {
                ch2 = cone;
                ch1 = czero;
                lb0 = lb;
                for (k=lb;k>=1;k--) {
                    ch0 = 2.0*k*ch1/z-ch2;
                    ch2 = ch1;
                    ch1 = ch0;
                }
                cp12 = ch0;
                cp22 = ch2;
                ch2 = czero;
                ch1 = cone;
                for (k=lb;k>=1;k--) {
                    ch0 = 2.0*k*ch1/z-ch2;
                    ch2 = ch1;
                    ch1 = ch0;
                }
                cp11 = ch0;
                cp21 = ch2;
                if (lb == nm)
                    cj[lb+1] = 2.0*lb*cj[lb]/z-cj[lb-1];
                if (abs(cj[0]) > abs(cj[1])) {
                    cy[lb+1] = (cj[lb+1]*cby0-2.0*cp11/(M_PI*z))/cj[0];
                    cy[lb] = (cj[lb]*cby0+2.0*cp12/(M_PI*z))/cj[0];
                }
                else {
                    cy[lb+1] = (cj[lb+1]*cby1-2.0*cp21/(M_PI*z))/cj[1];
                    cy[lb] = (cj[lb]*cby1+2.0*cp22/(M_PI*z))/cj[1];
                }
                cyl2 = cy[lb+1];
                cyl1 = cy[lb];
                for (k=lb-1;k>=0;k--) {
                    cylk = 2.0*(k+1.0)*cyl1/z-cyl2;
                    cy[k] = cylk;
                    cyl2 = cyl1;
                    cyl1 = cylk;
                }
                cyl1 = cy[lb];
                cyl2 = cy[lb+1];
                for (k=lb+1;k<n;k++) {
                    cylk = 2.0*k*cyl2/z-cyl1;
                    cy[k+1] = cylk;
                    cyl1 = cyl2;
                    cyl2 = cylk;
                }
                for (k=2;k<=nm;k++) {
                    wa = abs(cy[k]);
                    if (wa < abs(cy[k-1])) lb = k;
                }
            }
        }
        for (k=2;k<=nm;k++) {
            cyp[k] = cy[k-1]-(P)k*cy[k]/z;
        }
        return 0;
    }

    template<typename P>
    __host__ __device__
    int cbessjynb(int n, thrust::complex<P> z, int &nm, thrust::complex<P> *cj,
                  thrust::complex<P> *cy, thrust::complex<P> *cjp, thrust::complex<P> *cyp)
    {
        thrust::complex<P> cf,cf0,cf1,cf2,cbs,csu,csv,cs0,ce;
        thrust::complex<P> ct1,cp0,cq0,cp1,cq1,cu,cbj0,cby0,cbj1,cby1;
        thrust::complex<P> cyy,cbjk,ct2;
        thrust::complex<double> cii(0.0,1.0);
        thrust::complex<double> cone(1.0,0.0);
        thrust::complex<double> czero(0.0,0.0);
        P a0,y0;
        int k,m;
        static P a[] = {
                -0.7031250000000000e-1,
                0.1121520996093750,
                -0.5725014209747314,
                6.074042001273483};
        static P b[] = {
                0.7324218750000000e-1,
                -0.2271080017089844,
                1.727727502584457,
                -2.438052969955606e1};
        static P a1[] = {
                0.1171875,
                -0.1441955566406250,
                0.6765925884246826,
                -6.883914268109947};
        static P b1[] = {
                -0.1025390625,
                0.2775764465332031,
                -1.993531733751297,
                2.724882731126854e1};

        y0 = abs(imag(z));
        a0 = abs(z);
        nm = n;
        if (a0 < 1.0e-100) {
            for (k=0;k<=n;k++) {
                cj[k] = czero;
                cy[k] = thrust::complex<P> (-1e308,0);
                cjp[k] = czero;
                cyp[k] = thrust::complex<P>(1e308,0);
            }
            cj[0] = cone;
            cjp[1] = thrust::complex<P>(0.5,0.0);
            return 0;
        }
        if ((a0 <= 300.0) || (n > (int)(0.25*a0))) {
            if (n == 0) nm = 1;
            m = msta1(a0,200);
            if (m < nm) nm = m;
            else m = msta2(a0,nm,15);
            cbs = czero;
            csu = czero;
            csv = czero;
            cf2 = czero;
            cf1 = thrust::complex<P> (1.0e-100,0.0);
            for (k=m;k>=0;k--) {
                cf = 2.0*(k+1.0)*cf1/z-cf2;
                if (k <= nm) cj[k] = cf;
                if (((k & 1) == 0) && (k != 0)) {
                    if (y0 <= 1.0) {
                        cbs += 2.0*cf;
                    }
                    else {
                        cbs += (-1)*((k & 2)-1)*2.0*cf;
                    }
                    csu += (P)((-1)*((k & 2)-1))*cf/(P)k;
                }
                else if (k > 1) {
                    csv += (P)((-1)*((k & 2)-1)*k)*cf/(P)(k*k-1.0);
                }
                cf2 = cf1;
                cf1 = cf;
            }
            if (y0 <= 1.0) cs0 = cbs+cf;
            else cs0 = (cbs+cf)/cos(z);
            for (k=0;k<=nm;k++) {
                cj[k] /= cs0;
            }
            ce = log(0.5*z)+el;
            cy[0] = M_2_PI*(ce*cj[0]-4.0*csu/cs0);
            cy[1] = M_2_PI*(-cj[0]/z+(ce-1.0)*cj[1]-4.0*csv/cs0);
        }
        else {
            ct1 = z-M_PI_4;
            cp0 = cone;
            for (k=0;k<4;k++) {
                cp0 += a[k]*pow(z,-2.0*k-2.0);
            }
            cq0 = -0.125/z;
            for (k=0;k<4;k++) {
                cq0 += b[k] *pow(z,-2.0*k-3.0);
            }
            cu = sqrt(M_2_PI/z);
            cbj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1));
            cby0 = cu*(cp0*sin(ct1)+cq0*cos(ct1));
            cj[0] = cbj0;
            cy[0] = cby0;
            ct2 = z-0.75*M_PI;
            cp1 = cone;
            for (k=0;k<4;k++) {
                cp1 += a1[k]*pow(z,-2.0*k-2.0);
            }
            cq1 = 0.375/z;
            for (k=0;k<4;k++) {
                cq1 += b1[k]*pow(z,-2.0*k-3.0);
            }
            cbj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2));
            cby1 = cu*(cp1*sin(ct2)+cq1*cos(ct2));
            cj[1] = cbj1;
            cy[1] = cby1;
            for (k=2;k<=n;k++) {
                cbjk = 2.0*(k-1.0)*cbj1/z-cbj0;
                cj[k] = cbjk;
                cbj0 = cbj1;
                cbj1 = cbjk;
            }
        }
        cjp[0] = -cj[1];
        for (k=1;k<=nm;k++) {
            cjp[k] = cj[k-1]-(P)k*cj[k]/z;
        }
        if (abs(cj[0]) > 1.0)
            cy[1] = (cj[1]*cy[0]-2.0/(M_PI*z))/cj[0];
        for (k=2;k<=nm;k++) {
            if (abs(cj[k-1]) >= abs(cj[k-2]))
                cyy = (cj[k]*cy[k-1]-2.0/(M_PI*z))/cj[k-1];
            else
                cyy = (cj[k]*cy[k-2]-4.0*(k-1.0)/(M_PI*z*z))/cj[k-2];
            cy[k] = cyy;
        }
        cyp[0] = -cy[1];
        for (k=1;k<=nm;k++) {
            cyp[k] = cy[k-1]-(P)k*cy[k]/z;
        }

        return 0;
    }

    template<typename P>
    __host__ __device__
    int cbessjyva(P v, thrust::complex<P> z, P &vm, thrust::complex<P>*cjv,
                  thrust::complex<P>*cyv, thrust::complex<P>*cjvp, thrust::complex<P>*cyvp)
    {
        thrust::complex<P> z1,z2,zk,cjvl,cr,ca,cjv0,cjv1,cpz,crp;
        thrust::complex<P> cqz,crq,ca0,cck,csk,cyv0,cyv1,cju0,cju1,cb;
        thrust::complex<P> cs,cs0,cr0,cs1,cr1,cec,cf,cf0,cf1,cf2;
        thrust::complex<P> cfac0,cfac1,cg0,cg1,cyk,cp11,cp12,cp21,cp22;
        thrust::complex<P> ch0,ch1,ch2,cyl1,cyl2,cylk;
        thrust::complex<double> cii(0.0,1.0);
        thrust::complex<double> cone(1.0,0.0);
        thrust::complex<double> czero(0.0,0.0);

        P a0,v0,pv0,pv1,vl,ga,gb,vg,vv,w0,w1,ya0,yak,ya1,wa;
        int j,n,k,kz,l,lb,lb0,m;

        a0 = abs(z);
        z1 = z;
        z2 = z*z;
        n = (int)v;


        v0 = v-n;

        pv0 = M_PI*v0;
        pv1 = M_PI*(1.0+v0);
        if (a0 < 1.0e-100) {
            for (k=0;k<=n;k++) {
                cjv[k] = czero;
                cyv[k] = thrust::complex<P> (-1e308,0);
                cjvp[k] = czero;
                cyvp[k] = thrust::complex<P> (1e308,0);

            }
            if (v0 == 0.0) {
                cjv[0] = cone;
                cjvp[1] = thrust::complex<P> (0.5,0.0);
            }
            else {
                cjvp[0] = thrust::complex<P> (1e308,0);
            }
            vm = v;
            return 0;
        }
        if (z1.real() < 0.0) z1 = -z;
        if (a0 <= 12.0) {
            for (l=0;l<2;l++) {
                vl = v0+l;
                cjvl = cone;
                cr = cone;
                for (k=1;k<=40;k++) {
                    cr *= -0.25*z2/(k*(k+vl));
                    cjvl += cr;
                    if (abs(cr) < abs(cjvl)*eps) break;
                }
                vg = 1.0 + vl;
                ga = gamma(vg);
                ca = pow(0.5*z1,vl)/ga;
                if (l == 0) cjv0 = cjvl*ca;
                else cjv1 = cjvl*ca;
            }
        }
        else {
            if (a0 >= 50.0) kz = 8;
            else if (a0 >= 35.0) kz = 10;
            else kz = 11;
            for (j=0;j<2;j++) {
                vv = 4.0*(j+v0)*(j+v0);
                cpz = cone;
                crp = cone;
                for (k=1;k<=kz;k++) {
                    crp = -0.78125e-2*crp*(vv-pow(4.0*k-3.0,2.0))*
                          (vv-pow(4.0*k-1.0,2.0))/(k*(2.0*k-1.0)*z2);
                    cpz += crp;
                }
                cqz = cone;
                crq = cone;
                for (k=1;k<=kz;k++) {
                    crq = -0.78125e-2*crq*(vv-pow(4.0*k-1.0,2.0))*
                          (vv-pow(4.0*k+1.0,2.0))/(k*(2.0*k+1.0)*z2);
                    cqz += crq;
                }
                cqz *= 0.125*(vv-1.0)/z1;
                zk = z1-(0.5*(j+v0)+0.25)*M_PI;
                ca0 = sqrt(M_2_PI/z1);
                cck = cos(zk);
                csk = sin(zk);
                if (j == 0) {
                    cjv0 = ca0*(cpz*cck-cqz*csk);
                    cyv0 = ca0*(cpz*csk+cqz+cck);
                }
                else {
                    cjv1 = ca0*(cpz*cck-cqz*csk);
                    cyv1 = ca0*(cpz*csk+cqz*cck);
                }
            }
        }
        if (a0 <= 12.0) {
            if (v0 != 0.0) {
                for (l=0;l<2;l++) {
                    vl = v0+l;
                    cjvl = cone;
                    cr = cone;
                    for (k=1;k<=40;k++) {
                        cr *= -0.25*z2/(k*(k-vl));
                        cjvl += cr;
                        if (abs(cr) < abs(cjvl)*eps) break;
                    }
                    vg = 1.0-vl;
                    gb = gamma(vg);
                    cb = pow(2.0/z1,vl)/gb;
                    if (l == 0) cju0 = cjvl*cb;
                    else cju1 = cjvl*cb;
                }
                cyv0 = (cjv0*cos(pv0)-cju0)/sin(pv0);
                cyv1 = (cjv1*cos(pv1)-cju1)/sin(pv1);
            }
            else {
                cec = log(0.5*z1)+el;
                cs0 = czero;
                w0 = 0.0;
                cr0 = cone;
                for (k=1;k<=30;k++) {
                    w0 += 1.0/k;
                    cr0 *= -0.25*z2/(P)(k*k);
                    cs0 += cr0*w0;
                }
                cyv0 = M_2_PI*(cec*cjv0-cs0);
                cs1 = cone;
                w1 = 0.0;
                cr1 = cone;
                for (k=1;k<=30;k++) {
                    w1 += 1.0/k;
                    cr1 *= -0.25*z2/(k*(k+1.0));
                    cs1 += cr1*(2.0*w1+1.0/(k+1.0));
                }
                cyv1 = M_2_PI*(cec*cjv1-1.0/z1-0.25*z1*cs1);
            }
        }
        if (z.real() < 0.0) {
            cfac0 = exp(pv0*cii);
            cfac1 = exp(pv1*cii);
            if (z.imag() < 0.0) {
                cyv0 = cfac0*cyv0-(P)2.0*(thrust::complex<P>)cii*cos(pv0)*cjv0;
                cyv1 = cfac1*cyv1-(P)2.0*(thrust::complex<P>)cii*cos(pv1)*cjv1;
                cjv0 /= cfac0;
                cjv1 /= cfac1;
            }
            else if (z.imag() > 0.0) {
                cyv0 = cyv0/cfac0+(P)2.0*(thrust::complex<P>)cii*cos(pv0)*cjv0;
                cyv1 = cyv1/cfac1+(P)2.0*(thrust::complex<P>)cii*cos(pv1)*cjv1;
                cjv0 *= cfac0;
                cjv1 *= cfac1;
            }
        }
        cjv[0] = cjv0;
        cjv[1] = cjv1;
        if ((n >= 2) && (n <= (int)(0.25*a0))) {
            cf0 = cjv0;
            cf1 = cjv1;
            for (k=2;k<= n;k++) {
                cf = 2.0*(k+v0-1.0)*cf1/z-cf0;
                cjv[k] = cf;
                cf0 = cf1;
                cf1 = cf;
            }
        }
        else if (n >= 2) {
            m = msta1(a0,200);
            if (m < n) n = m;
            else  m = msta2(a0,n,15);
            cf2 = czero;
            cf1 = thrust::complex<P>(1.0e-45,0.0);
            for (k=m;k>=0;k--) {
                cf = 2.0*(v0+k+1.0)*cf1/z-cf2;
                if (k <= n) cjv[k] = cf;
                cf2 = cf1;
                cf1 = cf;
            }
            if (abs(cjv0) > abs(cjv1)) cs = cjv0/cf;
            else cs = cjv1/cf2;
            for (k=0;k<=n;k++) {
                cjv[k] *= cs;
            }
        }
        cjvp[0] = v0*cjv[0]/z-cjv[1];
        for (k=1;k<=n;k++) {
            cjvp[k] = -(k+v0)*cjv[k]/z+cjv[k-1];
        }
        cyv[0] = cyv0;
        cyv[1] = cyv1;
        ya0 = abs(cyv0);
        lb = 0;
        cg0 = cyv0;
        cg1 = cyv1;
        for (k=2;k<=n;k++) {
            cyk = 2.0*(v0+k-1.0)*cg1/z-cg0;
            yak = abs(cyk);
            ya1 = abs(cg0);
            if ((yak < ya0) && (yak< ya1)) lb = k;
            cyv[k] = cyk;
            cg0 = cg1;
            cg1 = cyk;
        }
        lb0 = 0;
        if ((lb > 4) && (z.imag() != 0.0)) {
            while(lb != lb0) {
                ch2 = cone;
                ch1 = czero;
                lb0 = lb;
                for (k=lb;k>=1;k--) {
                    ch0 = 2.0*(k+v0)*ch1/z-ch2;
                    ch2 = ch1;
                    ch1 = ch0;
                }
                cp12 = ch0;
                cp22 = ch2;
                ch2 = czero;
                ch1 = cone;
                for (k=lb;k>=1;k--) {
                    ch0 = 2.0*(k+v0)*ch1/z-ch2;
                    ch2 = ch1;
                    ch1 = ch0;
                }
                cp11 = ch0;
                cp21 = ch2;
                if (lb == n)
                    cjv[lb+1] = 2.0*(lb+v0)*cjv[lb]/z-cjv[lb-1];
                if (abs(cjv[0]) > abs(cjv[1])) {
                    cyv[lb+1] = (cjv[lb+1]*cyv0-2.0*cp11/(M_PI*z))/cjv[0];
                    cyv[lb] = (cjv[lb]*cyv0+2.0*cp12/(M_PI*z))/cjv[0];
                }
                else {
                    cyv[lb+1] = (cjv[lb+1]*cyv1-2.0*cp21/(M_PI*z))/cjv[1];
                    cyv[lb] = (cjv[lb]*cyv1+2.0*cp22/(M_PI*z))/cjv[1];
                }
                cyl2 = cyv[lb+1];
                cyl1 = cyv[lb];
                for (k=lb-1;k>=0;k--) {
                    cylk = 2.0*(k+v0+1.0)*cyl1/z-cyl2;
                    cyv[k] = cylk;
                    cyl2 = cyl1;
                    cyl1 = cylk;
                }
                cyl1 = cyv[lb];
                cyl2 = cyv[lb+1];
                for (k=lb+1;k<n;k++) {
                    cylk = 2.0*(k+v0)*cyl2/z-cyl1;
                    cyv[k+1] = cylk;
                    cyl1 = cyl2;
                    cyl2 = cylk;
                }
                for (k=2;k<=n;k++) {
                    wa = abs(cyv[k]);
                    if (wa < abs(cyv[k-1])) lb = k;
                }
            }
        }
        cyvp[0] = v0*cyv[0]/z-cyv[1];
        for (k=1;k<=n;k++) {
            cyvp[k] = cyv[k-1]-(k+v0)*cyv[k]/z;
        }
        vm = n+v0;
        return 0;
    }

///Calculate the spherical bessel functions and their derivatives up to order v
/// When allocating arrays to store the resulting values, arrays must be of size [v+2]
    template<typename P>
    __host__ __device__
    int cbessjyva_sph(int v, thrust::complex<P> z, P &vm, thrust::complex<P>*cjv,
                      thrust::complex<P>*cyv, thrust::complex<P>*cjvp, thrust::complex<P>*cyvp)
    {
        //first, compute the bessel functions of fractional order
        cbessjyva<P>(v + 0.5, z, vm, cjv, cyv, cjvp, cyvp);

        if(z == 0){													//handle degenerate case of z = 0
            memset(cjv, 0, sizeof(P) * (v+1));
            cjv[0] = 1;
        }

        //iterate through each and scale
        for(int n = 0; n<=v; n++)
        {
            if(z != 0){												//handle degenerate case of z = 0
                cjv[n] = cjv[n] * sqrt(M_PI/(z * 2.0));
                cyv[n] = cyv[n] * sqrt(M_PI/(z * 2.0));
            }

            cjvp[n] = -1.0 / (z * 2.0) * cjv[n] + cjvp[n] * sqrt(M_PI / (z * 2.0));
            cyvp[n] = -1.0 / (z * 2.0) * cyv[n] + cyvp[n] * sqrt(M_PI / (z * 2.0));
        }

        return 0;

    }

    template<typename P>
    __host__ __device__
    int chankelva_sph(const int v, const thrust::complex<P> z, thrust::complex<P>* chv) {
        P vm;
        int returnVal;
        thrust::complex<P>* cjv  = new thrust::complex<P>[v+2];
        thrust::complex<P>* cyv  = new thrust::complex<P>[v+2];
        thrust::complex<P>* cjvp = new thrust::complex<P>[v+2];
        thrust::complex<P>* cyvp = new thrust::complex<P>[v+2];

        returnVal = cbessjyva_sph<P>(v, z, vm, cjv, cyv, cjvp, cyvp);

        for (int i = 0; i <= v; i++) {
            chv[i] = cjv[i] + thrust::complex<P> (0, 1) * cyv[i];
        }

        delete[] cjv; delete[] cyv; delete[] cjvp; delete[] cyvp;
        return returnVal;
    }

    template<typename P>
    __host__ __device__
    int chankelvap_sph(const int v, const thrust::complex<P> z, thrust::complex<P>* chvp) {
        // to get derivative of upto order v we need to calculate hv upto order v+1
        int vTemp = v + 1;
        P vm;
        int returnVal;
        thrust::complex<P>* chv  = new thrust::complex<P>[vTemp + 1];

        returnVal = chankelva_sph<P>(vTemp, z, chv);
        chvp[0] = - chv[1];
        for (int i = 1; i <= v; i++) {
            chvp[i] = thrust::complex<P>(i, 0) * chv[i-1] - thrust::complex<P>(i+1, 0) * chv[i+1];
            chvp[i] = chvp[i] / thrust::complex<P> (2*i+1, 0);
        }

        return returnVal;
    }
}	//end namespace rts


#endif //PROJECT_BESSEL_CUH
