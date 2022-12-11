//
// Created by ragib1481 on 11/28/22.
//

#ifndef PROJECT_HELPER_CUH
#define PROJECT_HELPER_CUH

#include <thrust/host_vector.h>
#include <fstream>
#include <string>

#include "Sphere.cuh"
#include "Vec2.cuh"
#include "Vec3.cuh"


template <typename T>
class Parameters {
public:
    Vec2<T> sampling;                               // sampling interval along first and second dimension
    T xStart;                                       // start of x
    T xEnd;                                         // end of x
    T zStart;                                       // start of z
    T zEnd;                                         // end of z
    T na;                                           // numerical aperture of focusing optics
    T k;                                            // wavenumber for simulation
    Vec2<T> kUnitVector;                            // direction of k
    Vec2<T> focus;                                  // focal point
    thrust::host_vector<Sphere<T>> spheres;         // list of spheres in the mesh

    // default constructor
    Parameters() {
        xStart          = 0;
        xEnd            = 0;
        zStart          = 0;
        zEnd            = 0;
        na              = 0;
        k               = 0;
        sampling        = 0;
        kUnitVector     = 0;
    }

    // construct object from specified parameter file
    Parameters(std::string filename) {
        std::string discard;
        std::ifstream infile;
        infile.open(filename);
        if (infile.is_open()) {
            infile >> discard;
            infile >> xStart;
            infile >> discard;
            infile >> xEnd;

            infile >> discard;
            infile >> zStart;
            infile >> discard;
            infile >> zEnd;

            T samplingTemp;
            infile >> discard;
            infile >> samplingTemp;     sampling.setElement(0, samplingTemp);
            infile >> samplingTemp;     sampling.setElement(1, samplingTemp);

            infile >> discard;
            infile >> na;

            infile >> discard;
            infile >> k;                                                            // load wavelength
            k = 2 * M_PI / k;                                                       // wavelength to wave-number

            T kUnitVectorTemp;
            infile >> discard;
            infile >> kUnitVectorTemp;
            kUnitVector.setElement(0, kUnitVectorTemp);
            infile >> kUnitVectorTemp;
            kUnitVector.setElement(1, kUnitVectorTemp);

            T focusTemp;
            infile >> discard;
            infile >> focusTemp;
            focus.setElement(0, focusTemp);
            infile >> focusTemp;
            focus.setElement(1, focusTemp);

            T temp1, temp2, r;
            infile >> discard;
            while (infile.peek() !=EOF) {
                // load sphere center
                infile >> temp1;
                infile >> temp2;
                Vec2<T> c(temp1, temp2);

                // load radius;
                infile >> r;

                // load refractive index
                infile >> temp1;
                infile >> temp2;
                thrust::complex<T> n(temp1, temp2);

                Sphere<T> sp(c, r, n, k);
                spheres.push_back(sp);
            }

            infile.close();
        }
    }

    // copy constructor
    Parameters(const Parameters& params) {
        xStart           = params.xStart;
        xEnd             = params.xEnd;
        zStart           = params.zStart;
        zEnd             = params.zEnd;
        sampling         = params.sampling;
        na               = params.na;
        k                = params.k;
        kUnitVector      = params.kUnitVector;
        focus            = params.focus;
        spheres          = params.spheres;
    }

    // overload assignment operator
    Parameters& operator=(const Parameters& params) {
        xStart           = params.xStart;
        xEnd             = params.xEnd;
        zStart           = params.zStart;
        zEnd             = params.zEnd;
        sampling         = params.sampling;
        na               = params.na;
        k                = params.k;
        kUnitVector      = params.kUnitVector;
        focus            = params.focus;
        spheres          = params.spheres;

        return *this;
    }

    void print() {
        std::cout << "x:                    " << xStart << " " << xEnd << std::endl;
        std::cout << "z:                    " << zStart << " " << zEnd << std::endl;
        std::cout << "sampling:             " << sampling << std::endl;
        std::cout << "numerical aperture:   " << na << std::endl;
        std::cout << "wavelength:           " << 2*M_PI/k << std::endl;
        std::cout << "k dir:                " << kUnitVector << std::endl;
        std::cout << "focus:                " << focus << std::endl;

        for (int i =0; i < spheres.size(); i++) {
            std::cout << spheres[i].center() << " " << spheres[i].radius() << " "
                      << spheres[i].refractiveIndex() << std::endl;
        }
    }
};

template <typename T>
thrust::host_vector<T> arange(const T start, const T end, const T dx) {
    /* returns array containing values [start, end) with dx increment
     */
    thrust::host_vector<T> values;
    T presentValue = start;
    while (presentValue < end) {
        values.push_back(presentValue);
        presentValue += dx;
    }
    return values;
}

thrust::host_vector<double> linspace(const double start, const double end, const int num) {
    /* returns `num equidistant values from `start` to `end`
     */
    thrust::host_vector<double> values;
    double presentVal = start;
    double dx = (end - start) / (double)(num - 1);
    while (values.size() < num) {
        values.push_back(presentVal);
        presentVal += dx;
    }
    return values;
}

thrust::host_vector<double> fftfreq(const int n, const double d, const double mult = 1.0) {
    /* Returns discrete frequency components given the number of samples and sampling.
     * n -> number of samples
     * d -> sampling
     * mult -> scaling factors for each frequency value
     */

    thrust::host_vector<double> freq(n);
    if (n%2 == 0) {
        for (int i = 0; i < freq.size(); i++) {
            freq[i] = mult * (-n/2 + i) / (d * n);
        }
    }
    else {
        for (int i = 0; i < freq.size(); i++) {
            freq[i] = mult * (-(n-1)/2 + i) / (d * n);
        }
    }
    return freq;
}

template <typename T>
void fftshift(thrust::host_vector<T>& signal) {
    /* shift zero frequency component to the center of the array
     * or undo the shift if the zero frequency is already at the center.
     * The function performs the operation ***inplace***
     * WARNING: ***this function only works for even number of samples***
     */
    T temp;
    int i = 0;
    int j = signal.size() / 2;

    while (i < signal.size()/2) {
        temp = signal[i];
        signal[i] = signal[j];
        signal[j] = temp;
        i++; j++;
    }
}

template <typename T>
void fftshift2d(thrust::host_vector<T>& data, const int width, const int height) {
    /* shift zero frequency component to the center of the 2d array
     * or undo the shift if the zero frequency is already at the center.
     * The function performs the operation ***inplace***
     * WARNING: ***this function only works for even number of samples***
     */
    T temp;
    int j;
    int k;

    // shift along the width
    for (int i = 0; i < height; i++) {
        j = 0;
        k = width / 2;

        while ( j < width / 2) {
            temp = data[i * width + j];
            data[i * width + j] = data[i * width + k];
            data[i * width + k] = temp;
            j++; k++;
        }
    }

    // shift along the height
    for (int i = 0; i < width; i++) {
        j = 0;
        k = height / 2;

        while(j < height / 2) {
            temp = data[i + j * width];
            data[i + j * width] = data[i + k * width];
            data[i + k * width] = temp;
            j++; k++;
        }
    }
}

#endif //PROJECT_HELPER_CUH
