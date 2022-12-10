//
// Created by ragib1481 on 11/28/22.
//

#ifndef PROJECT_HELPER_CUH
#define PROJECT_HELPER_CUH

#include <thrust/host_vector.h>
#include <fstream>
#include <string>


template <typename T>
struct doublet {
    T x = 0; T y = 0;
};

class Parameters {
public:
    doublet<unsigned short> dim;                        // dimension of sample along x and y
    doublet<double> sampling;                            // sampling interval along x and Y
    doublet<unsigned short> numSamples;                 // samples along x and y
    double na;                                           // numerical aperture of focusing optics
    double k;                                            // wavenubmer for simulation
    double zStart;
    double zEnd;
    double zSampling;

    // default constructor
    Parameters() {
        dim.x           = 0;
        dim.y           = 0;
        sampling.x      = 0;
        sampling.y      = 0;
        numSamples.x    = 0;
        numSamples.y    = 0;
        na              = 0;
        k               = 0;
        zStart          = 0;
        zEnd            = 0;
        zSampling       = 0;
    }

    // construct object from specified parameter file
    Parameters(std::string filename) {
        std::string discard;
        std::ifstream infile;
        infile.open(filename);
        if (infile.is_open()) {
            infile >> discard;
            infile >> dim.x;
            infile >> dim.y;

            infile >> discard;
            infile >> sampling.x;
            infile >> sampling.y;

            numSamples.x = (unsigned short) ((double)dim.x / sampling.x);
            numSamples.y = (unsigned short) ((double)dim.y / sampling.y);

            infile >> discard;
            infile >> na;

            infile >> discard;
            infile >> k;                                                            // load wavelength
            k = 2 * M_PI / k;                                                       // wavelength to wave-number

            infile >> discard;
            infile >> zStart;

            infile >> discard;
            infile >> zEnd;

            infile >> discard;
            infile >> zSampling;

            infile.close();
        }
    };

    // copy constructor
    Parameters(const Parameters& params) {
        dim              = params.dim;
        sampling         = params.sampling;
        numSamples       = params.numSamples;
        na               = params.na;
        k                = params.k;
        zStart           = params.zStart;
        zEnd             = params.zEnd;
        zSampling        = params.zSampling;
    }

    // overload assignment operator
    Parameters& operator=(const Parameters& params) {
        dim                 = params.dim;
        sampling            = params.sampling;
        numSamples          = params.numSamples;
        na                  = params.na;
        k                   = params.k;
        zStart              = params.zStart;
        zEnd                = params.zEnd;
        zSampling           = params.zSampling;

        return *this;
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

__host__ __device__
int factorial(int n) {
    int val = 1;
    for (int i = 2; i <= n; i++) {
        val *= i;
    }
    return val;
}

#endif //PROJECT_HELPER_CUH
