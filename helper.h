//
// Created by ragib1481 on 11/28/22.
//

#ifndef PROJECT_HELPER_H
#define PROJECT_HELPER_H

#include <thrust/host_vector.h>

thrust::host_vector<double> arange(const double start, const double end, const double dx) {
    /* returns array containing values [start, end) with dx increment
     */
    thrust::host_vector<double> values;
    double presentValue = start;
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

#endif //PROJECT_HELPER_H
