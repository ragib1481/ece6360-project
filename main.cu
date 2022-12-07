#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <complex.h>
#include <chrono>

#include "helper.h"
#include "myMath.h"

using namespace std;
using namespace std::chrono;


int main(int argc, char* argv[]) {
    // parse command line arguments
    if (argc != 2)
        return 1;
    string filename(argv[1]);

    // define simulation parameters.
    Parameters parameters("./" + filename);

    cout << "Bessel Function" << endl;

    for (int n = 0; n < 5; n++) {
        for (float x = 0.01; x <= 10.0; x++) {
            cout << n << " " << x << " " << redefined::spBesselN<float>(n, x) - std::sph_bessel<float>(n, x) << endl;
        }
        cout << endl;
    }
    cout << endl;
    cout << endl;

    cout << "Legendre Polynomial" << endl;
    for (unsigned int n = 0; n <= 1; n++) {
        for (float x = 0.0; x <= 10.0; x++) {
            cout << redefined::legendre<float>(n, x) - std::legendre<float>(n, x) << ", ";
        }
        cout << endl;
    }
    cout << endl;
    cout << endl;

    cout << "hankel function" << endl;

    return 0;
}
