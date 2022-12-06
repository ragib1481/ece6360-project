#include <iostream>
#include <thrust/host_vector.h>
#include <complex.h>
#include <chrono>

#include "helper.h"

using namespace std;
using namespace std::chrono;


int main(int argc, char* argv[]) {
    // parse command line arguments
    if (argc != 2)
        return 1;
    string filename(argv[1]);

    // define simulation parameters.
    Parameters parameters("./" + filename);

    return 0;
}
