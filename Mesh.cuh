//
// Created by ragib1481 on 12/7/22.
//

#ifndef PROJECT_MESH_CUH
#define PROJECT_MESH_CUH

#include <thrust/host_vector.h>
#include "Vec3.cuh"
#include "Vec2.cuh"
#include "helper.cuh"
#include "Sphere.cuh"


namespace mesh {

    namespace cpu {

        template <typename T>
        class CartesianMesh {
            /* This data structure defines a mesh in the cartesian co-ordinate system.
             * The constructor takes the limits of the mesh and sampling parameters and defines a mesh.
             * The x direction is taken along the width and the y/z is taken along the height.
             * Each element is stored as a Point in a one dimensional vector in the row major order.
             * */
            unsigned int width, height;

            thrust::host_vector<Vec2<T>> grid;

        public:
            CartesianMesh() {
                width  = static_cast<unsigned int>(0);
                height = static_cast<unsigned int>(0);
            }

            CartesianMesh(T xInit, T zInit, T xEnd, T zEnd, T xSampling, T zSampling) {

                width  = static_cast<unsigned int>((xEnd - xInit) / xSampling);
                height = static_cast<unsigned int>((zEnd - zInit) / zSampling);

                grid.resize(width * height);

                auto x = arange<T>(xInit, xEnd, xSampling);
                auto z = arange<T>(zInit, zEnd, zSampling);

                for (unsigned int j = 0; j < height; j++) {
                    for(unsigned int i = 0; i < width; i++) {
                        Vec2<T> p(x[i], z[j]);
                        grid[j * width + i] = p;
                    }
                }
            }

            unsigned int getWidth() { return width;}
            unsigned int getHeight() { return height;}

            Vec2<T> getPoint(unsigned int j, unsigned int i) {
                /* returns the coordinate point at row/y/z = j, column/x = i*/
                return grid[j * width + i];
            }
        };

    }
} // mesh

#endif //PROJECT_MESH_CUH
