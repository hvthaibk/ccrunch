/*
 * Copyright (C) 2011 Florian Rathgeber, florian.rathgeber@gmail.com
 *
 * This code is licensed under the MIT License.
 * See the FindCUDA.cmake script for the text of the license.
 *
 * Based on code by Christopher Bruns published on Stack Overflow (CC-BY):
 * http://stackoverflow.com/questions/2285185
 */

#include <stdio.h>
#include <cuda_runtime.h>

int main()
{
    int                     deviceCount, device, major = 9999, minor = 9999;
    int                     gpuDeviceCount = 0;
    struct cudaDeviceProp   properties;

    if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
        // ref: http://stackoverflow.com/questions/23067320/no-cuda-capable-device-is-detected-using-ubuntu-12-04-4-server
        printf("running cudaGetDeviceCount() failed!\n");
        return 1;
    }

    /* machines with no GPUs can still report one emulation device */
    for (device = 0; device < deviceCount; ++device) {
        cudaGetDeviceProperties(&properties, device);

        if (properties.major != 9999) {/* 9999 means emulation only */
            ++gpuDeviceCount;

            /*  get minimum compute capability of all devices */
            if (major > properties.major) {
                major = properties.major;
                minor = properties.minor;
            } else if (minor > properties.minor) {
                minor = properties.minor;
            }
        }
    }

    /* don't just return the number of gpus, because other runtime cuda
       errors can also yield non-zero return values */
    if (gpuDeviceCount > 0) {
        /* this output will be parsed by FindCUDA.cmake */
        printf("%d%d", major, minor);
        return 0; /* success */
    }
    return 1; /* failure */
}
