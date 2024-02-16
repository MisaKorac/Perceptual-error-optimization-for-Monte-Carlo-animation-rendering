# Perceptual error optimization for Monte Carlo animation rendering

Source code of "Perceptual error optimization for Monte Carlo animation rendering". This example is for demonstration purposes and has not been optimized. It is fully CPU based with a parallelization via OpenMP.

In case of problems or bugs don't hesitate to contact me I will do my best to solve it.

Dependancies:
=============
 + OpenMP (`brew install libomp` on macOS)

Code compilation:
=================

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make
    cd ..

Example:
=================

    ./build/SpatioTemporalFSOT -n 4096 -p 8000 --nbproj 3000 --tileSize 64 --frames 10 -d 2 -o ./results/SpatioTemporalFSOT_64_64_4spp_10_frames.dat
