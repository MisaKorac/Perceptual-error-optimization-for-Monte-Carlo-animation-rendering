# Perceptual error optimization for Monte Carlo animation rendering

Open-source code for ["Perceptual error optimization for Monte Carlo animation rendering"](https://misakorac.com/publications/2023_korac_perceptual_error_optimization_for_animation_rendering.html). The code is parallelized on CPU using OpenMP and tested on Ubuntu 22.04 LTS and Windows 11 23H2. It is built upon that of ["Filtered Sliced Optimal Transport
"](https://github.com/iribis/filtered-sliced-optimal-transport) and ["Sliced Optimal Transport Sampling"](https://github.com/loispaulin/Sliced-Optimal-Transport-Sampling), and we would like to take this opportunity to thank the authors.

![Teaser](/Assets/Teaser.jpg)

Note that the code is for demonstration purposes only. Still, in case of any issues or questions, do not hesitate to contact us: (Miša Korać: misa.korac@dfki.de, korac@cg.uni-saarland.de or Corentin Salaün: csalaun@mpi-inf.mpg.de).

Dependancies:
=============
 + OpenMP

Code compilation:
=================

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    cmake --build .

Toy example:
=================

    ./SpatioTemporalFSOT -n 4096 -p 8000 --nbproj 3000 --tileSize 64 --frames 10 -d 2 -o tile_64x64x10_1spp.dat

Code outputs .dat file and .h file which can be directly included inside a C/C++ code. For loading .dat files, please check demonstrative pre-computed tiles at [the project webpage](https://misakorac.com/publications/2023_korac_perceptual_error_optimization_for_animation_rendering.html).
