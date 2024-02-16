#!/bin/sh

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make

cd .. 

# example of generating 4096 2D points from a 64^2 toroidal tile (1spp) over 10 frames using 8000 Adam step  and 3000 projection per step
./build/SpatioTemporalFSOT -n 4096 -p 8000 --nbproj 3000 --tileSize 64 --frames 10 -d 2 -o ./results/SpatioTemporalFSOT_64_64_4spp_10_frames.dat
