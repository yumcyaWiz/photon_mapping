# photon_mapping

minimal but extensible implementation of photon mapping in C++.

## Requirements

* C++ (20>=)
* CMake (3.20>=)
* [spdlog](https://github.com/gabime/spdlog)
* OpenMP

## Build

CMake options.

|Option|Description|
|:--|:--|
|BUILD_TESTS|build tests|

```
mkdir build
cd build
cmake ..
```

## Structure

|Name|Description|
|:--|:--|
|`include/camera.h`|ray generation from camera|
|`include/core.h`|math, basic data types|
|`include/image.h`|image. PPM output is supported.|
|`include/integrator.h`|implement photon mapping, path tracing(for reference)|
|`include/photon_map.`|implement photon map with kdtree|
|`include/primitive.h`|primitive object|
|`include/sampler.h`|random number generation, sampling utilities|
|`include/scene.h`|scene object, implement ray-scene intersection|
|`include/shape.h`|implement ray-primitive intersection, sampling a point on primitive|

## References

* Jensen, Henrik Wann. Realistic image synthesis using photon mapping. AK Peters/crc Press, 2001.
* https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/Stochastic_Progressive_Photon_Mapping# 
* http://www.cs.cmu.edu/afs/cs/academic/class/15462-s12/www/lec_slides/lec18.pdf