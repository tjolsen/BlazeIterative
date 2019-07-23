## Blaze-Iterative
===============

This is a set of iterative linear system solvers intended for use
with the [Blaze library](https://bitbucket.org/blaze-lib/blaze/src/master/), a high-performance C++ linear algebra library.
The API is currently based on a tag-dispatch system to choose a particular algorithm.

### Currently implemented algorithms:
 #### [Conjugate Gradient](https://github.com/STEllAR-GROUP/BlazeIterative/blob/master/docs/Conjugate%20Gradient.md) 
 #### BiCGSTAB
 #### Preconditioned BiCGSTAB
 #### [Arnoldi](https://github.com/STEllAR-GROUP/BlazeIterative/blob/master/docs/Arnoldi.md)

### Planned algorithms:
#### [Preconditioned CG](https://github.com/STEllAR-GROUP/BlazeIterative/blob/master/docs/Precondition%20Conjugate%20Gradient.md)
#### (Preconditioned) BiCGSTAB(l)
#### [GMRES](https://github.com/STEllAR-GROUP/BlazeIterative/blob/master/docs/GMRES.md)
#### [Lanczos](https://github.com/STEllAR-GROUP/BlazeIterative/blob/master/docs/Lanczos.md)

### Potential algorithms (if sufficient interest):
- LSQR
- LSMR


Please open an issue to discuss/request features.


Installation
------------
Installing BlazeIterative is done via the usual CMake procedure:
- create an out-of-tree build directory (highly recommended)
- run cmake
- build and install using your favorite build systems

An example of the process on a linux-like system might be (from the BlazeIterative root):
```bash
mkdir build && cd build
ccmake .. # This is where you set your install location, defaults to /usr/local on my machine
make && make install
```

The library headers will install to `${CMAKE_INSTALL_PREFIX}/include/BlazeIterative`.
The cmake files that allow you to `find_package(BlazeIterative)` are installed in
`${CMAKE_INSTALL_PREFIX}/share/BlazeIterative/cmake`.
If `CMAKE_INSTALL_PREFIX` is on your `CMAKE_PREFIX_PATH`, then `find_package(BlazeIterative)` should
work without any issues (tested with cmake 3.9).


Using BlazeIterative in your projects
-------------------------------------
BlazeIterative is developed using CMake, so it is easy to use in another CMake project.
The installation process exports the target: "BlazeIterative::BlazeIterative", which
should be linked against. This will take care of linking against Blaze, but you still
must find the blaze package in your own project (if anyone knows how to solve the transitive
dependency problem in cmake, please teach me!).
This is a minimal example of a CMakeLists.txt for a project using BlazeIterative.

```cmake
project(MyReallyComplicatedProject)

find_package(blaze)
find_package(BlazeIterative)

add_executable(myCoolExecutable main.cpp)
target_link_libraries(myCoolExecutable PUBLIC BlazeIterative::BlazeIterative)
```

