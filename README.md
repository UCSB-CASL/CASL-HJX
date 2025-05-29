# CASL-HJX

CASL-HJX is a C++ library for solving Hamilton-Jacobi equations and related partial differential equations using high-order numerical methods.

[![License: Academic](https://img.shields.io/badge/License-Academic-green.svg)](LICENSE)
[![C++](https://img.shields.io/badge/C%2B%2B-17%2B-blue.svg)](https://isocpp.org/)

## Features

- Hamilton-Jacobi-Bellman equations for optimal control problems
- Linear quadratic regulator (LQR) with SIMD optimization  
- Level-set methods for interface tracking
- High-order WENO and ENO spatial discretization
- TVD Runge-Kutta time integration
- IMEX methods for stiff problems
- MPI and OpenMP parallelization

The library supports various PDE types including advection, diffusion, Burgers equation, and advection-diffusion systems.

## Example: LQR Optimal Control

<p align="center">
  <img src="CASLProjects/projectLQR2D/Results/lqr_hero_professional.gif" alt="LQR optimal control demonstration" width="700"/>
</p>

## Installation

```bash
git clone https://github.com/UCSB-CASL/CASL-HJX.git
cd CASL-HJX
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

**Requirements:** C++17 compiler, CMake 3.15+, OpenMP 4.0+

## Quick Start

```bash
# Run LQR example
cd CASLProjects/projectLQR2D
./projectLQR2D 80 160

# Run other examples
cd ../projectAdvection && ./projectAdvection
cd ../projectBurgers && ./projectBurgers
```

## Documentation

- User Manual Paper: [https://arxiv.org/abs/2505.11527](https://arxiv.org/abs/2505.11527)
- Examples: `examples/` directory

## Authors

Faranak Rajabi, Jacob Fingerman, Andrew Wang, Jeff Moehlis, Frederic Gibou

**Computational Applied Systems Laboratory (CASL)**  
Department of Mechanical Engineering  
University of California, Santa Barbara  
[https://sites.engineering.ucsb.edu/~fgibou/index.html](https://sites.engineering.ucsb.edu/~fgibou/index.html)

## Citation

```bibtex
@article{rajabi2025caslhjx,
  title={CASL-HJX: A Comprehensive Guide to Solving Deterministic and 
         Stochastic Hamilton-Jacobi Equations},
  author={Rajabi, Faranak and Fingerman, Jacob and Wang, Andrew and 
          Moehlis, Jeff and Gibou, Frederic},
  journal={Computer Physics Communications},
  year={2025},
  note={Submitted}
}
```

## License

Academic research license. See LICENSE file for details.