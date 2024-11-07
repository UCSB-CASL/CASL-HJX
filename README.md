# CASLHJB2D

A C++ library for solving Hamilton-Jacobi-Bellman equations in 2D, with specific applications for Linear Quadratic Regulator (LQR) and other control problems.

## Features

- 2D Hamilton-Jacobi-Bellman equation solver
- Linear Quadratic Regulator (LQR) implementation
- Laplacian solver with dimensional reduction
- MATLAB export capabilities
- Flexible grid management system
- Configurable solver options

## Project Structure

```
CASLHJB2D/                  
├── CMakeLists.txt          # Main CMake configuration
├── CASLCommonLibrary/      # Core library components
│   ├── CASLArray2D.cpp     # 2D array implementation
│   ├── CaslGrid2D.h        # Grid management
│   ├── CaslOptions.h       # Configuration options
│   ├── CASLCppToMatlab2D.cpp     # MATLAB export utilities
│   ├── CASLHamiltonian2D.cpp     # 2D Hamiltonian solver
│   ├── CASLHamiltonJacobi.cpp    # Hamilton-Jacobi solver
│   └── CASLHamiltonJacobi2D.cpp  # 2D HJB implementation
└── CASLProjects/           # Project implementations
    ├── projectLQR2D/       # Linear Quadratic Regulator
    │   ├── CMakeLists.txt
    │   ├── main.cpp
    │   └── projectLQR2D_lib/
    │       ├── CaslLQRSystemDynamics.h
    │       └── CaslHamiltonianLQR2D.h
    ├── projectLaplacian2D_1D/     # 2D to 1D Laplacian solver
    ├── projectStochasticHH2D/     # Stochastic HJB (planned)
    └── projectDeterministicHH2D/  # Deterministic HJB (planned)
```

## Prerequisites

- CMake 3.15 or higher
- C++17 compatible compiler
- For macOS: LLVM/Clang

## Installation

### macOS

1. Install required dependencies:
```bash
brew install llvm cmake
```

2. Clone and build the project:
```bash
git clone https://github.com/yourusername/CASLHJB2D.git
cd CASLHJB2D
mkdir build && cd build
cmake ..
make
```

## Usage

### Running Projects

Each project in the `CASLProjects` directory can be built and run independently:

```bash
cd build/CASLProjects/projectLQR2D
./projectLQR2D
```

### Output Structure

Project outputs are organized as follows:

```
projectLQR2D/
└── __Output/
    └── LQR2D_[gridsize]/
        └── phi/
```

## Components

### Core Library (CASLCommonLibrary)

- **CASLArray2D**: Efficient 2D array implementation
- **CaslGrid2D**: Advanced grid management system
- **CASLHamiltonJacobi2D**: HJB equation solver
- **CASLCppToMatlab2D**: MATLAB integration utilities

### Projects

1. **projectLQR2D**
   - Implements Linear Quadratic Regulator in 2D
   - Includes system dynamics and Hamiltonian solvers
   - Optimized for control problems

2. **projectLaplacian2D_1D**
   - Handles 2D to 1D dimensional reduction
   - Implements various boundary conditions
   - Efficient Laplacian solving algorithms

3. **projectStochasticHH2D** (In Development)
   - Will implement stochastic Hamilton-Jacobi equations
   - Planned support for uncertainty handling

4. **projectDeterministicHH2D** (In Development)
   - Will implement deterministic Hamilton-Jacobi equations
   - Focus on precise control solutions

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

Faranak Rajabi

## Contact

For questions and support, please open an issue in the GitHub repository.
