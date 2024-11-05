# CASLHJB2D

A C++ library for solving Hamilton-Jacobi-Bellman equations in 2D, with specific applications for Linear Quadratic Regulator (LQR) and other control problems.

Here's the edited and merged version for your README.md:

```markdown
## Project Structure
```
CASLHJB2D/                  # Root directory
├── CMakeLists.txt         # Main CMake file
├── CASLCommonLibrary/     # Core library code
│   ├── CASLArray2D.cpp           # 2D array implementation
│   ├── CaslGrid2D.h              # Grid management
│   ├── CaslOptions.h             # Configuration options
│   ├── CASLCppToMatlab2D.cpp     # MATLAB export utilities
│   ├── CASLHamiltonian2D.cpp     # 2D Hamiltonian solver
│   ├── CASLHamiltonJacobi.cpp    # Hamilton-Jacobi solver
│   └── CASLHamiltonJacobi2D.cpp  # 2D HJB implementation
└── CASLProjects/          # Individual projects
    ├── projectLQR2D/             # Linear Quadratic Regulator
    │   ├── CMakeLists.txt        # Project-specific CMake
    │   ├── main.cpp              # Main implementation
    │   └── projectLQR2D_lib/     # LQR specific library
    │       ├── CaslLQRSystemDynamics.h
    │       └── CaslHamiltonianLQR2D.h
    ├── projectLaplacian2D_1D/    # 2D to 1D Laplacian solver
    ├── projectStochasticHH2D/    # Stochastic HJB (planned)
    └── projectDeterministicHH2D/ # Deterministic HJB (planned)

# Build artifacts (ignored by git)
├── cmake-build-debug/     # Debug build directory
└── cmake-build-release/   # Release build directory
```

## Prerequisites

- CMake (version 3.15 or higher)
- C++17 compatible compiler
- For macOS: LLVM/Clang

## Building the Project

### macOS

1. Install dependencies:
```bash
brew install llvm cmake
```

2. Clone and build:
```bash
git clone https://github.com/yourusername/CASLHJB2D.git
cd CASLHJB2D
mkdir build && cd build
cmake ..
make
```

### Running Tests

Each project in the `CASLProjects` directory can be built and run independently:

```bash
cd build/CASLProjects/projectLQR2D
./projectLQR2D
```

## Project Components

### CASLCommonLibrary
Core components used across all projects:
- `CASLArray2D`: 2D array implementation
- `CaslGrid2D`: 2D grid management
- `CASLHamiltonJacobi2D`: HJB equation solver

### Projects

1. **projectLQR2D**: Linear Quadratic Regulator in 2D
   - Solves optimal control problems
   - Implements system dynamics and Hamiltonian

2. **projectLaplacian2D_1D**: Laplacian solver
   - Handles 2D to 1D dimensional reduction
   - Implements boundary conditions

3. **projectStochasticHH2D**: (In Development)
   - Will implement stochastic Hamilton-Jacobi equations

4. **projectDeterministicHH2D**: (In Development)
   - Will implement deterministic Hamilton-Jacobi equations

## Output Structure
Each project creates its output in an `__Output` directory:
```
projectLQR2D/
└── __Output/
    └── LQR2D_[gridsize]/
        └── phi/
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Author

Faranak Rajabi
