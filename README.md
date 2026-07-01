# Quantum Error Correction (QEC) with QuEST

This project implements quantum error correction simulations using the QuEST quantum state simulator.

The project prompts the user to enter text. The text is then converted into bits, which are individually encoded into the Steane [[7,1,3]] and Shor [[9,1,3]] quantum error-correcting codes. These encoded states are then passed through a noise channel with a per-qubit error rate specified by the user. The corresponding error correction protocols are then applied in an attempt to recover the original information.

The project also computes the von Neumann entropy of qubit subsystems for both the [[7,1,3]] and [[9,1,3]] codes.

The project is written in C and uses CMake for the build system.

-------------------------------------

## Dependencies

This project depends on:

- [QuEST](https://github.com/QuEST-Kit/QuEST) 
- LAPACK
- LAPACKE
- BLAS
- CMake ≥ 3.16
- C compiler (gcc/clang)

--------------------------------------

## Clone the repository


```bash
git clone https://github.com/Ciaran-DalCais/Quantum-Error-Correction.git
mkdir external && cd external && mkdir QuEST && cd QuEST 
git clone https://github.com/QuEST-Kit/QuEST.git    
--------------------------------------

## Build Instructions

From the project root: 

```bash

cmake -S . -B build
cmake --build build

---------------------------------------

## Run the simulation

```bash
./qec_sim


---------------------------------------

## The Project Structure

QEC/
├── external/
│   └── QuEST/        # git submodule
├── include/          # project headers
├── src/              # source code
├── CMakeLists.txt
└── README.md


-----------------------------------------

## License

This project is licensed under the MIT License. See the LICENSE file for details.


