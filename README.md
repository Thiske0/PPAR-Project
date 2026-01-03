# PPAR Project

Made by Pascal Kessler and Mathis Poppe

## Installation

### Windows

1. Install MSYS2 (<https://www.msys2.org/>)
2. Open the mingw64 shell, located under "C:/msys64/mingw64.exe"
3. Run the following commands

    ```bash
    pacman -Sy
    pacman -Syu
    pacman -Sy
    pacman -S mingw-w64-x86_64-toolchain
    pacman -S mingw-w64-x86_64-openmp
    pacman -S mingw-w64-x86_64-msmpi
    pacman -S mingw-w64-x86_64-gnuplot
    pacman -S mingw-w64-x86_64-mpich
    ```

4. Add "C:\msys64\mingw64\bin" to PATH
5. Reopen VS-code

### Mac

1. Install gcc and gnuplot
2. Install MPI and OpenMP if not already present
3. Reopen VS-code

### Linux

1. Install gcc and gnuplot
2. Install MPI and OpenMP if not already present
3. Reopen VS-code

## Running the project

The project is meant to be compiled and run on Paradoxe in Rennes.
Succesfull compilation is not guaranteed on other systems due to possible lack of libraries.

### manual compilation

1. run the command specified at the top of each file to compile the program
2. run the executable using the required arguments as listed in the second line

### using script

1. run either the `./run.sh` or the `./bench.sh` script
