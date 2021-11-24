# Installation Instructions

## Cygwin Installation - Windows only

This stage is only for Windows users.

Go to Cygwin's official site:

https://www.cygwin.com

Download the executable and run the Installation process.

During the Installation, install the following Cygwin packages:

    gcc, BLAS, LAPACK.

After the installation is finished, open the Cygwin terminal, and follow the instructions below.

Note that in Cygwin the root directory of a particular drive (say "C"): 

    "C:\" 
is

    "/cygdrive/c/"

## [Optional] Open a projects directory

To create a particular directory under which the git repos will be downloaded - e.g. the direction 'projects', use the following commands in the terminal:

    mkdir projects
    cd projects

## ITensor installation

Run the following commands in the terminal (on Windows, replace the "cp" command with "copy"). This assumes that you have git locally configured:

    mkdir itensor3
    cd itensor3
    git clone https://github.com/ITensor/ITensor.git .
    cp options.mk.sample options.mk

Open options.mk in a text editor and if necessary, comment/uncomment the lines configuring the compilation and BLAS/LAPACK platform according to the operating system (Linux, Mac OS or Cygwin on Windows), and save. By default, no editing is required on Linux or Windows (using cygwin), with the basic BLAS/LAPACK installation.

Alternatively, within the lindbladmpo repository (see below) there is an options.mk file that can also be copied to the itensor directory, where the choice of operating system is done once by defining a variable at the top of the file (or setting it from the command line invokation of make).

Run the following command in the Linux/Mac OS terminal, or on Windows - navigate inside the Cygwin64 Terminal to the same directory and type:

    make

## Building the lindbladmpo solver

In order to leave the ITensor library and download the lindbladmpo repo and build the solver, run the following commands in the terminal:

    cd ..
    mkdir lindbladmpo
    cd lindbladmpo
    git clone https://github.com/haggaila/lindbladmpo.git .
    cd src
    make

The above script will work for the default locations. In order to override where ITensor has been built, use:

    make LIBRARY_DIR="..."

where above, replace three dots with the path for the itensor library, for example "/cygdrive/c/Users/Goran/itensor3".

In order to run a first test of the solver (with output to the console, and result files saved), type in the Linux terminal:

    cd ../bin/
    ./lindbladmpo

Or, on Windows, navigate inside the Cygwin64 Terminal on Windows to /bin directory as above, and type:

    ./lindbladmpo

## Note on multithreading

The ITensor documentation and make options discuss a few possibilities for multithreading under Linux. The simplest option (which performs very well) is to install the OpenBLAS library, see https://www.openblas.net/ for further information.
