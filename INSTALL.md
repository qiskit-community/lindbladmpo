# Installation Instructions

## Cygwin Installation

This stage is only for Windows users.

Go to Cygwin's official site:

https://www.cygwin.com

Download the executable and run the Installation process.

During the Installation, install the following Cygwin packages:

    gcc ,BLAS, LAPACK.

After finishing the Installation, open Cygwin terminal.

Note that in Cygwin the root directory of a particular drive (say "C"): 

    "C:/" 
is

    "/cygdrive/c/"

## [Optional] Open a projects directory

Run the following commands in your terminal, to create a directory under which the git repos will be downloaded:

    mkdir projects
    cd projects

## ITensor installation

Run the following commands in your terminal:

    mkdir itensor3
    cd itensor3
    git clone https://github.com/ITensor/ITensor.git .
    cp options.mk.sample options.mk

If you run Linux or Windows, you must un/comment as follows:

Linux/Windows only: Open options.mk in editor and comment out two MACOS lines, uncomment Linux BLAS, save.

Run the following command in your Linux terminal, or on Windows, navigate inside the Cygwin64 Terminal to the same directory and type:

    make

## Building the lindbladmpo solver

In order to leave the ITensor library and download the lindbladmpo library and build the solver, run the following commands in your terminal:

    cd ..
    mkdir lindbladmpo
    cd lindbladmpo
    git clone git@github.com:haggaila/lindbladmpo.git .
    cd src
    make

The above script will work for the default locations. In order to override where ITensor has been built, use:

    make LIBRARY_DIR="..."  - here insert the path for the itensor library, for example "/cygdrive/c/Users/Goran/itensor3"

In order to run a first test of the solver, type in the Linux terminal:

    cd ../bin/
    ./lindbladmpo

On Windows, navigate inside the Cygwin64 Terminal on Windows to /bin directory as above, and type:

    ./lindbladmpo
