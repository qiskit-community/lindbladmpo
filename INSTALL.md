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

    "C:\" 
is

    "/cygdrive/c/"

## [Optional] Open a projects directory

Run the following commands in the terminal, to create a directory under which the git repos will be downloaded:

    mkdir projects
    cd projects

## ITensor installation

Run the following commands in the terminal (on Windows, replace the "cp" command with "copy"). This assumes that you have a git locally configured:

    mkdir itensor3
    cd itensor3
    git clone https://github.com/ITensor/ITensor.git .
    cp options.mk.sample options.mk

Linux/Windows only: Open options.mk in editor and comment out two MACOS lines, uncomment Linux BLAS, and save.

Run the following command in the Linux terminal, or on Windows, navigate inside the Cygwin64 Terminal to the same directory and type:

    make

## Building the lindbladmpo solver

In order to leave the ITensor library and download the lindbladmpo repo and build the solver, run the following commands in the terminal:

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

Or, on Windows, navigate inside the Cygwin64 Terminal on Windows to /bin directory as above, and type:

    ./lindbladmpo
