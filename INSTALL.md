# Installation Instructions

## Cygwin Installation - Windows only

This step is for Windows users only.

Go to Cygwin's [official site](https://www.cygwin.com). Download the executable and run the installation process, making sure to install the following Cygwin packages:

    gcc, BLAS, LAPACK.

After the installation is done open the Cygwin64 Terminal App and follow the instructions below.

Note that in Cygwin the root directory of a particular drive (say "C"): 

    "C:\" 
is

    "/cygdrive/c/"

## [Optional] Open a projects directory

To create a particular directory under which the git repos will be downloaded, e.g., the folder 'projects', use the following commands in the terminal:

    mkdir projects
    cd projects

## ITensor installation

Execute the following commands in the terminal. This assumes that you have git locally configured:

    mkdir itensor3
    cd itensor3
    git clone https://github.com/ITensor/ITensor.git .

The above clone command uses https, while you can also clone using ssh if that is how you usually work with github (see the [ITensor](https://github.com/ITensor/ITensor) repository).

Now that ITensor is cloned locally, its binaries should be built. The default settings are configured for _Linux_ and _Windows_ (using Cygwin), with the basic BLAS/LAPACK installation.

On _Mac OS_ execute:

    cp options.mk.sample options.mk
    make

On _Linux_ or on _Windows_ (in the Cygwin64 Terminal App), there are two options:
- Use the first command as above to copy `options.mk.sample` into `options.mk` but then first edit `options.mk` in a text editor, commenting out the lines configuring the compilation and BLAS/LAPACK platform unsuitable for you platform and uncommenting the ones for your platform (the file is well documented). Save and exit the editor, typing `make` in the command line to build the binaries.
- Alternatively, within the lindbladmpo repository (see below) there is an `options.mk` file that can also be copied to the ITensor directory, where the choice of operating system is done once by defining the variable `OS_TARGET` at the top of the file, or setting it from the command-line invokation. See the details below and repeat the build of ITensor binaries.

## Building the lindbladmpo solver

In order to leave the ITensor folder and clone the lindbladmpo repository execute the following commands:

    cd ..
    mkdir lindbladmpo
    cd lindbladmpo
    git clone https://github.com/qiskit-community/lindbladmpo.git .

In order to build the solver binaries from the sources directory type:

    cd src
    make

The above script will work for the default locations, and on _Linux_ only (see below for _Mac OS_ and _Windows_). In order to override where ITensor has been built, use:

    make ITENSOR3_DIR="<path>"

where above, replace <path> with the path to the ITensor library, for example "/cygdrive/c/Users/Goran/itensor3".

In order to build on _Mac OS_ execute:

    make OS_TARGET=MACOS
    
while in order to build on _Windows_ execute:

    make OS_TARGET=WINDOWS

The default `OS_TARGET` parameter can also be set directly at the top of the `options.mk` file.

Finally, in order to run a simple test of the solver (with output to the console, and result files saved), execute:

    cd ../bin/
    ./lindbladmpo

Great! Now we suggest opening your favorite Python IDE and following the tutorial and examples - go back to the [Table of Contents](README.md#table-of-contents).

## Note on multithreading

The ITensor documentation and make options discuss a few possibilities for multithreading under Linux. The simplest option (which performs very well) is to install the OpenBLAS library, see https://www.openblas.net/ for further information.
