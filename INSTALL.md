# Installation Instructions

## Cygwin Installation

This stage is only for Windows users.

Go to Cygwin's official site:

https://www.cygwin.com

Download the executable and Run the Installation process.

During the Installation, install the following Cygwin packages:

    gcc ,BLAS, LAPACK.

After finishing the Installation, open Cygwin terminal.

Note that in Cygwin the directory 

    "C:/Users/" 
is

    "/cygdrive/c/Users/"

## Open a project directory

Run the following commands in your terminal:

    mkdir gitprojects
    cd gitprojects`

## ITensor installation

Run the following commands in your terminal:

    mkdir itensor3
    cd itensor3
    git clone https://github.com/ITensor/ITensor.git 
    cp options.mk.sample options.mk

If you run Linux or Windows, you must un/comment as follows:

* Linux: Open options.mk in editor and comment out two MACOS lines, uncomment Linux BLAS, save.

* Windows: Open options.mk in editor and comment out two MACOS lines save.


Run the following command in your terminal:

    make

## Building the Library

Nevigate to the path you want to install the library and run the following commands in your terminal:

    mkdir lindbladmpo
    cd lindbladmpo/
    git clone git@github.com:haggaila/lindbladmpo.git .
    cd src
    make LIBRARY_DIR="..."  - here insert the path for the itensor library, for example "/cygdrive/c/Users/Goran/itensor3"
    cd ../bin/
    ./lindbladmpo
