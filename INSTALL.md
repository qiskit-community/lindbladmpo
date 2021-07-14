# Installation Instructions

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

Run the following commands in your terminal:

    cd ..
    mkdir lindbladmpo
    cd lindbladmpo/
    git clone git@github.com:haggaila/lindbladmpo.git .
    cd src
    make
    cd ../bin/
    ./lindbladmpo

## Cygwin Installation

Go to Cygwin's official site:

https://www.cygwin.com

Download the executable and preform a fresh install.