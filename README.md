fusiondynamics
==============

fusiondynamics is a physics student project that aims to study the fusion dynamics of heavy ions.

Installation
------------

    git clone https://github.com/fusiondynamics/fusiondynamics
    cd fusiondynamics/source/libs/liblangevin/testing
    make

Usage
-----

Assuming you are in the testing directory:

    ./test_fusion > data/mydatafile.data

There are 3 columns in the data file: dimensionless time, dimensionless position, and fusion probability. You can exploit these data using the gnuplot script results.gplt.

TODO
----

* First derivative of fusion probability
* Two-dimensional approach
* Code/program optimization
