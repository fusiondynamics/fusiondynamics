#fusiondynamics

fusiondynamics is a physics student project that aims to study the fusion dynamics of heavy ions.

##Installation

    git clone https://github.com/fusiondynamics/fusiondynamics
    cd fusiondynamics/source/libs/liblangevin/testing
    make

##Usage

Assuming you are in the testing directory:

    ./turfu_langevin > data/mydatafile.data

<!--There are 3 columns in the data file: dimensionless time, dimensionless position, and fusion probability. You can exploit these data using the gnuplot script results.gplt.-->
There are 2 columns in the data file: dimensionless time, and fusion probability. You can exploit these data using the gnuplot script results.gplt.

##Stuff that is not working at the moment

<!--* For low kinetic energy values, probability data are too high to be physically explained. We are currently working on this issue.-->
* Probability data are messed up. We are currently working on this issue.

##Possible enhancements

We are looking to add a fourth column in data files that will represent the first derivative of probability. 
