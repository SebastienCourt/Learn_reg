Sebastien Court, University of Innsbruck, Copyright 2023.

0/ This C++ code was written for performing the numerical experiments that are presented in the following preprint: 
https://arxiv.org/abs/2109.11658

Requirements: Linux system, libraries "qhull" and "MUMPS" installed.

1/ For compiling and running the code under Linux, first install the Getfem++ library:
        1.1/ Download: http://download-mirror.savannah.gnu.org/releases/getfem/stable/getfem-5.3.tar.gz
        1.2/ Uncompress the file: "tar xvfz getfem-5.3.tar.gz"
        1.3/ Create the folder "learn_reg" in "./getfem-5.3/contrib/" 
        1.4/ Copy the "Makefile.am" in "./getfem-5.3/contrib/learn_reg/"
        1.5/ Modify the file Makefile.am located in ./getfem-5.3/contrib/ by adding the sub-folder "learn_reg"
        1.6/ Modify the file "./getfem-5.3/configure.ac" by adding around line 1200 the string 
                "contrib/learn_reg/Makefile                                      \"

2/ Install Getfem as follows: In "./getfem-5.3/"
        2.1/ run successively "aclocal", "autoconf", "autoheader" and "automake --add-missing"
        2.2/ run "./configure --enable-qhull"
        2.3/ run "make"
        2.4/ run "sudo make install"

3/ Select the code "learn_reg" that you want to test:
        - either the one which tries to re-discover the L2-norm as regularizer (section 6.2 of the preprint).
        Then copy both files located from the folder "L2" to "/getfem-5.3/contrib/learn_reg/"
        - or the one which trains a NN that compensates the effect of a random noise (section 6.3 of the preprint).
        Then copy both files located from the folder "noise3" to "/getfem-5.3/contrib/learn_reg/"
        - or the one which trains a multi-dimensional NN which regularizes the problem of identification of several coefficients in an elliptic PDE (section 6.4 of the preprint), in order to compensate the effect of a random noise.
        Then copy both files located from folder "2D-codes" to "/getfem-5.3/contrib/learn_reg/"

4/  Compile the code:
        4.1/ First create a "MATLAB" folder in "./getfem-5.3/contrib/learn_reg/"
        4.2/ Run "make learn_reg"

5/ Execute the code by running "./learn_reg learn_reg.param". The results are contained in the "MATLAB" folder.
