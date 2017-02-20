# exact_diagonalization
An "exact diagonalization" program - mapping many-particle quantum mechanical Hamiltonians onto large-scale sparse eigenvalue problems

INSTALL:
autoconf
automake
chmod a+x configure

./configure --enable-blas --enable-lapack --enable-speed-optimization

DEPENDENCIES: boost program options

Google sparse-hash (if using speed optimization)
