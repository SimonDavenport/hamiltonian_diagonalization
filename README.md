# exact_diagonalization
An "exact diagonalization" program - mapping many-particle quantum mechanical Hamiltonians onto large-scale sparse eigenvalue problems. This code was used for calculating results presented in this research paper: https://arxiv.org/abs/1504.07232.

INSTALL:
autoconf
automake
chmod +x configure

./configure --enable-blas --enable-lapack --enable-speed-optimization

DEPENDENCIES: boost program options

Google sparse-hash (if using speed optimization)
