bin_PROGRAMS = \
fqhe_sphere_spinless_fermions_diagonalization \
optical_flux_lattice_diagonalization \
optical_flux_lattice_single_particle_model \
optical_flux_lattice_sql_gen

##  Optical flux lattice diagonalization

optical_flux_lattice_diagonalization_SOURCES = \
optical_flux_lattice/main/optical_flux_lattice_diagonalization.cpp \
optical_flux_lattice/optical_flux_lattice_hamiltonian/optical_flux_lattice_hamiltonian.cpp \
optical_flux_lattice/optical_flux_lattice_hamiltonian/linear_momentum_constraints_2D.cpp \
optical_flux_lattice/observables_manager/observables_manager.cpp \
optical_flux_lattice/optical_flux_lattice_hamiltonian/optical_flux_lattice_hamiltonian_data.cpp \
optical_flux_lattice/single_particle_hamiltonian/single_particle_hamiltonian.cpp \
optical_flux_lattice/single_particle_hamiltonian/single_particle_hamiltonian_array.cpp \
optical_flux_lattice/single_particle_hamiltonian/single_particle_parameters.cpp \
optical_flux_lattice/optical_flux_lattice_hamiltonian/lookup_tables.cpp \
hamiltonians/hamiltonian_data.cpp \
basis/fermion_fock_basis.cpp  \
utilities/wrappers/mpi_wrapper.cpp \
utilities/wrappers/arpack_wrapper.cpp \
utilities/wrappers/fftw_wrapper.cpp \
utilities/wrappers/sqlite_wrapper.cpp \
utilities/mathematics/binary_number_tools.cpp \
utilities/mathematics/binomials.cpp \
utilities/mathematics/murmur_hash2.cpp \
utilities/mathematics/murmur_hash3.cpp

optical_flux_lattice_diagonalization_LINK = $(FC) $(LDFLAGS) -o optical_flux_lattice_diagonalization

##  Optical flux lattice single particle model

optical_flux_lattice_single_particle_model_SOURCES = \
optical_flux_lattice/main/optical_flux_lattice_single_particle_model.cpp \
optical_flux_lattice/single_particle_hamiltonian/single_particle_hamiltonian.cpp \
optical_flux_lattice/single_particle_hamiltonian/single_particle_hamiltonian_array.cpp \
optical_flux_lattice/single_particle_hamiltonian/single_particle_parameters.cpp \
utilities/wrappers/mpi_wrapper.cpp \
utilities/wrappers/fftw_wrapper.cpp \
utilities/wrappers/sqlite_wrapper.cpp \
utilities/mathematics/murmur_hash2.cpp \
utilities/mathematics/murmur_hash3.cpp

optical_flux_lattice_single_particle_model_LINK = $(FC) $(LDFLAGS) -o optical_flux_lattice_single_particle_model

##  Optical flux lattice SQL database generator

optical_flux_lattice_sql_gen_SOURCES = \
optical_flux_lattice/main/optical_flux_lattice_sql_gen.cpp \
utilities/wrappers/sqlite_wrapper.cpp

##  FQHE sphere spinless fermions diagonalization

fqhe_sphere_spinless_fermions_diagonalization_SOURCES = \
fqhe_sphere/main/fqhe_sphere_spinless_fermions_diagonalization.cpp \
fqhe_sphere/pseudopotential_hamiltonian/pseudopotential_hamiltonian.cpp \
fqhe_sphere/pseudopotential_hamiltonian/two_level_pseudopotential_hamiltonian.cpp \
fqhe_sphere/pseudopotential_hamiltonian/pseudopotential_hamiltonian_data.cpp \
fqhe_sphere/pseudopotential_hamiltonian/lookup_tables.cpp \
hamiltonians/hamiltonian_data.cpp \
basis/fermion_fock_basis.cpp \
basis/two_level_fermion_fock_basis.cpp \
utilities/wrappers/mpi_wrapper.cpp \
utilities/wrappers/arpack_wrapper.cpp \
utilities/mathematics/binary_number_tools.cpp \
utilities/mathematics/binomials.cpp \
utilities/mathematics/clebsch_gordan_coefficients.cpp \
utilities/mathematics/murmur_hash2.cpp \
utilities/mathematics/murmur_hash3.cpp


fqhe_sphere_spinless_fermions_diagonalization_LINK = $(FC) $(LDFLAGS) -o fqhe_sphere_spinless_fermions_diagonalization


