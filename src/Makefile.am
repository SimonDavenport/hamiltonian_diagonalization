bin_PROGRAMS = \
fqhe_sphere_spinless_fermions_diagonalization \
optical_flux_lattice_diagonalization \
optical_flux_lattice_noninteracting_model \
optical_flux_lattice_sql_gen

##  Optical flux lattice diagonalization

optical_flux_lattice_diagonalization_SOURCES = \
optical_flux_lattice/main/optical_flux_lattice_diagonalization.cpp \
optical_flux_lattice/interacting_ofl_model/interacting_ofl_model.cpp \
optical_flux_lattice/interacting_ofl_model/linear_momentum_constraints_2D.cpp \
optical_flux_lattice/observables_manager/observables_manager.cpp \
optical_flux_lattice/interacting_ofl_model/interacting_ofl_model_data.cpp \
optical_flux_lattice/noninteracting_ofl_model/noninteracting_ofl_model.cpp \
optical_flux_lattice/noninteracting_ofl_model/noninteracting_ofl_model_grid.cpp \
optical_flux_lattice/noninteracting_ofl_model/noninteracting_ofl_model_data.cpp \
optical_flux_lattice/interacting_ofl_model/term_tables.cpp \
optical_flux_lattice/interacting_ofl_model/term_hash_tables.cpp \
optical_flux_lattice/program_options/interacting_ofl_model_options.cpp \
hamiltonians/hamiltonian_data.cpp \
basis/fermion_fock_basis.cpp  \
program_options/general_options.cpp  \
utilities/wrappers/mpi_wrapper.cpp \
utilities/wrappers/arpack_wrapper.cpp \
utilities/wrappers/fftw_wrapper.cpp \
utilities/wrappers/sqlite_wrapper.cpp \
utilities/mathematics/binary_number_tools.cpp \
utilities/mathematics/binomials.cpp \
utilities/mathematics/murmur_hash2.cpp

optical_flux_lattice_diagonalization_LINK = $(FC) $(LDFLAGS) -o optical_flux_lattice_diagonalization

##  Optical flux lattice single particle model

optical_flux_lattice_noninteracting_model_SOURCES = \
optical_flux_lattice/main/optical_flux_lattice_noninteracting_model.cpp \
optical_flux_lattice/noninteracting_ofl_model/noninteracting_ofl_model.cpp \
optical_flux_lattice/noninteracting_ofl_model/noninteracting_ofl_model_grid.cpp \
optical_flux_lattice/noninteracting_ofl_model/noninteracting_ofl_model_data.cpp \
program_options/general_options.cpp  \
utilities/wrappers/mpi_wrapper.cpp \
utilities/wrappers/fftw_wrapper.cpp \
utilities/wrappers/sqlite_wrapper.cpp \
utilities/mathematics/murmur_hash2.cpp

optical_flux_lattice_noninteracting_model_LINK = $(FC) $(LDFLAGS) -o optical_flux_lattice_noninteracting_model

##  Optical flux lattice SQL database generator

optical_flux_lattice_sql_gen_SOURCES = \
optical_flux_lattice/main/optical_flux_lattice_sql_gen.cpp \
utilities/wrappers/sqlite_wrapper.cpp

##  FQHE sphere spinless fermions diagonalization

fqhe_sphere_spinless_fermions_diagonalization_SOURCES = \
fqhe_sphere/main/fqhe_sphere_spinless_fermions_diagonalization.cpp \
fqhe_sphere/pseudopotential_model/pseudopotential_model.cpp \
fqhe_sphere/pseudopotential_model/two_level_pseudopotential_model.cpp \
fqhe_sphere/pseudopotential_model/pseudopotential_model_data.cpp \
fqhe_sphere/pseudopotential_model/term_tables.cpp \
fqhe_sphere/interaction_pseudopotentials/coulomb_pseudopotentials.cpp \
hamiltonians/hamiltonian_data.cpp \
basis/fermion_fock_basis.cpp \
basis/two_level_fermion_fock_basis.cpp \
program_options/general_options.cpp  \
utilities/wrappers/mpi_wrapper.cpp \
utilities/wrappers/arpack_wrapper.cpp \
utilities/mathematics/binary_number_tools.cpp \
utilities/mathematics/binomials.cpp \
utilities/mathematics/clebsch_gordan_coefficients.cpp \
utilities/mathematics/murmur_hash2.cpp


fqhe_sphere_spinless_fermions_diagonalization_LINK = $(FC) $(LDFLAGS) -o fqhe_sphere_spinless_fermions_diagonalization


