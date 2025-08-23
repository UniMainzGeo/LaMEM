/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//......................   Default solver options   .........................
//---------------------------------------------------------------------------
#ifndef __options_h__
#define __options_h__

//-----------------------------------------------------------------------------

struct FB;
struct FDSTAG;

//-----------------------------------------------------------------------------
// Default solver options
//-----------------------------------------------------------------------------

PetscErrorCode solverOptionsSetDefaults(FB *fb);

PetscErrorCode solverOptionsSetRequired();

PetscErrorCode set_tolerances(const char *prefix, PetscScalar tolerances[3]);

PetscErrorCode get_num_mg_levels(
		FDSTAG  *fs,
		PetscInt &num_mg_levels);

PetscErrorCode get_coarse_reduction_factor(
		FDSTAG  *fs,
		PetscInt num_mg_levels,
		PetscInt coarse_cells_per_cpu,
		PetscInt &reduction_factor);

PetscErrorCode set_default_smoother(
		const char *smoother_type,
		char       *smoother_ksp,
		char       *smoother_pc);

PetscErrorCode set_smoother_options(
		const char *prefix,
		const char *smoother_ksp,
		const char *smoother_pc,
		PetscScalar smoother_damping,
		PetscScalar smoother_omega,
		PetscInt    smoother_num_sweeps,
		PetscInt    subdomain_overlap,
		PetscInt    subdomain_num_cells,
		PetscInt    num_local_cells);

PetscErrorCode set_subdomain_options(
		const char *prefix,
		const char *smoother_pc,
		PetscInt    subdomain_overlap,
		PetscInt    subdomain_num_cells,
		PetscInt    num_local_cells);

PetscErrorCode set_integer_option(const char *key, const PetscInt val, const char *prefix = NULL);

PetscErrorCode set_scalar_option(const char *key, const PetscScalar val, const char *prefix = NULL);

PetscErrorCode set_string_option(const char *key, const char *val, const char *prefix = NULL);

PetscErrorCode set_empty_option(const char *key, const char *prefix = NULL);

//-----------------------------------------------------------------------------
#endif
