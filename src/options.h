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

struct SolOptDB
{
	// option database
	PetscInt    skip_defaults                  =  0;                         // specify solver options explicitly
	PetscInt    view_solvers                   =  0;                         // show linear solver configuration
	PetscInt    monitor_solvers                =  0;                         // show linear iteration convergence
	PetscInt    set_linear_problem             =  0;                         // linear problem flag (skip nonlinear iteration)
	PetscScalar nonlinear_tolerances[3]        =  { 1e-5, -1.0, 50.0  };     // rtol, atol, maxit (-1 = automatic setting, only for atol)
	PetscScalar linear_tolerances[3]           =  { 1e-6, -1.0, 200.0 };     // rtol, atol, maxit (-1 = automatic setting, only for atol)
	PetscScalar picard_to_newton[4]            =  { 1e-2,  5.0, 1.2, 20.0 }; // picard rtol, picard minit, newton rtol, newton maxit
	PetscInt    use_line_search                =  1;
	PetscInt    use_eisenstat_walker           =  0;
	PetscInt    use_mat_free_jac               =  0;
	char        stokes_solver[_str_len_]       = "block_direct";             // [block_direct, coupled_mg, block_mg, wbfbt]
	char        direct_solver_type[_str_len_]  = "superlu_dist";             // [mumps, superlu_dist, lu]
	PetscScalar block_tolerances[2]            = { 1e-2, 30 } ;              // rtol, maxit (fgmres settings for block solves in block_mg and wbfbt)
	PetscScalar penalty                        =  1e3;                       // (only for block_direct)
	PetscInt    num_mg_levels                  = -1;                         // (-1 = automatic setting)
	PetscInt    num_mat_free_levels            =  0;                         // (only for coupled_mg)
	char        smoother_type[_str_len_]       = "heavy";                    // [light (richardson + jacobi), intermediate (chebyshev + sor), heavy (gmres + bjacobi)]
	char        smoother_ksp[_str_len_]        = {'\0'};                     // [richardson, chebyshev, gmres]
	char        smoother_pc[_str_len_]         = {'\0'};                     // [jacobi, sor, bjacobi, asm]
	PetscScalar smoother_damping               =  0.5;                       // (only for richardson)
	PetscScalar smoother_omega                 =  1.0;                       // (only for sor)
	PetscInt    smoother_num_sweeps            =  20;                        // maxit
	PetscInt    coarse_reduction_factor        = -1;                         // (-1 = automatic setting)
	PetscInt    coarse_cells_per_cpu           =  2048;                      // (-1 = all cpus are used by coarse solve)
	char        coarse_solver[_str_len_]       = "direct";                   // [direct, hypre, bjacobi, asm] (hypre, bjacobi and asm use fgmres)
	PetscScalar coarse_tolerances[2]           = { 1e-2, 30 } ;              // rtol, maxit (fgmres settings for hypre, bjacobi and asm)
	PetscInt    subdomain_overlap              =  1;                         // (only for asm)
	PetscInt    subdomain_ilu_levels           =  0;                         // (only for bjacobi and asm)
	PetscInt    subdomain_num_cells            = -1;                         // (only for bjacobi and asm) (-1 = one subdomain per cpu)
	char        init_thermal_solver[_str_len_] = "mg";                       // [mg, default]
	PetscScalar thermal_tolerances[3]          =  { 1e-8, -1.0, 500.0 };     // rtol, atol, maxit (-1 = automatic setting, only for atol)
	// computational parameters
	PetscInt    coarse_num_local_blocks                      =  0;
	PetscInt    levels_num_local_blocks[_max_num_mg_levels_] = {0};
	PetscInt    levels_num_blocks_constant                   =  0;
};

//-----------------------------------------------------------------------------

PetscErrorCode solverOptionsSetDefaults(FB *fb);

PetscErrorCode solverOptionsReadFromFile(FB *fb, SolOptDB &opt);

PetscErrorCode solverOptionsCheck(SolOptDB &opt);

PetscErrorCode get_num_mg_levels(SolOptDB &opt, FDSTAG *fs);

PetscErrorCode get_coarse_reduction_factor(
		SolOptDB &opt,
		PetscInt  coarse_num_local_cells);

PetscErrorCode get_num_local_blocks(
		SolOptDB &opt,
		PetscInt  levels_num_local_cells[],
		PetscInt  coarse_num_local_cells);

PetscErrorCode set_default_smoother(SolOptDB &opt);

PetscErrorCode set_smoother_options(
		SolOptDB   &opt,
		const char *prefix,
		PetscInt    num_local_blocks);

PetscErrorCode set_subdomain_options(
		SolOptDB   &opt,
		const char *prefix,
		const char *pc_type,
		PetscInt    num_local_blocks);

PetscErrorCode set_coarse_options(
		SolOptDB   &opt,
		const char *mg_prefix);

PetscErrorCode set_levels_options(
		SolOptDB   &opt,
		const char *mg_prefix);

PetscErrorCode set_custom_mg_options(
		SolOptDB   &opt,
		const char *prefix);

PetscErrorCode set_standard_mg_options(SolOptDB &opt, const char *prefix);

PetscErrorCode set_ksp_solver(const char *prefix, const char *type, PetscScalar rtol, PetscScalar maxit);

PetscErrorCode set_tolerances(const char *prefix, PetscScalar tolerances[3]);

PetscErrorCode set_integer_option(const char *key, const PetscInt val, const char *prefix = NULL);

PetscErrorCode set_scalar_option(const char *key, const PetscScalar val, const char *prefix = NULL);

PetscErrorCode set_string_option(const char *key, const char *val, const char *prefix = NULL);

PetscErrorCode set_empty_option(const char *key, const char *prefix = NULL);

//-----------------------------------------------------------------------------
// Driver routines
//-----------------------------------------------------------------------------

PetscErrorCode setSolverOptions(FB *fb);

PetscErrorCode PetscOptionsReadFromFile(FB *fb);

//-----------------------------------------------------------------------------

#endif
