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
#include "LaMEM.h"
#include "options.h"
#include "parsing.h"
#include "scaling.h"
#include "fdstag.h"
//-----------------------------------------------------------------------------
PetscErrorCode solverOptionsSetDefaults(FB *fb)
{
	// set "best-guess" solver options to help an inexperienced user
	// all options can be overridden by the usual PETSC options

	Scaling  scal_obj, *scal(&scal_obj);
	FDSTAG   fs_obj,   *fs  (&fs_obj);
	PetscInt complete_build, skip_defaults;

	PetscFunctionBeginUser;

	// read solver options block
	PetscCall(FBFindBlocks(fb, _OPTIONAL_, "<SolverOptionsStart>", "<SolverOptionsEnd>"));

	if(fb->nblocks)
	{
		if(fb->nblocks > 1)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many solver options blocks. Only one is allowed");
		}

		// check if no defaults should be applied
		skip_defaults =  0;

		PetscCall(getIntParam(fb, _OPTIONAL_, "skip_defaults", &skip_defaults, 1, 1));

		if(skip_defaults)
		{
			PetscCall(FBFreeBlocks(fb));

			PetscFunctionReturn(0);
		}
	}

	// set defaults
	PetscInt    set_linear_problem             =  0;                         // linear problem flag (skip nonlinear iteration)
	PetscScalar nonlinear_tolerances[ ]        =  { 1e-5, -1.0, 50.0  };     // rtol, atol, maxit (-1 = automatic setting, only for atol)
	PetscScalar linear_tolerances [ ]          =  { 1e-6, -1.0, 200.0 };     // rtol, atol, maxit (-1 = automatic setting, only for atol)
	PetscScalar picard_to_newton[ ]            =  { 1e-2,  5.0, 1.2, 20.0 }; // picard rtol, picard minit, newton rtol, newton maxit
	PetscInt    use_line_search                =  1;
	PetscInt    use_eisenstat_walker           =  0;
	PetscInt    use_mat_free_jac               =  0;
	char        stokes_solver[_str_len_]       = "block_direct";             // [block_direct, coupled_mg, block_mg, wbfbt]
	char        direct_solver_type[_str_len_]  = "superlu_dist";             // [mumps, superlu_dist, lu]
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
	char        coarse_solver[_str_len_]       = "direct";                   // [direct, hypre, bjacobi, asm] (hypre, bjacobi and asm use gmres)
	PetscScalar coarse_tolerances[]            = { 1e-3, 100 } ;             // rtol, maxit (gmres settings for hypre, bjacobi and asm)
	PetscInt    subdomain_overlap              =  1;                         // (only for asm)
	PetscInt    subdomain_ilu_levels           =  0;                         // (only for bjacobi and asm)
	PetscInt    subdomain_num_cells            = -1;                         // (-1 = automatic setting) (one subdomain per cpu)
	char        init_thermal_solver[_str_len_] = "mg";                       // [mg, default]
	PetscScalar thermal_tolerances[ ]          =  { 1e-8, -1.0, 500.0 };     // rtol, atol, maxit (-1 = automatic setting, only for atol)

	// read solver options block
	if(fb->nblocks)
	{
		// read simplified solver options
		PetscCall(getIntParam   (fb, _OPTIONAL_, "set_linear_problem",      &set_linear_problem,      1, 1));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "nonlinear_tolerances",     nonlinear_tolerances,    3, 1.0));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "linear_tolerances",        linear_tolerances,       3, 1.0));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "picard_to_newton",         picard_to_newton,        4, 1.0));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "use_line_search",         &use_line_search,         1, 1));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "use_eisenstat_walker",    &use_eisenstat_walker,    1, 1));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "use_mat_free_jac",        &use_mat_free_jac,        1, 1));
		PetscCall(getStringParam(fb, _OPTIONAL_, "stokes_solver",            stokes_solver,           NULL));
		PetscCall(getStringParam(fb, _OPTIONAL_, "direct_solver_type",       direct_solver_type,      NULL));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "penalty",                 &penalty,                 1, 1.0));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "num_mg_levels",           &num_mg_levels,           1, _max_num_mg_levels_));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "num_mat_free_levels",     &num_mat_free_levels,     1, _max_num_mat_free_levels_));
		PetscCall(getStringParam(fb, _OPTIONAL_, "smoother_type",            smoother_type,           NULL));
		PetscCall(getStringParam(fb, _OPTIONAL_, "smoother_ksp",             smoother_ksp,            NULL));
		PetscCall(getStringParam(fb, _OPTIONAL_, "smoother_pc",              smoother_pc,             NULL));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "smoother_damping",        &smoother_damping,        1, 1.0));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "smoother_omega",          &smoother_omega,          1, 1.0));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "smoother_num_sweeps",     &smoother_num_sweeps,     1, 1000));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "coarse_reduction_factor", &coarse_reduction_factor, 1, 1024));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "coarse_cells_per_cpu",    &coarse_cells_per_cpu,    1, 32768));
		PetscCall(getStringParam(fb, _OPTIONAL_, "coarse_solver",            coarse_solver,           NULL));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "coarse_tolerances",        coarse_tolerances,       2, 1.0));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "subdomain_overlap",       &subdomain_overlap,       1, 10));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "subdomain_ilu_levels",    &subdomain_ilu_levels,    1, 8));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "subdomain_num_cells",     &subdomain_num_cells,     1, 32768));
		PetscCall(getStringParam(fb, _OPTIONAL_, "init_thermal_solver",      init_thermal_solver,     NULL));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "thermal_tolerances",       thermal_tolerances,      3, 1.0));
	}

	// clear options block
	PetscCall(FBFreeBlocks(fb));

	// check parameters
	if(!(!strcmp(stokes_solver, "block_direct")
	||   !strcmp(stokes_solver, "coupled_mg")
	||   !strcmp(stokes_solver, "block_mg")
	||   !strcmp(stokes_solver, "wbfbt")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect Stokes solver type (stokes_solver): %s", stokes_solver);
	}

	if(!(!strcmp(direct_solver_type, "superlu_dist")
	||   !strcmp(direct_solver_type, "mumps")
	||   !strcmp(direct_solver_type, "lu")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect direct solver type (direct_solver_type): %s", direct_solver_type);
	}

	if(!(!strcmp(smoother_type, "light")
	||   !strcmp(smoother_type, "intermediate")
	||   !strcmp(smoother_type, "heavy")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect smoother type (smoother_type): %s", smoother_type);
	}

	// set default smoother if not specified
	PetscCall(set_default_smoother(smoother_type, smoother_ksp, smoother_pc));

	if(!(!strcmp(smoother_ksp, "richardson")
	||   !strcmp(smoother_ksp, "chebyshev")
	||   !strcmp(smoother_ksp, "gmres")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect smoother solver type (smoother_ksp): %s", smoother_ksp);
	}

	if(!(!strcmp(smoother_pc, "jacobi")
    ||   !strcmp(smoother_pc, "sor")
	||   !strcmp(smoother_pc, "bjacobi")
	||   !strcmp(smoother_pc, "asm")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect smoother preconditioner type (smoother_pc): %s", smoother_pc);
	}

	if(!(!strcmp(coarse_solver, "direct")
	||   !strcmp(coarse_solver, "hypre")
	||   !strcmp(coarse_solver, "bjacobi")
	||   !strcmp(coarse_solver, "asm")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect coarse solver type (coarse_solver): %s", coarse_solver);
	}

	if(!(!strcmp(init_thermal_solver, "mg")
	||   !strcmp(init_thermal_solver, "default")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect initial thermal solver type (init_thermal_solver): %s", init_thermal_solver);
	}

	//=================
	// NONLINEAR SOLVER
	//=================

	PetscCall(PetscOptionsInsertString(NULL, "-snes_monitor"));

	if(set_linear_problem)
	{
		PetscCall(PetscOptionsInsertString(NULL, "-snes_type ksponly"));
	}
	else
	{
		PetscCall(PetscOptionsInsertString(NULL, "-snes_max_funcs 1000000000"));

		PetscCall(set_tolerances("snes_", nonlinear_tolerances));

		if(use_line_search)
		{
			PetscCall(PetscOptionsInsertString(NULL, "-snes_linesearch_type l2"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_linesearch_max_it 5"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_linesearch_maxstep 1.0"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_linesearch_minlambda 0.05"));

		}

		if(use_eisenstat_walker)
		{
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_version"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_rtol0"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_rtolmax"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_gamma"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_alpha"));
		}

		PetscCall(set_scalar_option ("snes_picard_rtol",            picard_to_newton[0]));
		PetscCall(set_integer_option("snes_picard_minit", (PetscInt)picard_to_newton[1]));
		PetscCall(set_scalar_option ("snes_newton_rtol",            picard_to_newton[2]));
		PetscCall(set_integer_option("snes_newton_maxit", (PetscInt)picard_to_newton[3]));
	}

	//==============
	// LINEAR SOLVER
	//==============

	PetscCall(PetscOptionsInsertString(NULL, "-js_ksp_type fgmres"));
	PetscCall(PetscOptionsInsertString(NULL, "-js_ksp_monitor"));
	PetscCall(PetscOptionsInsertString(NULL, "-js_ksp_converged_reason"));

	PetscCall(set_tolerances("js_ksp_", linear_tolerances));

	if(use_mat_free_jac)
	{
		PetscCall(PetscOptionsInsertString(NULL, "-js_mat_free"));
	}

	//==========
	// MULTIGRID
	//==========

	// set multigrid flag
	if((!strcmp(stokes_solver,       "coupled_mg")
	||  !strcmp(stokes_solver,       "block_mg")
	||  !strcmp(stokes_solver,       "wbfbt")
	||  !strcmp(init_thermal_solver, "mg")))
	{
		// clear objects
		PetscCall(PetscMemzero(scal, sizeof(Scaling)));
		PetscCall(PetscMemzero(fs,   sizeof(FDSTAG)));

		// create scaling object
		PetscCall(ScalingCreate(scal, fb));

		// set link
		fs->scal = scal;

		// create grid (incomplete)
		PetscCall(FDSTAGCreate(fs, fb, complete_build = 0));

		// select number of multigrid levels
		PetscCall(get_num_mg_levels(fs, num_mg_levels));

		// select coarse solve reduction factor
		PetscCall(get_coarse_reduction_factor(fs, num_mg_levels, coarse_cells_per_cpu, coarse_reduction_factor));
/*

		subdomain_num_cells


		smoother_pc      bjacobi, asm
		coarse_solver    bjacobi, asm


*/

		// destroy grid object
		PetscCall(FDSTAGDestroy(fs));



	}




	//===============
	// PRECONDITIONER
	//===============

	if(!strcmp(stokes_solver, "block_direct"))
	{
		PetscCall(set_string_option("jp_type",                      "bf"));
		PetscCall(set_scalar_option("jp_pgamma",                    penalty));
		PetscCall(set_string_option("bf_vs_type",                   "user"));
		PetscCall(set_string_option("vs_ksp_type",                  "preonly"));
		PetscCall(set_string_option("vs_pc_type",                   "lu"));
		PetscCall(set_string_option("vs_pc_factor_mat_solver_type", direct_solver_type));
	}
	else if(!strcmp(stokes_solver, "coupled_mg"))
	{
		PetscCall(set_string_option ("jp_type", "mg"));

		// muligrid defaults
		// PetscCall(set_mg_options("gmg", nlevels, nsweeps, damping));

		if(num_mat_free_levels)
		{
			PetscCall(set_integer_option("gmg_mat_free_levels", num_mat_free_levels));

		}
	}
	else if(!strcmp(stokes_solver, "block_mg"))
	{
		PetscCall(set_string_option("jp_type",       "bf"));
		PetscCall(set_string_option("bf_vs_type",    "mg"));
		PetscCall(set_string_option("bf_schur_type", "inv_eta"));
		PetscCall(set_string_option("vs_ksp_type",   "preonly"));

		// muligrid defaults
		// PetscCall(set_mg_options("gmg", nlevels, nsweeps, damping));
	}
	else if(!strcmp(stokes_solver, "wbfbt"))
	{
		PetscCall(set_string_option("jp_type",       "bf"));
		PetscCall(set_string_option("bf_vs_type",    "mg"));
		PetscCall(set_string_option("bf_schur_type", "wbfbt"));
		PetscCall(set_string_option("vs_ksp_type",   "preonly"));
		PetscCall(set_string_option("ks_ksp_type",   "preonly"));

		// muligrid defaults
//		PetscCall(set_mg_options("gmg", nlevels, nsweeps, damping));
//		PetscCall(set_mg_options("ks",  nlevels, nsweeps, damping));
	}

/*
	-sub_pc_type ilu
	-sub_pc_factor_levels 1
	-sub_pc_factor_mat_ordering_type nd
	-sub_pc_factor_reuse_ordering
*/



	//===============
	// THERMAL SOLVER
	//===============

	PetscCall(set_tolerances("ts_ksp", thermal_tolerances));





	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode solverOptionsSetRequired()
{
	PetscFunctionBeginUser;

	// MAT
	PetscCall(PetscOptionsInsertString(NULL, "-mat_product_algorithm scalable"));
	PetscCall(PetscOptionsInsertString(NULL, "-matmatmatmult_via scalable"));
	PetscCall(PetscOptionsInsertString(NULL, "-matmatmult_via scalable"));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_tolerances(const char *prefix, PetscScalar tolerances[3])
{
	PetscFunctionBeginUser;

	PetscCall(set_scalar_option ("rtol", tolerances[0], prefix));

	if(tolerances[1] != -1.0)
	{
		PetscCall(set_scalar_option ("atol", tolerances[1], prefix));
	}
	else
	{
		PetscCall(set_empty_option ("atol_auto", prefix));
	}

	PetscCall(set_integer_option("max_it", (PetscInt)tolerances[2], prefix));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_default_smoother(
		const char *smoother_type,
		char       *smoother_ksp,
		char       *smoother_pc)
{
	PetscFunctionBeginUser;

	char smoother_ksp_default[_str_len_] = {'\0'};
	char smoother_pc_default [_str_len_] = {'\0'};

	if(!strcmp(smoother_type, "light"))
	{
		strcpy(smoother_ksp_default, "richardson");
		strcpy(smoother_pc_default,  "jacobi");
	}
	else if(!strcmp(smoother_type, "intermediate"))
	{
		strcpy(smoother_ksp_default, "chebyshev");
		strcpy(smoother_pc_default,  "sor");
	}
	else if(!strcmp(smoother_type, "heavy"))
	{
		strcpy(smoother_ksp_default, "gmres");
		strcpy(smoother_pc_default,  "bjacobi");
	}

	if(!strlen(smoother_ksp)) {  PetscCall(PetscStrncpy(smoother_ksp, smoother_ksp_default, _str_len_)); }
	if(!strlen(smoother_pc))  {  PetscCall(PetscStrncpy(smoother_pc,  smoother_pc_default,  _str_len_)); }

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_smoother_options(
		const char *prefix,
		const char *smoother_ksp,
		const char *smoother_pc,
		PetscScalar smoother_damping,
		PetscScalar smoother_omega,
		PetscInt    smoother_num_sweeps,
		PetscInt    subdomain_overlap,
		PetscInt    subdomain_num_cells,
		PetscInt    num_local_cells)
{
	PetscFunctionBeginUser;

	//====
	// KSP
	//====

	PetscCall(set_string_option ("ksp_type",   smoother_ksp,        prefix));
	PetscCall(set_integer_option("ksp_max_it", smoother_num_sweeps, prefix));

	if(!strcmp(smoother_ksp, "richardson"))
	{
		PetscCall(set_scalar_option("ksp_richardson_scale", smoother_damping, prefix));
	}
	else if(!strcmp(smoother_ksp, "gmres"))
	{
		PetscCall(set_integer_option("ksp_gmres_restart", smoother_num_sweeps, prefix));
	}

	//===
	// PC
	//===

	PetscCall(set_string_option ("pc_type", smoother_pc, prefix));

	if(!strcmp(smoother_pc, "sor"))
	{
		PetscCall(set_scalar_option("pc_sor_omega", smoother_omega, prefix));
	}
	else if(!strcmp(smoother_pc, "bjacobi")
	||      !strcmp(smoother_pc, "asm"))
	{
		PetscCall(set_subdomain_options(prefix, smoother_pc, subdomain_overlap, subdomain_num_cells, num_local_cells));
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_subdomain_options(
		const char *prefix,
		const char *smoother_pc,
		PetscInt    subdomain_overlap,
		PetscInt    subdomain_num_cells,
		PetscInt    num_local_cells)
{

	PetscInt num_local_blocks;

	PetscFunctionBeginUser;

	// compute number of local blocks
	if(subdomain_num_cells != 1)
	{
		num_local_blocks = PetscCeilInt(num_local_cells, subdomain_num_cells);
	}
	else
	{
		num_local_blocks = 1;
	}

	if(!strcmp(smoother_pc, "bjacobi"))
	{
		PetscCall(set_integer_option("pc_bjacobi_local_blocks", num_local_blocks, prefix));
	}
	else if(!strcmp(smoother_pc, "asm"))
	{
		PetscCall(set_integer_option("pc_asm_local_blocks", num_local_blocks,  prefix));
		PetscCall(set_integer_option("pc_asm_overlap",      subdomain_overlap, prefix));
		PetscCall(set_string_option ("pc_asm_type",         "restrict",        prefix));
		PetscCall(set_string_option ("pc_asm_local_type",   "additive",        prefix));
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode get_num_mg_levels(
		FDSTAG  *fs,
		PetscInt &num_mg_levels)
{
	// select number of multigrid levels

	PetscInt ncors;

	PetscFunctionBeginUser;

	// get maximum possible number of coarsening steps
	PetscCall(FDSTAGCheckMG(fs, ncors));

	if(num_mg_levels != -1)
	{
		// check user-specified number of levels
		if(num_mg_levels < 2 || num_mg_levels > ncors + 1)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect # of multigrid levels specified. Requested: %lld. Max. possible: %lld", (LLD)num_mg_levels, (LLD)(ncors + 1));
		}
	}
	else
	{
		num_mg_levels = ncors + 1;
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode get_coarse_reduction_factor(
		FDSTAG  *fs,
		PetscInt num_mg_levels,
		PetscInt coarse_cells_per_cpu,
		PetscInt &coarse_reduction_factor)
{
	// get number of processors for coarse grid solve

	PetscMPIInt size;
	PetscInt    nx, ny, nz, Nx, Ny, Nz, ncells;
	PetscInt    total_num_cpu, coarse_num_cpu, targer_factor, lower_factor, upper_factor;

	PetscFunctionBeginUser;

	// get number of ranks
	MPI_Comm_size(PETSC_COMM_WORLD, &size);

	total_num_cpu = (PetscInt)size;

	if(coarse_reduction_factor != -1)
	{
		// check user-specified reduction factor
		if(total_num_cpu % coarse_reduction_factor)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect reduction factor specified (coarse_reduction_factor): %lld", (LLD)coarse_reduction_factor);
		}
	}
	else if(total_num_cpu == 1)
	{
		// sequential case
		coarse_reduction_factor = 1;
	}
	else if(coarse_cells_per_cpu == -1)
	{
		// all processors are used for coarse solve
		coarse_reduction_factor = 1;
	}
	else
	{
		// get coarse grid size
		PetscCall(FDSTAGGetCoarseGridSize(fs, num_mg_levels, nx, ny, nz, Nx, Ny, Nz));

		// get total number of cells in the coarse grid
		ncells = Nx*Ny*Nz;

		// compute target number of processors
		coarse_num_cpu = PetscCeilInt(ncells, coarse_cells_per_cpu);

		// correct target number of processors
		if(coarse_num_cpu > total_num_cpu) { coarse_num_cpu = total_num_cpu; }

		// compute target reduction factor
		targer_factor = PetscCeilInt(total_num_cpu, coarse_num_cpu);

		// get lower estimate
		lower_factor = targer_factor; while(total_num_cpu % lower_factor) { lower_factor--; }

		// get upper estimate
		upper_factor = targer_factor; while(total_num_cpu % upper_factor) { upper_factor++; }

		// select optimal
		if((targer_factor - lower_factor) <= (upper_factor - targer_factor))
		{
			coarse_reduction_factor = lower_factor;
		}
		else
		{
			coarse_reduction_factor = upper_factor;
		}
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------







	// compute local grid size on all levels except the coarse


//	PetscInt levels_num_local_cells[num_mg_levels];

//	PetscCall(FDSTAGGetLevelsLocalGridSize(fs, num_mg_levels, levels_num_local_cells));




/*

PetscErrorCode set_coarse_options(
		const char *prefix,
		char        coarse_solver[],
		char        direct_solver_type[],
		PetscScalar coarse_tolerances[],
		PetscInt    coarse_reduction_factor,
		PetscInt    subdomain_overlap,
		PetscInt    subdomain_num_cells,
		PetscInt    num_local_cells)
{

	PetscFunctionBeginUser;


	PetscCall(set_scalar_option("ksp_gmres_restart", smoother_num_sweeps, prefix));
	PetscCall(set_integer_option("ksp_gmres_restart", smoother_num_sweeps, prefix));
	PetscCall(set_string_option ("pc_type", smoother_pc, prefix));


	//	ierr = FDSTAGGetCoarseGridSize(fs, nlevels, nx, ny, nz, Nx, Ny, Nz); CHKERRQ(ierr);


//===
// PC
//===

	coarse_num_cpu        = -1                    # (-1 = automatic setting) (-2 = all cpus are used by coarse solve)

	coarse_cells_per_cpu  =  2048                 # (only required for automatic setting of coarse_num_cpu parameter)

	coarse_solver         =  direct               # [direct, hypre, bjacobi, asm]

	coarse_tolerances     =  1e-3 100             # rtol, maxit (only for bjacobi and asm, since they use gmres)
PetscCall(set_string_option ("pc_type", smoother_pc, prefix));

if(!strcmp(smoother_pc, "sor"))
{
	PetscCall(set_scalar_option("pc_sor_omega", smoother_omega, prefix));
}
else if(!strcmp(smoother_pc, "bjacobi")


		PetscCall(set_integer_option("pc_asm_overlap",      subdomain_overlap, prefix));
		PetscCall(set_string_option ("pc_asm_type",         "restrict",        prefix));


	-gmg_mg_coarse_ksp_type preonly
	-gmg_mg_coarse_pc_type telescope
	-gmg_mg_coarse_pc_telescope_reduction_factor 2

	[A] direct solver
	-gmg_mg_coarse_telescope_ksp_type preonly
	-gmg_mg_coarse_telescope_pc_type lu
	-gmg_mg_coarse_telescope_pc_factor_mat_solver_type superlu_dist



	[B] hypre boomeramg
	-gmg_mg_coarse_telescope_ksp_type preonly


	pc_type hypre
	pc_hypre_type boomeramg
	pc_hypre_boomeramg_agg_nl             5
	pc_hypre_boomeramg_max_iter           1
	pc_hypre_boomeramg_relax_weight_all   0.8
	pc_hypre_boomeramg_strong_threshold   0.9
	pc_hypre_boomeramg_smooth_type        Euclid
	pc_hypre_boomeramg_eu_bj
	pc_hypre_boomeramg_coarsen_type       HMIS
	pc_mg_galerkin_mat_product_algorithm  hypre
	pc_hypre_boomeramg_coarsen_type       HMIS
	pc_mg_galerkin_mat_product_algorithm  hypre




	PetscCall(set_string_option ("pc_type",                             "hypre",     prefix));
	PetscCall(set_string_option ("pc_hypre_type",                       "boomeramg", prefix));
	PetscCall(set_integer_option("pc_hypre_boomeramg_agg_nl",           5,           prefix));
	PetscCall(set_integer_option("pc_hypre_boomeramg_max_iter",         1,           prefix));
	PetscCall(set_scalar_option ("pc_hypre_boomeramg_relax_weight_all", 0.8,         prefix));
	PetscCall(set_integer_option("pc_hypre_boomeramg_grid_sweeps_all",  10,          prefix));
	PetscCall(set_scalar_option ("pc_hypre_boomeramg_strong_threshold", 0.9,         prefix));
	PetscCall(set_string_option ("pc_hypre_boomeramg_smooth_type",      "Euclid",    prefix));
	PetscCall(set_empty_option  ("pc_hypre_boomeramg_eu_bj",                         prefix));
	PetscCall(set_string_option ("pc_hypre_boomeramg_coarsen_type",     "HMIS",      prefix));


	PetscFunctionReturn(0);

}

//-----------------------------------------------------------------------------
PetscErrorCode set_mg_options(const char *prefix, PetscInt nlevels, PetscInt nsweeps, PetscScalar damping)
{
	PetscFunctionBeginUser;

	PetscCall(set_string_option ("pc_type", "mg", prefix));
	PetscCall(set_integer_option("pc_mg_levels", nlevels, prefix));
	PetscCall(set_empty_option  ("pc_mg_galerkin", prefix));
	PetscCall(set_string_option ("pc_mg_type", "multiplicative", prefix));
	PetscCall(set_string_option ("pc_mg_cycle_type", "v", prefix));








	PetscFunctionReturn(0);
}



*/




PetscErrorCode set_integer_option(const char *key, const PetscInt val, const char *prefix)
{
	PetscFunctionBeginUser;
	char *opt;
	if(prefix) asprintf(&opt,"-%s_%s %lld", prefix, key, (LLD)val);
	else       asprintf(&opt,"-%s %lld",            key, (LLD)val);
	PetscCall(PetscOptionsInsertString(NULL, opt));
	free(opt);
	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_scalar_option(const char *key, const PetscScalar val, const char *prefix)
{
	PetscFunctionBeginUser;
	char *opt;
	if(prefix) asprintf(&opt,"-%s_%s %g", prefix, key, val);
	else       asprintf(&opt,"-%s %g",            key, val);
	PetscCall(PetscOptionsInsertString(NULL, opt));
	free(opt);
	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_string_option(const char *key, const char *val, const char *prefix)
{
	PetscFunctionBeginUser;
	char *opt;
	if(prefix) asprintf(&opt,"-%s_%s %s", prefix, key, val);
	else       asprintf(&opt,"-%s %s",            key, val);
	PetscCall(PetscOptionsInsertString(NULL, opt));
	free(opt);
	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_empty_option(const char *key, const char *prefix)
{
	PetscFunctionBeginUser;
	char *opt;
	if(prefix) asprintf(&opt,"-%s_%s", prefix, key);
	else       asprintf(&opt,"-%s",            key);
	PetscCall(PetscOptionsInsertString(NULL, opt));
	free(opt);
	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
