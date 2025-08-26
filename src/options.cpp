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
PetscErrorCode setSolverOptions(FB *fb)
{
	char *all_options;

	PetscFunctionBeginUser;

	// copy all command line and previously specified options to buffer
	PetscCall(PetscOptionsGetAll(NULL, &all_options));

	// remove command line options from database
	PetscCall(PetscOptionsClear(NULL));

	// set simplified solver options from the input file
	PetscCall(solverOptionsSetDefaults(fb));

	// load additional options from file
	PetscCall(PetscOptionsReadFromFile(fb));

	// push command line options to the end of database (priority)
	PetscCall(PetscOptionsInsertString(NULL, all_options));

	// clean
	PetscCall(PetscFree(all_options));

	// list entire option database
	PetscCall(PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode solverOptionsSetDefaults(FB *fb)
{
	// set "best-guess" solver options to help an inexperienced user
	// all options can be overridden by the usual PETSC options

	SolOptDB opt;
	Scaling   scal_obj, *scal(&scal_obj);
	FDSTAG    fs_obj,   *fs  (&fs_obj);
	PetscInt  act_temp_diff(0), act_steady_temp(0), complete_build, MG2D;
	PetscInt  levels_num_local_cells [_max_num_mg_levels_], coarse_num_local_cells;

	PetscFunctionBeginUser;

	// read file
	PetscCall(solverOptionsReadFromFile(fb, opt));

	// skip if explicitly requested
	if(opt.skip_defaults) {	PetscFunctionReturn(0); }

	// check options
	PetscCall(solverOptionsCheck(opt));

	//=================
	// NONLINEAR SOLVER
	//=================

	PetscCall(PetscOptionsInsertString(NULL, "-snes_monitor"));

	if(opt.set_linear_problem)
	{
		PetscCall(PetscOptionsInsertString(NULL, "-snes_type ksponly"));
	}
	else
	{
		PetscCall(PetscOptionsInsertString(NULL, "-snes_max_funcs 1000000000"));

		PetscCall(set_tolerances("snes", opt.nonlinear_tolerances));

		if(opt.use_line_search)
		{
			PetscCall(PetscOptionsInsertString(NULL, "-snes_linesearch_type l2"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_linesearch_max_it 5"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_linesearch_maxstep 1.0"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_linesearch_minlambda 0.05"));

		}

		if(opt.use_eisenstat_walker)
		{
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_version 3"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_rtol0 1e-2"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_rtolmax 1e-2"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_gamma 0.9"));
			PetscCall(PetscOptionsInsertString(NULL, "-snes_ksp_ew_alpha 2.0"));
		}

		PetscCall(set_scalar_option ("snes_picard_rtol",            opt.picard_to_newton[0]));
		PetscCall(set_integer_option("snes_picard_minit", (PetscInt)opt.picard_to_newton[1]));
		PetscCall(set_scalar_option ("snes_newton_rtol",            opt.picard_to_newton[2]));
		PetscCall(set_integer_option("snes_newton_maxit", (PetscInt)opt.picard_to_newton[3]));
	}

	//==============
	// LINEAR SOLVER
	//==============

	PetscCall(PetscOptionsInsertString(NULL, "-js_ksp_type fgmres"));
	PetscCall(PetscOptionsInsertString(NULL, "-js_ksp_converged_reason"));

	if(opt.monitor_solvers)
	{
		PetscCall(PetscOptionsInsertString(NULL, "-js_ksp_monitor"));
	}

	PetscCall(set_tolerances("js_ksp", opt.linear_tolerances));

	if(opt.use_mat_free_jac)
	{
		PetscCall(PetscOptionsInsertString(NULL, "-js_mat_free"));
	}

	//==========
	// MULTIGRID
	//==========

	if((!strcmp(opt.stokes_solver,       "coupled_mg")
	||  !strcmp(opt.stokes_solver,       "block_mg")
	||  !strcmp(opt.stokes_solver,       "wbfbt")
	||  !strcmp(opt.init_thermal_solver, "mg")))
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

		// get 2D coarsening flag
		PetscCall(FDSTAGCheckMG2D(fs, MG2D));

		// check restrictions
		if(MG2D)
		{
			// WARNING! 2D coarsening is not implemented for 3D multigrid in PETSc
			if(!strcmp(opt.stokes_solver,"wbfbt"))
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "wBFBT solver is not available for 2D grid (stokes_solver): %s", opt.stokes_solver);
			}
			if(!strcmp(opt.init_thermal_solver, "mg"))
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Multigrid temperature solver is not available for 2D grid (init_thermal_solver): %s", opt.init_thermal_solver);
			}
		}

		// select number of multigrid levels
		PetscCall(get_num_mg_levels(opt, fs));

		// compute local grid size on all levels
		PetscCall(FDSTAGGetLevelsLocalGridSize(fs, opt.num_mg_levels,
				levels_num_local_cells, coarse_num_local_cells));

		// select coarse solve reduction factor
		PetscCall(get_coarse_reduction_factor(opt, coarse_num_local_cells));

		// get number of local blocks per processor
		PetscCall(get_num_local_blocks(opt, levels_num_local_cells, coarse_num_local_cells));

		// destroy grid
		PetscCall(FDSTAGDestroy(fs));

		// select option for scalable triple matrix product (required for a)
		PetscCall(PetscOptionsInsertString(NULL, "-matmatmatmult_via scalable"));
	}

	//===============
	// PRECONDITIONER
	//===============

	if(!strcmp(opt.stokes_solver, "block_direct"))
	{
		PetscCall(set_string_option("jp_type",                      "bf"));
		PetscCall(set_scalar_option("jp_pgamma",                    opt.penalty));
		PetscCall(set_string_option("bf_vs_type",                   "user"));
		PetscCall(set_string_option("vs_ksp_type",                  "preonly"));
		PetscCall(set_string_option("vs_pc_type",                   "lu"));
		PetscCall(set_string_option("vs_pc_factor_mat_solver_type", opt.direct_solver_type));
	}
	else if(!strcmp(opt.stokes_solver, "coupled_mg"))
	{
		PetscCall(set_string_option ("jp_type", "mg"));

		if(opt.num_mat_free_levels)
		{
			PetscCall(set_integer_option("gmg_mat_free_levels", opt.num_mat_free_levels));
		}

		PetscCall(set_custom_mg_options(opt, "gmg"));
	}
	else if(!strcmp(opt.stokes_solver, "block_mg"))
	{
		PetscCall(set_string_option("jp_type",       "bf"));
		PetscCall(set_string_option("bf_vs_type",    "mg"));
		PetscCall(set_string_option("bf_schur_type", "inv_eta"));
		PetscCall(set_string_option("vs_ksp_type",   "preonly"));

		PetscCall(set_custom_mg_options(opt, "gmg"));
	}
	else if(!strcmp(opt.stokes_solver, "wbfbt"))
	{
		PetscCall(set_string_option("jp_type",       "bf"));
		PetscCall(set_string_option("bf_vs_type",    "mg"));
		PetscCall(set_string_option("bf_schur_type", "wbfbt"));
		PetscCall(set_string_option("vs_ksp_type",   "preonly"));
		PetscCall(set_string_option("ks_ksp_type",   "preonly"));

		PetscCall(set_custom_mg_options(opt, "gmg"));

		PetscCall(set_standard_mg_options(opt, "ks"));
	}

	//===============
	// THERMAL SOLVER
	//===============

	PetscCall(getIntParam(fb, _OPTIONAL_, "act_temp_diff",   &act_temp_diff,   1, 1));
	PetscCall(getIntParam(fb, _OPTIONAL_, "act_steady_temp", &act_steady_temp, 1, 1));

	// transient solver
	if(act_temp_diff)
	{
		PetscCall(set_string_option("ksp_type", "gmres", "ts"));

		PetscCall(set_tolerances("ts_ksp", opt.thermal_tolerances));

		if(opt.monitor_solvers) { PetscCall(set_empty_option  ("ksp_monitor",  "ts")); }
		if(opt.view_solvers)    { PetscCall(set_empty_option  ("ksp_view",     "ts")); }
	}

	if(act_steady_temp)
	{
		PetscCall(set_string_option ("ksp_type",             "gmres",                   "its"));
		PetscCall(set_scalar_option ("ksp_rtol",             opt.thermal_tolerances[0], "its"));
		PetscCall(set_integer_option("ksp_max_it", (PetscInt)opt.thermal_tolerances[2], "its"));

		if(!strcmp(opt.init_thermal_solver, "mg"))
		{
			PetscCall(set_standard_mg_options(opt, "its"));
		}

		if(opt.monitor_solvers) { PetscCall(set_empty_option  ("ksp_monitor",  "its")); }
		if(opt.view_solvers)    { PetscCall(set_empty_option  ("ksp_view",     "its")); }
	}


	set_empty_option("options_left");

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode solverOptionsReadFromFile(FB *fb, SolOptDB &opt)
{
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
		opt.skip_defaults =  0;

		PetscCall(getIntParam(fb, _OPTIONAL_, "skip_defaults", &opt.skip_defaults, 1, 1));

		if(opt.skip_defaults)
		{
			PetscCall(FBFreeBlocks(fb));

			PetscFunctionReturn(0);
		}

		// read simplified solver options
		PetscCall(getIntParam   (fb, _OPTIONAL_, "view_solvers",            &opt.view_solvers,            1, 1));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "monitor_solvers",         &opt.monitor_solvers,         1, 1));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "set_linear_problem",      &opt.set_linear_problem,      1, 1));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "nonlinear_tolerances",     opt.nonlinear_tolerances,    3, 1.0));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "linear_tolerances",        opt.linear_tolerances,       3, 1.0));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "picard_to_newton",         opt.picard_to_newton,        4, 1.0));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "use_line_search",         &opt.use_line_search,         1, 1));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "use_eisenstat_walker",    &opt.use_eisenstat_walker,    1, 1));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "use_mat_free_jac",        &opt.use_mat_free_jac,        1, 1));
		PetscCall(getStringParam(fb, _OPTIONAL_, "stokes_solver",            opt.stokes_solver,           NULL));
		PetscCall(getStringParam(fb, _OPTIONAL_, "direct_solver_type",       opt.direct_solver_type,      NULL));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "penalty",                 &opt.penalty,                 1, 1.0));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "num_mg_levels",           &opt.num_mg_levels,           1, _max_num_mg_levels_));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "num_mat_free_levels",     &opt.num_mat_free_levels,     1, _max_num_mat_free_levels_));
		PetscCall(getStringParam(fb, _OPTIONAL_, "smoother_type",            opt.smoother_type,           NULL));
		PetscCall(getStringParam(fb, _OPTIONAL_, "smoother_ksp",             opt.smoother_ksp,            NULL));
		PetscCall(getStringParam(fb, _OPTIONAL_, "smoother_pc",              opt.smoother_pc,             NULL));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "smoother_damping",        &opt.smoother_damping,        1, 1.0));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "smoother_omega",          &opt.smoother_omega,          1, 1.0));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "smoother_num_sweeps",     &opt.smoother_num_sweeps,     1, 1000));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "coarse_reduction_factor", &opt.coarse_reduction_factor, 1, 1024));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "coarse_cells_per_cpu",    &opt.coarse_cells_per_cpu,    1, 32768));
		PetscCall(getStringParam(fb, _OPTIONAL_, "coarse_solver",            opt.coarse_solver,           NULL));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "coarse_tolerances",        opt.coarse_tolerances,       2, 1.0));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "subdomain_overlap",       &opt.subdomain_overlap,       1, 10));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "subdomain_ilu_levels",    &opt.subdomain_ilu_levels,    1, 8));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "subdomain_num_cells",     &opt.subdomain_num_cells,     1, 32768));
		PetscCall(getStringParam(fb, _OPTIONAL_, "init_thermal_solver",      opt.init_thermal_solver,     NULL));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "thermal_tolerances",       opt.thermal_tolerances,      3, 1.0));
	}

	// clear options block
	PetscCall(FBFreeBlocks(fb));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode solverOptionsCheck(SolOptDB &opt)
{
	PetscFunctionBeginUser;

	// check parameters
	if(!(!strcmp(opt.stokes_solver, "block_direct")
	||   !strcmp(opt.stokes_solver, "coupled_mg")
	||   !strcmp(opt.stokes_solver, "block_mg")
	||   !strcmp(opt.stokes_solver, "wbfbt")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect Stokes solver type (stokes_solver): %s", opt.stokes_solver);
	}

	if(!(!strcmp(opt.direct_solver_type, "superlu_dist")
	||   !strcmp(opt.direct_solver_type, "mumps")
	||   !strcmp(opt.direct_solver_type, "lu")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect direct solver type (direct_solver_type): %s", opt.direct_solver_type);
	}

	if(!(!strcmp(opt.smoother_type, "light")
	||   !strcmp(opt.smoother_type, "intermediate")
	||   !strcmp(opt.smoother_type, "heavy")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect smoother type (smoother_type): %s", opt.smoother_type);
	}

	// set default smoother if not specified
	PetscCall(set_default_smoother(opt));

	if(!(!strcmp(opt.smoother_ksp, "richardson")
	||   !strcmp(opt.smoother_ksp, "chebyshev")
	||   !strcmp(opt.smoother_ksp, "gmres")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect smoother solver type (smoother_ksp): %s", opt.smoother_ksp);
	}

	if(!(!strcmp(opt.smoother_pc, "jacobi")
    ||   !strcmp(opt.smoother_pc, "sor")
	||   !strcmp(opt.smoother_pc, "bjacobi")
	||   !strcmp(opt.smoother_pc, "asm")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect smoother preconditioner type (smoother_pc): %s", opt.smoother_pc);
	}

	if(!(!strcmp(opt.coarse_solver, "direct")
	||   !strcmp(opt.coarse_solver, "hypre")
	||   !strcmp(opt.coarse_solver, "bjacobi")
	||   !strcmp(opt.coarse_solver, "asm")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect coarse solver type (coarse_solver): %s", opt.coarse_solver);
	}

	if(!(!strcmp(opt.init_thermal_solver, "mg")
	||   !strcmp(opt.init_thermal_solver, "default")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect initial thermal solver type (init_thermal_solver): %s", opt.init_thermal_solver);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode get_num_mg_levels(SolOptDB &opt, FDSTAG *fs)
{

	// select number of multigrid levels

	PetscInt ncors;

	PetscFunctionBeginUser;

	// get maximum possible number of coarsening steps
	PetscCall(FDSTAGCheckMG(fs, ncors));

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "General multigrid settings      :\n");

	if(opt.num_mg_levels != -1)
	{
		// check user-specified number of levels
		if(opt.num_mg_levels < 2 || opt.num_mg_levels > ncors + 1)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect # of multigrid levels specified. Requested: %lld. Max. possible: %lld", (LLD)opt.num_mg_levels, (LLD)(ncors + 1));
		}

		PetscPrintf(PETSC_COMM_WORLD, "   Number of multigrid levels   : %lld (prescribed)\n", (LLD)opt.num_mg_levels);
	}
	else
	{
		//==========================================================================================
		// WARNING!
		// for small problems try to respect "coarse_cells_per_cpu" constraint here (use less levels)
		//==========================================================================================

		opt.num_mg_levels = ncors + 1;

		PetscPrintf(PETSC_COMM_WORLD, "   Number of multigrid levels   : %lld (automatic)\n", (LLD)opt.num_mg_levels);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode get_num_local_blocks(
		SolOptDB &opt,
		PetscInt  levels_num_local_cells[],
		PetscInt  coarse_num_local_cells)
{
	PetscInt i, ncells, petsc_mg_level;

	PetscFunctionBeginUser;

	if(opt.subdomain_num_cells != -1)
	{
		// levels
		for(i = 0; i < opt.num_mg_levels - 1; i++)
		{
			opt.levels_num_local_blocks[i] = PetscCeilInt(levels_num_local_cells[i], opt.subdomain_num_cells);
		}

		// update number of local cells per aggregated cpu
		ncells = coarse_num_local_cells*opt.coarse_reduction_factor;

		//	coarse grid
		opt.coarse_num_local_blocks = PetscCeilInt(ncells, opt.subdomain_num_cells);
	}
	else
	{
		// levels
		for(i = 0; i < opt.num_mg_levels - 1; i++)
		{
			opt.levels_num_local_blocks[i] = 1;
		}

		//	coarse grid
		opt.coarse_num_local_blocks = 1;
	}

	// check constant number of blocks on all levels
	opt.levels_num_blocks_constant = opt.levels_num_local_blocks[0];

	for(i = 1; i < opt.num_mg_levels - 1; i++)
	{
		if(opt.levels_num_local_blocks[i] != opt.levels_num_blocks_constant)
		{
			opt.levels_num_blocks_constant = 0;

			break;
		}
	}

	if(opt.levels_num_blocks_constant)
	{
		PetscPrintf(PETSC_COMM_WORLD, "   Number of blocks on levels   : %lld\n", (LLD)opt.levels_num_blocks_constant);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "   Number of blocks on          :\n");

		for(i = 0, petsc_mg_level = opt.num_mg_levels-1; i < opt.num_mg_levels - 1; i++, petsc_mg_level--)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      PETSc level               : %lld -> %lld\n", (LLD)petsc_mg_level, (LLD)opt.levels_num_local_blocks[i]);
		}
	}

	PetscPrintf(PETSC_COMM_WORLD, "   Number of blocks coarse grid : %lld\n", (LLD)opt.coarse_num_local_blocks);

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode get_coarse_reduction_factor(
		SolOptDB &opt,
		PetscInt  coarse_num_local_cells)
{
	// get number of processors for coarse grid solve

	PetscMPIInt size;
	PetscInt    total_num_cpu, coarse_num_cpu, targer_factor, lower_factor, upper_factor;

	PetscFunctionBeginUser;

	// get number of ranks
	MPI_Comm_size(PETSC_COMM_WORLD, &size);

	total_num_cpu = (PetscInt)size;

	if(opt.coarse_reduction_factor != -1)
	{
		// check user-specified reduction factor
		if(total_num_cpu % opt.coarse_reduction_factor)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect reduction factor specified (coarse_reduction_factor): %lld", (LLD)opt.coarse_reduction_factor);
		}

		PetscPrintf(PETSC_COMM_WORLD, "   Coarse grid reduction factor : %lld (prescribed)\n", (LLD)opt.coarse_reduction_factor);
	}
	else if(total_num_cpu == 1)
	{
		// sequential case
		opt.coarse_reduction_factor = 1;

		PetscPrintf(PETSC_COMM_WORLD, "   Coarse grid reduction factor : %lld (sequential solve)\n", (LLD)opt.coarse_reduction_factor);
	}
	else if(opt.coarse_cells_per_cpu == -1)
	{
		// all processors are used for coarse solve
		opt.coarse_reduction_factor = 1;

		PetscPrintf(PETSC_COMM_WORLD, "   Coarse grid reduction factor : %lld (default)\n", (LLD)opt.coarse_reduction_factor);
	}
	else
	{
		// compute target number of processors
		coarse_num_cpu = PetscCeilInt(coarse_num_local_cells, opt.coarse_cells_per_cpu);

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
			opt.coarse_reduction_factor = lower_factor;
		}
		else
		{
			opt.coarse_reduction_factor = upper_factor;
		}

		PetscPrintf(PETSC_COMM_WORLD, "   Coarse grid reduction factor : %lld (automatic)\n", (LLD)opt.coarse_reduction_factor);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_default_smoother(SolOptDB &opt)
{
	PetscFunctionBeginUser;

	char smoother_ksp_default[_str_len_] = {'\0'};
	char smoother_pc_default [_str_len_] = {'\0'};

	if(!strcmp(opt.smoother_type, "light"))
	{
		strcpy(smoother_ksp_default, "richardson");
		strcpy(smoother_pc_default,  "jacobi");
	}
	else if(!strcmp(opt.smoother_type, "intermediate"))
	{
		strcpy(smoother_ksp_default, "chebyshev");
		strcpy(smoother_pc_default,  "sor");
	}
	else if(!strcmp(opt.smoother_type, "heavy"))
	{
		strcpy(smoother_ksp_default, "gmres");
		strcpy(smoother_pc_default,  "bjacobi");
	}

	if(!strlen(opt.smoother_ksp)) {  PetscCall(PetscStrncpy(opt.smoother_ksp, smoother_ksp_default, _str_len_)); }
	if(!strlen(opt.smoother_pc))  {  PetscCall(PetscStrncpy(opt.smoother_pc,  smoother_pc_default,  _str_len_)); }

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_smoother_options(
		SolOptDB   &opt,
		const char *prefix,
		PetscInt    num_local_blocks)
{
	PetscFunctionBeginUser;

	//====
	// KSP
	//====

	PetscCall(set_string_option ("ksp_type",   opt.smoother_ksp,        prefix));
	PetscCall(set_integer_option("ksp_max_it", opt.smoother_num_sweeps, prefix));

	if(!strcmp(opt.smoother_ksp, "richardson"))
	{
		PetscCall(set_scalar_option("ksp_richardson_scale", opt.smoother_damping, prefix));
	}
	else if(!strcmp(opt.smoother_ksp, "gmres"))
	{
		PetscCall(set_integer_option("ksp_gmres_restart", opt.smoother_num_sweeps, prefix));
	}

	//===
	// PC
	//===

	PetscCall(set_string_option("pc_type", opt.smoother_pc, prefix));

	if(!strcmp(opt.smoother_pc, "sor"))
	{
		PetscCall(set_scalar_option("pc_sor_omega", opt.smoother_omega, prefix));
	}
	else if(!strcmp(opt.smoother_pc, "bjacobi")
	||      !strcmp(opt.smoother_pc, "asm"))
	{
		PetscCall(set_subdomain_options(opt, prefix, opt.smoother_pc, num_local_blocks));
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_subdomain_options(
		SolOptDB   &opt,
		const char *prefix,
		const char *pc_type,
		PetscInt    num_local_blocks)
{
	PetscFunctionBeginUser;

	if(!strcmp(pc_type, "bjacobi"))
	{
		PetscCall(set_string_option ("pc_type",                "bjacobi",         prefix));
		PetscCall(set_integer_option("pc_bjacobi_local_blocks", num_local_blocks, prefix));
	}
	else if(!strcmp(pc_type, "asm"))
	{
		PetscCall(set_string_option ("pc_type",             "asm",                 prefix));
		PetscCall(set_integer_option("pc_asm_local_blocks", num_local_blocks,      prefix));
		PetscCall(set_integer_option("pc_asm_overlap",      opt.subdomain_overlap, prefix));
		PetscCall(set_string_option ("pc_asm_type",         "restrict",            prefix));
		PetscCall(set_string_option ("pc_asm_local_type",   "additive",            prefix));
	}

	//==========================================================================================
	// WARNING!
	// check wheather natural ordering gives better results compared to nested disection
	// provide better pre-allocation parameter if more then zero levels is selected
	//==========================================================================================

	PetscCall(set_string_option ("sub_ksp_type",                    "preonly",                 prefix));
	PetscCall(set_string_option ("sub_pc_type",                     "ilu",                     prefix));
	PetscCall(set_integer_option("sub_pc_factor_levels",             opt.subdomain_ilu_levels, prefix));
	PetscCall(set_string_option ("sub_pc_factor_mat_ordering_type",  "nd",                     prefix));
	PetscCall(set_empty_option  ("sub_pc_factor_reuse_ordering",                               prefix));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_coarse_options(
		SolOptDB   &opt,
		const char *mg_prefix)
{
	char *prefix;

	PetscFunctionBeginUser;

	// compile coarse solver prefix
	asprintf(&prefix,"%s_mg_coarse", mg_prefix);

	//==========
	// TELESCOPE
	//==========

	if(opt.coarse_reduction_factor > 1)
	{
		PetscCall(set_string_option ("ksp_type",                     "preonly",                    prefix));
		PetscCall(set_string_option ("pc_type",                      "telescope",                  prefix));
		PetscCall(set_integer_option("pc_telescope_reduction_factor", opt.coarse_reduction_factor, prefix));

		free(prefix);

		// compile telescope prefixs
		asprintf(&prefix,"%s_mg_coarse_telescope", mg_prefix);
	}

	//====
	// KSP
	//====

	if(!strcmp(opt.coarse_solver, "direct"))
	{
		PetscCall(set_string_option("ksp_type", "preonly", prefix));
	}
	else if((!strcmp(opt.coarse_solver, "hypre")
	||       !strcmp(opt.coarse_solver, "bjacobi")
	||       !strcmp(opt.coarse_solver, "asm")))
	{
		PetscCall(set_string_option ("ksp_type",             "gmres",                  prefix));
		PetscCall(set_scalar_option ("ksp_rtol",             opt.coarse_tolerances[0], prefix));
		PetscCall(set_integer_option("ksp_max_it", (PetscInt)opt.coarse_tolerances[1], prefix));
	}

	//===
	// PC
	//===

	if(!strcmp(opt.coarse_solver, "direct"))
	{
		PetscCall(set_string_option("pc_type ",                  "lu",                    prefix));
		PetscCall(set_string_option("pc_factor_mat_solver_type ", opt.direct_solver_type, prefix));
	}
	else if(!strcmp(opt.coarse_solver, "hypre"))
	{
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
	}
	else if(!strcmp(opt.coarse_solver, "bjacobi")
	||      !strcmp(opt.coarse_solver, "asm"))
	{
		PetscCall(set_subdomain_options(opt, prefix, opt.coarse_solver, opt.coarse_num_local_blocks));
	}

	free(prefix);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_levels_options(
		SolOptDB   &opt,
		const char *mg_prefix)
{
	char    *prefix;
	PetscInt i, petsc_mg_level;

	PetscFunctionBeginUser;

	if((!strcmp(opt.smoother_pc, "bjacobi")
	||  !strcmp(opt.smoother_pc, "asm"))
	&&  !opt.levels_num_blocks_constant)
	{
		for(i = 0, petsc_mg_level = opt.num_mg_levels-1; i < opt.num_mg_levels - 1; i++, petsc_mg_level--)
		{
			// compile level prefix
			asprintf(&prefix,"%s_mg_levels_%lld", mg_prefix, (LLD)petsc_mg_level);

			// set smoother options level-wise (different number of local subdomains (blocks) per processor)
			PetscCall(set_smoother_options(opt, prefix, opt.levels_num_local_blocks[i]));

			free(prefix);
		}
	}
	else
	{
		asprintf(&prefix,"%s_mg_levels", mg_prefix);

		// set same options for all levels
		PetscCall(set_smoother_options(opt, prefix, opt.levels_num_blocks_constant));

		free(prefix);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_custom_mg_options(
		SolOptDB   &opt,
		const char *prefix)
{
	PetscFunctionBeginUser;

	if(opt.view_solvers)
	{
		PetscCall(set_empty_option("pc_view", prefix));
	}

	// number of multigrid levels
	PetscCall(set_integer_option("pc_mg_levels", opt.num_mg_levels, prefix));

	// setup coarse solver
	PetscCall(set_coarse_options(opt, prefix));

	// setup level smoothers
	PetscCall(set_levels_options(opt, prefix));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode set_standard_mg_options(SolOptDB &opt, const char *prefix)
{
	PetscFunctionBeginUser;

	PetscCall(set_string_option ("pc_type", "mg",                   prefix));
	PetscCall(set_integer_option("pc_mg_levels", opt.num_mg_levels, prefix));
	PetscCall(set_empty_option  ("pc_mg_galerkin",                  prefix));
	PetscCall(set_string_option ("pc_mg_type", "multiplicative",    prefix));
	PetscCall(set_string_option ("pc_mg_cycle_type", "v",           prefix));

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
PetscErrorCode PetscOptionsReadFromFile(FB *fb)
{
	// * load additional options from input file
	// * push command line options to the end of database
	// (PETSc prioritizes options appearing LAST)

	PetscInt  jj, i, lnbeg, lnend;
	char     *line, **lines, *key, *val, *option;

	PetscFunctionBeginUser;

	// setup block access mode
	PetscCall(FBFindBlocks(fb, _OPTIONAL_, "<PetscOptionsStart>", "<PetscOptionsEnd>"));

	if(fb->nblocks > 1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many petsc options blocks. Only one is allowed");
	}

	// get line buffer
	line = fb->lbuf;

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		lines = FBGetLineRanges(fb, &lnbeg, &lnend);

		for(i = lnbeg; i < lnend; i++)
		{
			// copy line for parsing
			strcpy(line, lines[i]);

			// get key
			key = strtok(line, " ");

			if(!key) continue;

			// get value
			val = strtok(NULL, " ");

			if(!val) option = key;
			else     asprintf(&option, "%s %s", key, val);

			PetscCall(PetscOptionsInsertString(NULL, option));

			if(val) free(option);
		}

		fb->blockID++;
	}

	PetscCall(FBFreeBlocks(fb));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------


