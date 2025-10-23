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
//.....................   MATRIX AUXILIARY FUNCTIONS   ......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "matData.h"
#include "matrix.h"
//---------------------------------------------------------------------------
// MATRIX SERVICE FUNCTIONS
//---------------------------------------------------------------------------
PetscErrorCode MatAIJCreate(
	PetscInt m, PetscInt n,
	PetscInt d_nz, const PetscInt d_nnz[],
	PetscInt o_nz, const PetscInt o_nnz[],
	Mat *P)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// create matrix
	ierr = MatCreate(PETSC_COMM_WORLD, P); CHKERRQ(ierr);
	ierr = MatSetType((*P), MATAIJ); CHKERRQ(ierr);
	ierr = MatSetSizes((*P), m, n, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);

	// preallocate matrix
	ierr = MatSeqAIJSetPreallocation((*P), d_nz, d_nnz); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation((*P), d_nz, d_nnz, o_nz, o_nnz); CHKERRQ(ierr);

	// read custom options (required to resolve SuperLU_DIST issue)
	ierr = MatSetFromOptions((*P)); CHKERRQ(ierr);

	// throw an error if preallocation fails
	ierr = MatSetOption((*P), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetUp((*P)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatAIJCreateDiag(PetscInt m, PetscInt istart, Mat *P)
{
	PetscInt i, ii;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// preallocate
	ierr = MatAIJCreate(m, m, 1, NULL, 0, NULL, P); CHKERRQ(ierr);

	// put explicit zeroes on the diagonal
	for(i = 0; i < m; i++)
	{
		// get global row index
		ii = istart + i;

		ierr = MatSetValue((*P), ii, ii, 0.0, INSERT_VALUES); CHKERRQ(ierr);
	}

	ierr = MatSetFromOptions((*P)); CHKERRQ(ierr);

	// assemble
	ierr = MatAIJAssemble((*P), 0, NULL, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatAIJAssemble(Mat P, PetscInt numRows, const PetscInt rows[], PetscScalar diag)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = MatSetOption(P, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// freeze nonzero structure, constrain rows only locally
	ierr = MatSetOption(P, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(P, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);     CHKERRQ(ierr);
	ierr = MatSetOption(P, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);    CHKERRQ(ierr);

	// zero out constrained rows, form unit diagonal for the constrained block
	ierr = MatZeroRows(P, numRows, rows, diag, NULL, NULL); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode MatAIJSetNullSpace(Mat P, MatData *md)
{
	DOFIndex     *dof;
	MatNullSpace nullsp;                       // near null space
	Vec          nullsp_vecs[_max_nullsp_sz_]; // near null space vectors
	PetscScalar *v;
	PetscInt     i, j, sz, ln, iter, nullsp_sz, lbsz[_max_nullsp_sz_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access context
	dof = &md->fs->dof;

	// get number of vectors
	if     (md->idxmod == _IDX_COUPLED_) { nullsp_sz = 4; ln = dof->lnv + dof->lnp; }
	else if(md->idxmod == _IDX_BLOCK_)   { nullsp_sz = 3; ln = dof->lnv; }

	// set local block sizes & iterator
	iter    = 0;
	lbsz[0] = dof->lnvx;
	lbsz[1] = dof->lnvy;
	lbsz[2] = dof->lnvz;
	lbsz[3] = dof->lnp;

	// create near null space vectors
	for(i = 0; i < nullsp_sz; i++)
	{
		// create
		ierr = VecCreateMPI(PETSC_COMM_WORLD, ln, PETSC_DETERMINE, &nullsp_vecs[i]); CHKERRQ(ierr);
		ierr = VecSetFromOptions(nullsp_vecs[i]);   CHKERRQ(ierr);
		ierr = VecZeroEntries   (nullsp_vecs[i]);   CHKERRQ(ierr);

		// initialize
		ierr = VecZeroEntries (nullsp_vecs[i]);     CHKERRQ(ierr);
		ierr = VecGetArray    (nullsp_vecs[i], &v); CHKERRQ(ierr);

		for(j = 0, sz = lbsz[i]; j < sz; j++) v[iter++] = 1.0;

		ierr = VecRestoreArray(nullsp_vecs[i], &v); CHKERRQ(ierr);

		// normalize
		ierr = VecNormalize(nullsp_vecs[i], NULL); CHKERRQ(ierr);
	}

	// create near null space
	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, nullsp_sz, (const Vec*)nullsp_vecs, &nullsp); CHKERRQ(ierr);

	// attach near null space to the matrix
	ierr = MatSetNearNullSpace(P, nullsp); CHKERRQ(ierr);

	// clear storage
	ierr = MatNullSpaceDestroy(&nullsp); CHKERRQ(ierr);

	for(i = 0; i < nullsp_sz; i++)
	{
		ierr = VecDestroy(&nullsp_vecs[i]); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void getStiffMat(
	PetscScalar eta, PetscScalar diag,
	PetscScalar *v,  PetscScalar *cf,
	PetscScalar dx,  PetscScalar dy,   PetscScalar dz,
	PetscScalar fdx, PetscScalar fdy,  PetscScalar fdz,
	PetscScalar bdx, PetscScalar bdy,  PetscScalar bdz)
{
	// compute cell stiffness matrix with deviatoric projection

	PetscScalar E43 = 4.0*eta/3.0;
	PetscScalar E23 = 2.0*eta/3.0;

	//       vx_(i)               vx_(i+1)             vy_(j)               vy_(j+1)             vz_(k)               vz_(k+1)             p
	v[0]  =  E43/dx/bdx; v[1]  = -E43/dx/bdx; v[2]  = -E23/dy/bdx; v[3]  =  E23/dy/bdx; v[4]  = -E23/dz/bdx; v[5]  =  E23/dz/bdx; v[6]  =  cf[0]/bdx; // fx_(i)   [sxx]
	v[7]  = -E43/dx/fdx; v[8]  =  E43/dx/fdx; v[9]  =  E23/dy/fdx; v[10] = -E23/dy/fdx; v[11] =  E23/dz/fdx; v[12] = -E23/dz/fdx; v[13] = -cf[1]/fdx; // fx_(i+1) [sxx]
	v[14] = -E23/dx/bdy; v[15] =  E23/dx/bdy; v[16] =  E43/dy/bdy; v[17] = -E43/dy/bdy; v[18] = -E23/dz/bdy; v[19] =  E23/dz/bdy; v[20] =  cf[2]/bdy; // fy_(j)   [syy]
	v[21] =  E23/dx/fdy; v[22] = -E23/dx/fdy; v[23] = -E43/dy/fdy; v[24] =  E43/dy/fdy; v[25] =  E23/dz/fdy; v[26] = -E23/dz/fdy; v[27] = -cf[3]/fdy; // fy_(j+1) [syy]
	v[28] = -E23/dx/bdz; v[29] =  E23/dx/bdz; v[30] = -E23/dy/bdz; v[31] =  E23/dy/bdz; v[32] =  E43/dz/bdz; v[33] = -E43/dz/bdz; v[34] =  cf[4]/bdz; // fz_(k)   [szz]
	v[35] =  E23/dx/fdz; v[36] = -E23/dx/fdz; v[37] =  E23/dy/fdz; v[38] = -E23/dy/fdz; v[39] = -E43/dz/fdz; v[40] =  E43/dz/fdz; v[41] = -cf[5]/fdz; // fz_(k+1) [szz]
	v[42] =  1.0/dx;     v[43] = -1.0/dx;     v[44] =  1.0/dy;     v[45] = -1.0/dy;     v[46] =  1.0/dz;     v[47] = -1.0/dz;     v[48] =  diag;      // g
}
//---------------------------------------------------------------------------
void addDensGradStabil(
	PetscScalar fssa, PetscScalar *v,
	PetscScalar rho,  PetscScalar dt,   PetscScalar *grav,
	PetscScalar fdx,  PetscScalar fdy,  PetscScalar fdz,
	PetscScalar bdx,  PetscScalar bdy,  PetscScalar bdz)
{
	PetscScalar cf = -fssa*dt;

	// add stabilization terms
	v[0 ] -= cf*(rho*grav[0])/bdx;
	v[8 ] += cf*(rho*grav[0])/fdx;
	v[16] -= cf*(rho*grav[1])/bdy;
	v[24] += cf*(rho*grav[1])/fdy;
	v[32] -= cf*(rho*grav[2])/bdz;
	v[40] += cf*(rho*grav[2])/fdz;
}
//---------------------------------------------------------------------------
void getVelSchur(PetscScalar v[], PetscScalar d[], PetscScalar g[])
{
	PetscScalar k;

	// extract divergence operator
	d[0] = v[42]; d[1] = v[43]; d[2] = v[44]; d[3] = v[45]; d[4] = v[46]; d[5] = v[47];

	// extract gradient operator
	g[0] = v[6];  g[1] = v[13]; g[2] = v[20]; g[3] = v[27]; g[4] = v[34]; g[5] = v[41];

	// get penalty term
	k = -1.0/v[48];

	// compute velocity Schur complement
	v[0]  += k*g[0]*d[0]; v[1]  += k*g[0]*d[1]; v[2]  += k*g[0]*d[2]; v[3]  += k*g[0]*d[3]; v[4]  += k*g[0]*d[4]; v[5]  += k*g[0]*d[5];
	v[7]  += k*g[1]*d[0]; v[8]  += k*g[1]*d[1]; v[9]  += k*g[1]*d[2]; v[10] += k*g[1]*d[3]; v[11] += k*g[1]*d[4]; v[12] += k*g[1]*d[5];
	v[14] += k*g[2]*d[0]; v[15] += k*g[2]*d[1]; v[16] += k*g[2]*d[2]; v[17] += k*g[2]*d[3]; v[18] += k*g[2]*d[4]; v[19] += k*g[2]*d[5];
	v[21] += k*g[3]*d[0]; v[22] += k*g[3]*d[1]; v[23] += k*g[3]*d[2]; v[24] += k*g[3]*d[3]; v[25] += k*g[3]*d[4]; v[26] += k*g[3]*d[5];
	v[28] += k*g[4]*d[0]; v[29] += k*g[4]*d[1]; v[30] += k*g[4]*d[2]; v[31] += k*g[4]*d[3]; v[32] += k*g[4]*d[4]; v[33] += k*g[4]*d[5];
	v[35] += k*g[5]*d[0]; v[36] += k*g[5]*d[1]; v[37] += k*g[5]*d[2]; v[38] += k*g[5]*d[3]; v[39] += k*g[5]*d[4]; v[40] += k*g[5]*d[5];
}
//---------------------------------------------------------------------------
void getSubMat(PetscScalar v[],  PetscScalar a[], PetscScalar d[], PetscScalar g[])
{
	// extract divergence operator
	d[0] = v[42]; d[1] = v[43]; d[2] = v[44]; d[3] = v[45]; d[4] = v[46]; d[5] = v[47];

	// extract gradient operator
	g[0] = v[6];  g[1] = v[13]; g[2] = v[20]; g[3] = v[27]; g[4] = v[34]; g[5] = v[41];

	// extract velocity operator
	a[0]  = v[0];  a[1]  = v[1];  a[2]  = v[2];  a[3]  = v[3];  a[4]  = v[4];  a[5]  = v[5];
	a[6]  = v[7];  a[7]  = v[8];  a[8]  = v[9];  a[9]  = v[10]; a[10] = v[11]; a[11] = v[12];
	a[12] = v[14]; a[13] = v[15]; a[14] = v[16]; a[15] = v[17]; a[16] = v[18]; a[17] = v[19];
	a[18] = v[21]; a[19] = v[22]; a[20] = v[23]; a[21] = v[24]; a[22] = v[25]; a[23] = v[26];
	a[24] = v[28]; a[25] = v[29]; a[26] = v[30]; a[27] = v[31]; a[28] = v[32]; a[29] = v[33];
	a[30] = v[35]; a[31] = v[36]; a[32] = v[37]; a[33] = v[38]; a[34] = v[39]; a[35] = v[40];
}
//---------------------------------------------------------------------------
void getTwoPointConstr(PetscInt n, PetscInt idx[], PetscInt pdofidx[], PetscScalar cf[])
{
	// apply two-point constraints on the ghost nodes
	PetscInt j;

	for(j = 0; j < n; j++)
	{
		// detect ghost point (the one that has negative global index)
		if(idx[j] == -1)
		{
			// ghost point detected, check whether potential primary DOF is constrained
			if(cf[pdofidx[j]] != DBL_MAX)
			{
				// primary DOF is constrained, put single-point constraint on the ghost node
				cf[j]      =  0.0;
				pdofidx[j] = -1;
			}
			else
			{
				// primary DOF is unconstrained (active), check constraint type
				if(cf[j] != DBL_MAX) cf[j] = -1.0; // no-slip (arbitrary)
				else                 cf[j] =  1.0; // free-slip
			}
		}
		else
		{	// internal point detected, cannot have two-point constraints
			pdofidx[j] = -1;
		}
	}
}
//---------------------------------------------------------------------------
void constrLocalMat(PetscInt n, PetscInt pdofidx[], PetscScalar cf[], PetscScalar v[])
{
	//=========================================================================
	// Apply linear constraints to a local matrix
	//
	// Function supports:
	//  - single-point constraints (x = c)
	//  - two-point constraints    (x = cf*y + c)
	// where:
	//    x   - constrained DOF
	//    y   - primary DOF
	//    cf  - linear combination parameter
	//    c   - constant prescribed value (may be zero for homogeneous case)
	//
	// For the two-point constraints it must be guaranteed that the primary
	// DOF is an active (unconstrained) DOF. Otherwise the two-point constraint
	// simply degenerates to a single-point constraint.
	//
	// Active (unconstrained) DOF are marked by DBL_MAX constant:
	//
	// cf[j] != DBL_MAX:
	//    j is a constrained DOF
	//    cf[j]      - linear combination parameter (not used by single-point constraint)
	//    pdofidx[j] - relative index of the primary DOF (-1 for single-point constraint)
	//
	// otherwise:
	//    j is an unconstrained (active) DOF
	//=========================================================================

	PetscInt i, j, jj, jst;

	for(i = 0, jj = 0; i < n; i++)
	{
		// skip constrained rows (handled by MatZeroRows after assembly)
		if(cf[i] != DBL_MAX) { jj += n; continue; }

		// store pointer to row beginning
		jst = jj;

		// proceed with unconstrained row, check columns
		for(j = 0; j < n; j++, jj++)
		{
			// detect constrained columns
			if(cf[j] != DBL_MAX)
			{
				// compute linear combination of primary & constrained columns
				// NOTE: two-point constraints only!
				if(pdofidx[j] != -1)
				{
					v[jst + pdofidx[j]] += cf[j]*v[jj];
				}

				// zero out constrained column
				v[jj] = 0.0;
			}
		}
	}
}
//---------------------------------------------------------------------------
PetscErrorCode VecScatterBlockToMonolithic(Vec f, Vec g, Vec b, ScatterMode mode)
{
	// scatter block vectors to monolithic format forward & reverse

	PetscInt     fs,  gs,  bs;
	PetscScalar *fp, *gp;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// get sizes of the blocks
	ierr = VecGetLocalSize(f, &fs); CHKERRQ(ierr);
	ierr = VecGetLocalSize(g, &gs); CHKERRQ(ierr);
	ierr = VecGetLocalSize(b, &bs); CHKERRQ(ierr);

	if(bs != fs+gs)
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Block sizes don't match monolithic format");
	}

	// access vectors
	ierr = VecGetArray(f, &fp); CHKERRQ(ierr);
	ierr = VecGetArray(g, &gp); CHKERRQ(ierr);

	if(mode == SCATTER_FORWARD)
	{
		PetscScalar *bp;

		ierr = VecGetArray(b, &bp); CHKERRQ(ierr);

		// block-to-monolithic
		ierr = PetscMemcpy(bp,    fp, (size_t)fs*sizeof(PetscScalar)); CHKERRQ(ierr);
		ierr = PetscMemcpy(bp+fs, gp, (size_t)gs*sizeof(PetscScalar)); CHKERRQ(ierr);

		ierr = VecRestoreArray(b, &bp); CHKERRQ(ierr);

	}
	if(mode == SCATTER_REVERSE)
	{
		const PetscScalar *bp;

		ierr = VecGetArrayRead(b, &bp); CHKERRQ(ierr);

		// monolithic-to-block
		ierr = PetscMemcpy(fp, bp,    (size_t)fs*sizeof(PetscScalar)); CHKERRQ(ierr);
		ierr = PetscMemcpy(gp, bp+fs, (size_t)gs*sizeof(PetscScalar)); CHKERRQ(ierr);

		ierr = VecRestoreArrayRead(b, &bp); CHKERRQ(ierr);

	}

	// restore access
	ierr = VecRestoreArray(f, &fp); CHKERRQ(ierr);
	ierr = VecRestoreArray(g, &gp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


