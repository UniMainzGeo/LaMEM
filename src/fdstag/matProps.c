//---------------------------------------------------------------------------
//......................... MATERIAL INITIALIZATION .........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "solVar.h"
#include "fdstag.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "Parsing.h" // dependency on 'is_comment_line' and 'key_matches'
#include "matProps.h"
//---------------------------------------------------------------------------
// initialize MATERIAL PARAMETERS from file
#undef __FUNCT__
#define __FUNCT__ "MatPropInit"
PetscErrorCode MatPropInit(JacRes *jr, UserCtx *usr)
{
	FILE        *fp;
	Material_t   m;
	PetscInt     i, ils, ile, count, count1;
	PetscInt    *ls,*le;
	PetscInt    *check_phase;
	PetscBool    empty_phase;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// print overview of material parameters read from file
	PetscPrintf(PETSC_COMM_WORLD,"Phase material parameters read from %s: \n",usr->ParamFile);

	// allocate memory for arrays to store line info - with overhead
	ierr = PetscMalloc((size_t)(max_num_phases+_max_overhead_)*sizeof(PetscScalar), &ls); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(max_num_phases+_max_overhead_)*sizeof(PetscScalar), &le); CHKERRQ(ierr);

	// open file, count and get the positions of the material structures in file - check of file was done before
	fp = fopen(usr->ParamFile, "r" );

	// read number of cells
	getLineStruct(fp, ls, le, &count, &count1,"<MaterialStart>","<MaterialEnd>");

	// error checking
	if(count != count1  ) SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incomplete material structures! <MaterialStart>: %lld, <MaterialEnd>: %lld", (LLD)count, (LLD)count1);
	if(count > max_num_phases) SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "More phases found than maximum allowed! Actual: %lld, Max allowed: %lld", (LLD)count, (LLD)max_num_phases);

	// initialize phase checking array
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)count, &check_phase); CHKERRQ(ierr);
	ierr = PetscMemzero(check_phase, sizeof(PetscInt)*(size_t)count); CHKERRQ(ierr);
	for(i = 0; i < count; i++) check_phase[i] = 0;

	// read each individual phase
	for(i = 0; i < count; i++)
	{
		// set position in the file
		ils = ls[i];
		ile = le[i];

		// default values
		MatPropSet(&m);

		// read from file
		ierr = MatPropGetStruct( fp, &m, ils, ile, jr->scal.utype); CHKERRQ(ierr);

		// set material structure
		jr->phases[m.ID] = m;

		// for error check
		check_phase[m.ID] = 1;
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// check empty phases
	empty_phase = PETSC_FALSE;

	for(i = 0; i < count; i++)
	{
		if(!check_phase[i])
		{
			PetscPrintf(PETSC_COMM_WORLD, "Phase %lld is not initialized\n", (LLD)i);
			empty_phase = PETSC_TRUE;
		}
	}

	if(empty_phase == PETSC_TRUE)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect phase numbering");
	}

	jr->numPhases = count;

	// free memory
	fclose(fp);

	// free arrays
	ierr = PetscFree(ls         ); CHKERRQ(ierr);
	ierr = PetscFree(le         ); CHKERRQ(ierr);
	ierr = PetscFree(check_phase); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// set initial material properties
void MatPropSet(Material_t *m)
{
		m->rho                           = 1.0;
		// diffusion creep
		m->Bd                            = 0.0;
		m->Ed                            = 0.0;
		m->Vd                            = 0.0;
		// dislocation creep (power-law)
		m->Bn                            = 0.0;
		m->n                             = 1.0;
		m->En                            = 0.0;
		m->Vn                            = 0.0;
		// Peierls creep
		m->Bp                            = 0.0;
		m->Ep                            = 0.0;
		m->Vp                            = 0.0;
		m->taup                          = 0.0;
		m->gamma                         = 0.0;
		m->q                             = 0.0;
		// elasticity
		m->G                             = 0.0; // will make the model effectively viscous
		m->K                             = 0.0;
		m->Kp                            = 0.0;
		// plasticity
		m->ch                            = 0.0; // effectively switches off plasticity
		m->fr                            = 0.0;     // effectively switches off plasticity
		m->chSoftID                      = -1;     // no law
		m->frSoftID                      = -1;     // no law

		// temperature
		m->alpha                         = 0.0;     // coeff. of thermal expansion
		m->Cp                            = 0.0;     // heat capacity
		m->k                             = 0.0;     // thermal conductivity
		m->A                             = 0.0;     // radiogenic heat production
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPropGetStruct"
PetscErrorCode MatPropGetStruct( FILE *fp, Material_t *m, PetscInt ils, PetscInt ile, UnitsType utype)
{
	PetscInt    err;
	PetscInt    found, found1, found2, found3;
	PetscScalar eta, eta0, e0;

	PetscFunctionBegin;

	// Find material properties in file with error checking
	// id
	getMatPropInt   ( fp, ils, ile, "ID",                   &m->ID,      &found);
	if (!found) {err = 1; MatPropErrorCheck(m->ID,err);}

	// density
	getMatPropScalar( fp, ils, ile, "rho0",                 &m->rho,     &found);
	if (!found) {err = 2; MatPropErrorCheck(m->ID,err);}

	// diffusion creep
	getMatPropScalar( fp, ils, ile, "eta",                  &eta,       &found);
	getMatPropScalar( fp, ils, ile, "Bd",                   &m->Bd,      &found);

	// set diffusion creep constant
	if (!found) m->Bd = 1.0/(2.0*eta);

	// activation energy and volume for diff. creep
	getMatPropScalar( fp, ils, ile, "Ed",                   &m->Ed,      &found);
	getMatPropScalar( fp, ils, ile, "Vd",                   &m->Vd,      &found);

	// power-law (dislocation) creep
	getMatPropScalar( fp, ils, ile, "eta0",                 &eta0,      &found); found1 = found;
	getMatPropScalar( fp, ils, ile, "e0",                   &e0,        &found); found2 = found;
	getMatPropScalar( fp, ils, ile, "Bn",                   &m->Bn,      &found); found3 = found;
	getMatPropScalar( fp, ils, ile, "n",                    &m->n,       &found);

	// error check: some power-law params were specified, but not all important
	if (!found && found1 && found2 && !found3) {err = 40; MatPropErrorCheck(m->ID,err);} // eta0, e0, (Bn, n)
	if (!found &&                      found3) {err = 40; MatPropErrorCheck(m->ID,err);} // Bn, (n)
	if ( found &&!found1 && found2 && !found3) {err = 41; MatPropErrorCheck(m->ID,err);} // n, e0 (eta0, Bn)
	if ( found && found1 &&!found2 && !found3) {err = 41; MatPropErrorCheck(m->ID,err);} // n, eta0 (e0, Bn)
	if ( found &&!found1 &&!found2 && !found3) {err = 41; MatPropErrorCheck(m->ID,err);} // n, (eta0, e0, Bn)

	// convert power-law creep parameters to dislocation creep parameters
	if (!found3 && found1 && found2 && found3) m->Bn = pow(2.0*eta0, -m->n)*pow(e0, 1-m->n);

	// activation energy and volume for disl. creep
	getMatPropScalar( fp, ils, ile, "En",                   &m->En,      &found);
	getMatPropScalar( fp, ils, ile, "Vn",                   &m->Vn,      &found);

	// Peierls creep
	getMatPropScalar( fp, ils, ile, "Bp",                   &m->Bp,      &found); found1 = found;
	getMatPropScalar( fp, ils, ile, "taup",                 &m->taup,    &found);
	if (found && !found1) {err = 5; MatPropErrorCheck(m->ID,err);}
	getMatPropScalar( fp, ils, ile, "gamma",                &m->gamma,   &found);
	if (found && !found1) {err = 5; MatPropErrorCheck(m->ID,err);}
	getMatPropScalar( fp, ils, ile, "q",                    &m->q,       &found);
	if (found && !found1) {err = 5; MatPropErrorCheck(m->ID,err);}

	// activation energy and volume for Peierls creep
	getMatPropScalar( fp, ils, ile, "Ep",                   &m->Ep,      &found);
	if (found && !found1) {err = 5; MatPropErrorCheck(m->ID,err);}
	getMatPropScalar( fp, ils, ile, "Vp",                   &m->Vp,      &found);
	if (found && !found1) {err = 5; MatPropErrorCheck(m->ID,err);}

	// Bd, Bn, Bp = 0
	if ((m->Bd == 0) && (m->Bn == 0) && (m->Bp == 0)) {err = 3; MatPropErrorCheck(m->ID,err);}

	// elasticity
	getMatPropScalar( fp, ils, ile, "shear",                &m->G,       &found);
	getMatPropScalar( fp, ils, ile, "bulk",                 &m->K,       &found);
	getMatPropScalar( fp, ils, ile, "Kp",                   &m->Kp,      &found);

	// plasticity - Drucker-Prager
	getMatPropScalar( fp, ils, ile, "cohesion",             &m->ch,      &found); found1 = found;
	getMatPropScalar( fp, ils, ile, "friction",             &m->fr,      &found);
	if ( found &&!found1) {err = 6; MatPropErrorCheck(m->ID,err);}
	if (!found && found1) {err = 6; MatPropErrorCheck(m->ID,err);}

	// plasticity - strain weakening
	getMatPropInt   ( fp, ils, ile, "chSoftID",            &m->chSoftID,      &found);
	getMatPropInt   ( fp, ils, ile, "frSoftID",            &m->frSoftID,      &found);

	// energy
	getMatPropScalar( fp, ils, ile, "alpha",                &m->alpha,   &found);
	getMatPropScalar( fp, ils, ile, "cp",                   &m->Cp,      &found); found1 = found;
	getMatPropScalar( fp, ils, ile, "k",                    &m->k,       &found);
	if ( found &&!found1) {err = 7; MatPropErrorCheck(m->ID,err);}
	if (!found && found1) {err = 7; MatPropErrorCheck(m->ID,err);}

	getMatPropScalar( fp, ils, ile, "A",                    &m->A,       &found);

	MatPropPrint(m, eta, utype);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// error messages - please see documentation with more details
#undef __FUNCT__
#define __FUNCT__ "MatPropErrorCheck"
PetscErrorCode MatPropErrorCheck(PetscInt id, PetscInt err)
{
	if (err==1 ) SETERRQ (PETSC_COMM_WORLD, PETSC_ERR_USER, "MatPropError(1): No phase ID specified! ");
	if (err==2 ) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatPropError(2): No reference density specified for phase ID = %lld! ",(LLD)id);
	if (err==3 ) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatPropError(3): No creep law is specified for phase ID = %lld! ",(LLD)id);
	if (err==40) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatPropError(40): Power-law: Either (eta0, e0) or (Bn) were specified but not (n) for phase ID = %lld! ",(LLD)id);
	if (err==41) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatPropError(41): Power-law: (n) was specified but not (eta0, e0) or (Bn) for phase ID = %lld! ",(LLD)id);
	if (err==5 ) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatPropError(5): Peierls creep: Bp was not defined, but some other Peierls parameters were for phase ID = %lld! ",(LLD)id);
	if (err==6 ) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatPropError(6): Cohesion and friction angle were not defined together for phase ID = %lld! ",(LLD)id);
	if (err==7 ) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatPropError(7): Cp and k were not defined together for phase ID = %lld! ",(LLD)id);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// print essential material properties
void MatPropPrint(Material_t *m, PetscScalar eta, UnitsType utype)
{
	if      (utype == _NONE_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: rho = %g [-], eta0 = %g [-]\n", (LLD)(m->ID), m->rho, eta);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (diff ) Bd = %g [-], Ed = %g [-], Vd = %g [-] \n", (LLD)(m->ID), m->Bd, m->Ed, m->Vd);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (disl ) Bn = %g [-], En = %g [-], Vn = %g [-], n = %g [-] \n", (LLD)(m->ID), m->Bn, m->En, m->Vn, m->n);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (peirl) Bp = %g [-], Ep = %g [-], Vp = %g [-], taup = %g [-], gamma = %g [-], q = %g [-] \n", (LLD)(m->ID), m->Bp, m->Ep, m->Vp, m->taup, m->gamma, m->q);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (elast) G = %g [-], K = %g [-], Kp = %g [-] \n", (LLD)(m->ID), m->G, m->K, m->Kp);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (plast) cohesion = %g [-], friction angle = %g [-] \n", (LLD)(m->ID),m->ch, m->fr);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (sweak) cohesion SoftLaw = %lld [-], friction SoftLaw = %lld [-] \n", (LLD)(m->ID),(LLD)m->chSoftID, (LLD)m->frSoftID);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (temp ) alpha = %g [-], cp = %g [-], k = %g [-], A = %g [-] \n", (LLD)(m->ID),m->alpha, m->Cp,m->k, m->A);
		PetscPrintf(PETSC_COMM_WORLD,"    \n");
	}
	else if (utype == _SI_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: rho = %g [kg/m3], eta0 = %g [Pa.s]\n", (LLD)(m->ID), m->rho, eta);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (diff ) Bd = %g [1/(Pa.s)], Ed = %g [J/mol], Vd = %g [m3/mol] \n", (LLD)(m->ID), m->Bd, m->Ed, m->Vd);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (disl ) Bn = %g [1/(Pa^n.s)], En = %g [J/mol], Vn = %g [m3/mol], n = %g [-] \n", (LLD)(m->ID), m->Bn, m->En, m->Vn, m->n);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (peirl) Bp = %g [1/s], Ep = %g [J/mol], Vp = %g [m3/mol], taup = %g [Pa], gamma = %g [-], q = %g [-] \n", (LLD)(m->ID), m->Bp, m->Ep, m->Vp, m->taup, m->gamma, m->q);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (elast) G = %g [Pa], K = %g [Pa], Kp = %g [-] \n", (LLD)(m->ID), m->G, m->K, m->Kp);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (plast) cohesion = %g [Pa], friction angle = %g [deg] \n", (LLD)(m->ID),m->ch, m->fr);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (sweak) cohesion SoftLaw = %lld [-], friction SoftLaw = %lld [-] \n", (LLD)(m->ID),(LLD)m->chSoftID, (LLD)m->frSoftID);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (temp ) alpha = %g [1/K], cp = %g [J/kg/K], k = %g [W/m/K], A = %g [W/m3] \n", (LLD)(m->ID),m->alpha, m->Cp,m->k, m->A);
		PetscPrintf(PETSC_COMM_WORLD,"    \n");
	}
	else if (utype == _GEO_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: rho = %g [kg/m3], eta0 = %g [Pa.s]\n", (LLD)(m->ID), m->rho, eta);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (diff ) Bd = %g [1/(Pa.s)], Ed = %g [J/mol], Vd = %g [m3/mol] \n", (LLD)(m->ID), m->Bd, m->Ed, m->Vd);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (disl ) Bn = %g [1/(Pa^n.s)], En = %g [J/mol], Vn = %g [m3/mol], n = %g [-] \n", (LLD)(m->ID), m->Bn, m->En, m->Vn, m->n);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (peirl) Bp = %g [1/s], Ep = %g [J/mol], Vp = %g [m3/mol], taup = %g [MPa], gamma = %g [-], q = %g [-] \n", (LLD)(m->ID), m->Bp, m->Ep, m->Vp, m->taup, m->gamma, m->q);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (elast) G = %g [MPa], K = %g [MPa], Kp = %g [-] \n", (LLD)(m->ID), m->G, m->K, m->Kp);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (plast) cohesion = %g [MPa], friction angle = %g [deg] \n", (LLD)(m->ID),m->ch, m->fr);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (sweak) cohesion SoftLaw = %lld [-], friction SoftLaw = %lld [-] \n", (LLD)(m->ID),(LLD)m->chSoftID, (LLD)m->frSoftID);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (temp ) alpha = %g [1/K], cp = %g [J/kg/K], k = %g [W/m/K], A = %g [W/m3] \n", (LLD)(m->ID),m->alpha, m->Cp,m->k, m->A);
		PetscPrintf(PETSC_COMM_WORLD,"    \n");
	}
}
//---------------------------------------------------------------------------
// initialize SOFTENING LAWS from file
#undef __FUNCT__
#define __FUNCT__ "MatSoftInit"
PetscErrorCode MatSoftInit(JacRes *jr, UserCtx *usr)
{
	FILE        *fp;
	Soft_t       soft;
	PetscInt     i, ils, ile, count, count1;
	PetscInt    *ls,*le;
	PetscInt    *check_soft;
	PetscBool    empty_soft;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// print overview of softening laws from file
	PetscPrintf(PETSC_COMM_WORLD,"Softening laws read from %s: \n",usr->ParamFile);

	// allocate memory for arrays to store line info - with overhead
	ierr = PetscMalloc((size_t)(max_num_soft+_max_overhead_)*sizeof(PetscScalar), &ls); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)(max_num_soft+_max_overhead_)*sizeof(PetscScalar), &le); CHKERRQ(ierr);

	// open file, count and get the positions of the material structures in file - check of file was done before
	fp = fopen(usr->ParamFile, "r" );

	// read number of cells
	getLineStruct(fp, ls, le, &count, &count1,"<SofteningStart>","<SofteningEnd>");

	// error checking
	if(count != count1  ) SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incomplete softening-law structures! <SofteningStart>: %lld, <SofteningEnd>: %lld", (LLD)count, (LLD)count1);
	if(count > max_num_soft) SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "More softening laws found than max allowed! Actual: %lld, Max allowed: %lld", (LLD)count, (LLD)max_num_soft);

	// initialize phase checking array
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)count, &check_soft); CHKERRQ(ierr);
	ierr = PetscMemzero(check_soft, sizeof(PetscInt)*(size_t)count); CHKERRQ(ierr);
	for(i = 0; i < count; i++) check_soft[i] = 0;

	// read each individual phase
	for(i = 0; i < count; i++)
	{
		// set position in the file
		ils = ls[i];
		ile = le[i];

		// read from file or set default values zero
		ierr = MatSoftGetStruct( fp, &soft, ils, ile); CHKERRQ(ierr);

		// set structure
		jr->matSoft[soft.ID] = soft;

		// for error check
		check_soft[soft.ID] = 1;
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// check empty soft laws
	empty_soft = PETSC_FALSE;

	for(i = 0; i < count; i++)
	{
		if(!check_soft[i])
		{
			PetscPrintf(PETSC_COMM_WORLD, "Softening law %lld is not initialized\n", (LLD)i);
			empty_soft = PETSC_TRUE;
		}
	}

	if(empty_soft == PETSC_TRUE)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect softening laws numbering");
	}

	// store softening laws
	jr->numSoft = count;

	// free memory
	fclose(fp);

	// free arrays
	ierr = PetscFree(ls         ); CHKERRQ(ierr);
	ierr = PetscFree(le         ); CHKERRQ(ierr);
	ierr = PetscFree(check_soft);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Find softening laws in file
#undef __FUNCT__
#define __FUNCT__ "MatSoftGetStruct"
PetscErrorCode MatSoftGetStruct( FILE *fp, Soft_t *s, PetscInt ils, PetscInt ile)
{
	PetscInt    found;
	PetscFunctionBegin;

	// softening law id
	getMatPropInt   ( fp, ils, ile, "softID",               &s->ID,       &found);
	if (!found) SETERRQ (PETSC_COMM_WORLD, PETSC_ERR_USER, "MatSoftError: No softening law ID specified! ");

	// softening ratio
	getMatPropScalar( fp, ils, ile, "A",                    &s->A,        &found);

	// softening begin & end
	getMatPropScalar( fp, ils, ile, "APS1",                 &s->APS1,     &found);
	getMatPropScalar( fp, ils, ile, "APS2",                 &s->APS2,     &found);

	PetscPrintf(PETSC_COMM_WORLD,"    SoftLaw [%lld]: A = %g, APS1 = %g, APS2 = %g \n", (LLD)(s->ID), s->A, s->APS1, s->APS2);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// set material softening parameters in phases
#undef __FUNCT__
#define __FUNCT__ "MatSoftSet"
PetscErrorCode MatSoftSet(JacRes *jr)
{
	PetscInt     i, numPhases, numSoft;
	Soft_t      *matSoft, msoft;
	Material_t  *phases;

	PetscFunctionBegin;

	// access phase properties in user context variables
	numPhases       = jr->numPhases;
	numSoft         = jr->numSoft;
	phases          = jr->phases;
	matSoft         = jr->matSoft;

	// set zero values if there are no laws specified
	msoft.ID   = -1;
	msoft.A    = 0.0;
	msoft.APS1 = 0.0;
	msoft.APS2 = 0.0;

	// read phase parameters
	for(i = 0; i < numPhases; i++)
	{
		// cohesion softening law
		if(phases[i].chSoftID>-1)
		{
			if (phases[i].chSoftID > numSoft-1)
				{
				SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatSoftError: Softening law specified (%lld) is non-existent for phase ID: %lld! ",(LLD)phases[i].chSoftID, (LLD)i);
				}
			phases[i].chSoft = &matSoft[phases[i].chSoftID];
		}
		else
		{
			phases[i].chSoft = &msoft;
		}

		// friction softening law
		if(phases[i].frSoftID>-1) {
			if (phases[i].frSoftID > numSoft-1)
				{
				SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "MatSoftError: Softening law specified (%lld) is non-existent for phase ID: %lld! ",(LLD)phases[i].frSoftID, (LLD)i);
				}
			phases[i].frSoft = &matSoft[phases[i].frSoftID];
		}
		else
		{
			phases[i].frSoft = &msoft;
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// get the positions in the file for the material structures
void getLineStruct( FILE *fp, PetscInt *ls, PetscInt *le, PetscInt *count, PetscInt *count1, const char key[], const char key_end[])
{
	char line[MAX_LINE_LEN], next[MAX_LINE_LEN];
	PetscInt comment, start_match, end_match;
	PetscInt i = 0, ii = 0;

	// reset to start of file
	rewind(fp);

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );

		// get rid of white space
		trim(line);

		// if line is blank
		if( strlen(line) == 0 ) { continue; }

		// is first character a comment ?
		comment = is_comment_line( line );
		if( comment == _TRUE ) { continue; }

		// find start of material structure
		start_match = material_key_matches( key, line );
		if( start_match == _FALSE ) { continue; }

		// mark file position and scan until end_match is found
		ls[i++] = (PetscInt)ftell( fp );

		// find end of material structure
		end_match = material_key_matches( key_end, line );

		while( end_match == _FALSE ) {
			fgets( next, MAX_LINE_LEN-1, fp );

			// get rid of white space
			trim(next);

			// if line is blank
			if( strlen(next) == 0 ) { continue; }

			// is first character a comment ?
			comment = is_comment_line( next);
			if( comment == _TRUE ) { continue; }

			end_match = material_key_matches( key_end, next );
			if( end_match == _FALSE ) { continue; }

			// mark file position
			le[ii++] = (PetscInt)ftell( fp );
		}
	}

	*count  = i;
	*count1 = ii;
}
//---------------------------------------------------------------------------
// get integer within specified positions of the file
void getMatPropInt( FILE *fp, PetscInt ils, PetscInt ile, const char key[], PetscInt *value, PetscInt *found )
{
	char line[MAX_LINE_LEN];
	PetscInt comment, pos;
	PetscInt match, int_val;

	// init flag
	*found = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );
		pos = (PetscInt)ftell( fp );

		// search only within specified positions of the file
		if ((pos > ils) && (pos < ile))
		{
			// get rid of white space
			trim(line);

			// if line is blank
			if( strlen(line) == 0 ) { continue; }

			// is first character a comment ?
			comment = is_comment_line( line );
			if( comment == _TRUE ) {   continue;  }

			match = key_matches( key, line );
			if( match == _FALSE ) {   continue;   }

			// strip word and equal sign
			strip(line);

			int_val = (PetscInt)strtol( line, NULL, 0 );

			*value = int_val;
			*found = _TRUE;
			return;
		}
	}
}
//---------------------------------------------------------------------------
// get scalar within specified positions of the file
void getMatPropScalar( FILE *fp, PetscInt ils, PetscInt ile, const char key[], PetscScalar *value, PetscInt *found )
{
	char          line[MAX_LINE_LEN];
	PetscInt      comment, pos;
	PetscInt      match;
	PetscScalar   double_val;

	// init flag
	*found = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) ) {
		fgets( line, MAX_LINE_LEN-1, fp );
		pos = (PetscInt)ftell( fp );

		// search only within specified positions of the file
		if ((pos > ils) && (pos < ile))
		{
			// get rid of white space
			trim(line);

			// if line is blank
			if( strlen(line) == 0 ) { continue; }

			// is first character a comment ?
			comment = is_comment_line( line );
			if( comment == _TRUE ) {   continue;  }

			match = key_matches( key, line );
			if( match == _FALSE ) {   continue;   }

			// strip word and equal sign
			strip(line);

			double_val = (PetscScalar)strtod( line, NULL );

			*value = double_val;
			*found = _TRUE;
			return;
		}
	}
}
//---------------------------------------------------------------------------
