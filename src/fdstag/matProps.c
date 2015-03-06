//---------------------------------------------------------------------------
//.................. MATERIAL PARAMETERS READING ROUTINES....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "solVar.h"
#include "fdstag.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "Parsing.h"
#include "matProps.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPropInit"
PetscErrorCode MatPropInit(JacRes *jr, UserCtx *usr)
{
	// initialize MATERIAL PARAMETERS from file

	FILE      *fp;
	PetscInt  *ls, *le;
	PetscInt   i, count_starts, count_ends;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// print overview of material parameters read from file
	PetscPrintf(PETSC_COMM_WORLD,"Phase material parameters read from %s: \n", usr->ParamFile);

	// clear memory
	ierr = PetscMemzero(jr->phases, sizeof(Material_t)*(size_t)max_num_phases); CHKERRQ(ierr);

	// initialize ID for consistency checks
	for(i = 0; i < max_num_phases; i++)
	{
		jr->phases[i].ID = -1;
	}

	// allocate memory for arrays to store line info
	ierr = makeIntArray(&ls, NULL, max_num_phases); CHKERRQ(ierr);
	ierr = makeIntArray(&le, NULL, max_num_phases); CHKERRQ(ierr);

	// open file, count and get the positions of the material structures in file - check of file was done before
	fp = fopen(usr->ParamFile, "r" );

	// read number of entries
	getLineStruct(fp, ls, le, max_num_phases, &count_starts, &count_ends, "<MaterialStart>","<MaterialEnd>");

	// error checking
	if(count_starts > max_num_phases || count_ends > max_num_phases)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many material structures specified! Max allowed: %lld", (LLD)max_num_phases);
	}
	if(count_starts != count_ends)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incomplete material structures! <MaterialStart> & <MaterialEnd> don't match");
	}

	// store actual number of materials
	jr->numPhases = count_starts;

	// read each individual phase
	for(i = 0; i < jr->numPhases; i++)
	{
		ierr = MatPropGetStruct(fp,
			jr->numPhases, jr->phases,
			jr->numSoft, jr->matSoft,
			ls[i], le[i], jr->scal.utype); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// free memory
	fclose(fp);

	// free arrays
	ierr = PetscFree(ls); CHKERRQ(ierr);
	ierr = PetscFree(le); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatPropGetStruct"
PetscErrorCode MatPropGetStruct(FILE *fp,
	PetscInt numPhases, Material_t *phases,
	PetscInt numSoft,   Soft_t     *matSoft,
	PetscInt ils, PetscInt ile, UnitsType utype)
{
	// read material properties from file with error checking
	// WARNING! This function assumes correctly defined softening parameters

	Material_t *m;
	PetscScalar eta, eta0, e0;
	PetscInt    ID, chSoftID, frSoftID, found;

	PetscFunctionBegin;

	// phase ID
	getMatPropInt(fp, ils, ile, "ID", &ID, &found);

	// error checking
	if(!found)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No phase ID specified! ");
	}
	if(ID > numPhases - 1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect phase numbering!");
	}

	// initialize additional parameters
	eta      =  0.0;
	eta0     =  0.0;
	e0       =  0.0;
	chSoftID = -1;
	frSoftID = -1;

	// get pointer to specified phase
	m = phases + ID;

	// check ID
	if(m->ID != -1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect phase numbering!");
	}

	// set ID
	m->ID = ID;

	//============================================================
	// density
	//============================================================
	getMatPropScalar(fp, ils, ile, "rho0",      &m->rho,   &found);
	//============================================================
	// Newtonian linear diffusion creep
	//============================================================
	getMatPropScalar(fp, ils, ile, "eta",       &eta,      &found);
	getMatPropScalar(fp, ils, ile, "Bd",        &m->Bd,    &found);
	getMatPropScalar(fp, ils, ile, "Ed",        &m->Ed,    &found);
	getMatPropScalar(fp, ils, ile, "Vd",        &m->Vd,    &found);
	//============================================================
	// power-law (dislocation) creep
	//============================================================
	getMatPropScalar(fp, ils, ile, "eta0",      &eta0,     &found);
	getMatPropScalar(fp, ils, ile, "e0",        &e0,       &found);
	getMatPropScalar(fp, ils, ile, "Bn",        &m->Bn,    &found);
	getMatPropScalar(fp, ils, ile, "n",         &m->n,     &found);
	getMatPropScalar(fp, ils, ile, "En",        &m->En,    &found);
	getMatPropScalar(fp, ils, ile, "Vn",        &m->Vn,    &found);
	//============================================================
	// Peierls creep
	//============================================================
	getMatPropScalar(fp, ils, ile, "Bp",        &m->Bp,    &found);
	getMatPropScalar(fp, ils, ile, "taup",      &m->taup,  &found);
	getMatPropScalar(fp, ils, ile, "gamma",     &m->gamma, &found);
	getMatPropScalar(fp, ils, ile, "q",         &m->q,     &found);
	getMatPropScalar(fp, ils, ile, "Ep",        &m->Ep,    &found);
	getMatPropScalar(fp, ils, ile, "Vp",        &m->Vp,    &found);
	//============================================================
	// elasticity
	//============================================================
	getMatPropScalar(fp, ils, ile, "shear",     &m->G,     &found);
	getMatPropScalar(fp, ils, ile, "bulk",      &m->K,     &found);
	getMatPropScalar(fp, ils, ile, "Kp",        &m->Kp,    &found);
	//============================================================
	// plasticity (Drucker-Prager)
	//============================================================
	getMatPropScalar(fp, ils, ile, "cohesion",  &m->ch,    &found);
	getMatPropScalar(fp, ils, ile, "friction",  &m->fr,    &found);
	getMatPropInt   (fp, ils, ile, "chSoftID",  &chSoftID, &found);
	getMatPropInt   (fp, ils, ile, "frSoftID",  &frSoftID, &found);
	//============================================================
	// energy
	//============================================================
	getMatPropScalar(fp, ils, ile, "alpha",     &m->alpha, &found);
	getMatPropScalar(fp, ils, ile, "cp",        &m->Cp,    &found);
	getMatPropScalar(fp, ils, ile, "k",         &m->k,     &found);
	getMatPropScalar(fp, ils, ile, "A",         &m->A,     &found);
	//============================================================

	// check softening laws
	if(chSoftID > numSoft-1)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect cohesion softening law specified for phase %lld", (LLD)ID);
	}
	if(frSoftID > numSoft-1)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect friction softening law specified for phase %lld", (LLD)ID);
	}

	if(m->fr && !m->ch)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Nonzero cohesion must be specified for phase %lld", (LLD)ID);
	}

	// set pointers to softening laws
	if(chSoftID != -1) m->chSoft = matSoft + chSoftID;
	if(frSoftID != -1) m->frSoft = matSoft + frSoftID;

	// check strain-rate dependent creep
	if((!eta0 && e0) || (eta0 && !e0))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "eta0 & e0 must be specified simultaneously for phase %lld", (LLD)ID);
	}

	// check power-law exponent
	if(!m->n && ((eta0 && e0) || m->Bn))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Power-law exponent must be specified for phase %lld", (LLD)ID);
	}

	// check Peierls creep
	if(m->Bp && (!m->taup || !m->gamma || !m->q || !m->Ep || !m->Vp))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "All Peierls creep parameters must be specified simultaneously for phase %lld", (LLD)ID);
	}

	// recompute creep parameters
	if(eta)        m->Bd = 1.0/(2.0*eta);
	if(eta0 && e0) m->Bn = pow (2.0*eta0, -m->n)*pow(e0, 1 - m->n);

	// check that at least one essential deformation mechanism is specified
	if(!m->Bd && !m->Bn && !m->G)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "At least one of the parameter (set) Bd (eta), Bn (eta0, e0), G must be specified for phase %lld", (LLD)ID);
	}

	// print
	// THIS PART IS A BIT UGLY, CONSIDER STRUCTURING

	if(utype == _NONE_)
	{
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: rho = %g [-], eta0 = %g [-]\n", (LLD)(m->ID), m->rho, eta);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (diff ) Bd = %g [-], Ed = %g [-], Vd = %g [-] \n", (LLD)(m->ID), m->Bd, m->Ed, m->Vd);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (disl ) Bn = %g [-], En = %g [-], Vn = %g [-], n = %g [-] \n", (LLD)(m->ID), m->Bn, m->En, m->Vn, m->n);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (peirl) Bp = %g [-], Ep = %g [-], Vp = %g [-], taup = %g [-], gamma = %g [-], q = %g [-] \n", (LLD)(m->ID), m->Bp, m->Ep, m->Vp, m->taup, m->gamma, m->q);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (elast) G = %g [-], K = %g [-], Kp = %g [-] \n", (LLD)(m->ID), m->G, m->K, m->Kp);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (plast) cohesion = %g [-], friction angle = %g [-] \n", (LLD)(m->ID),m->ch, m->fr);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (sweak) cohesion SoftLaw = %lld [-], friction SoftLaw = %lld [-] \n", (LLD)(m->ID),(LLD)chSoftID, (LLD)frSoftID);
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
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (sweak) cohesion SoftLaw = %lld [-], friction SoftLaw = %lld [-] \n", (LLD)(m->ID),(LLD)chSoftID, (LLD)frSoftID);
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
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (sweak) cohesion SoftLaw = %lld [-], friction SoftLaw = %lld [-] \n", (LLD)(m->ID),(LLD)chSoftID, (LLD)frSoftID);
		PetscPrintf(PETSC_COMM_WORLD,"    Phase [%lld]: (temp ) alpha = %g [1/K], cp = %g [J/kg/K], k = %g [W/m/K], A = %g [W/m3] \n", (LLD)(m->ID),m->alpha, m->Cp,m->k, m->A);
		PetscPrintf(PETSC_COMM_WORLD,"    \n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatSoftInit"
PetscErrorCode MatSoftInit(JacRes *jr, UserCtx *usr)
{
	// initialize SOFTENING LAWS from file

	FILE        *fp;
	PetscInt    *ls,*le;
	PetscInt     i, count_starts, count_ends;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// print overview of softening laws from file
	PetscPrintf(PETSC_COMM_WORLD,"Softening laws read from %s: \n",usr->ParamFile);

	// clear memory
	ierr = PetscMemzero(jr->matSoft, sizeof(Soft_t)*(size_t)max_num_soft); CHKERRQ(ierr);

	// initialize ID for consistency checks
	for(i = 0; i < max_num_soft; i++)
	{
		jr->matSoft[i].ID = -1;
	}

	// allocate memory for arrays to store line info
	ierr = makeIntArray(&ls, NULL, max_num_soft); CHKERRQ(ierr);
	ierr = makeIntArray(&le, NULL, max_num_soft); CHKERRQ(ierr);

	// open file, count and get the positions of the material structures in file - check of file was done before
	fp = fopen(usr->ParamFile, "r");

	// read number of entries
	getLineStruct(fp, ls, le, max_num_soft, &count_starts, &count_ends, "<SofteningStart>", "<SofteningEnd>");

	// error checking
	if(count_starts > max_num_soft || count_ends > max_num_soft)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many softening laws specified! Max allowed: %lld", (LLD)max_num_soft);
	}
	if(count_starts != count_ends)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incomplete material structures! <SofteningStart> & <SofteningEnd> don't match");
	}

	// store actual number of softening laws
	jr->numSoft = count_starts;

	// read each individual softening law
	for(i = 0; i < jr->numSoft; i++)
	{
		ierr = MatSoftGetStruct(fp, jr->numSoft, jr->matSoft, ls[i], le[i]); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// close file
	fclose(fp);

	// free arrays
	ierr = PetscFree(ls); CHKERRQ(ierr);
	ierr = PetscFree(le); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MatSoftGetStruct"
PetscErrorCode MatSoftGetStruct(FILE *fp,
	PetscInt numSoft, Soft_t *matSoft,
	PetscInt ils, PetscInt ile)
{
	// read softening law from file

	Soft_t   *s;
	PetscInt  found, ID;

	PetscFunctionBegin;

	// softening law ID
	getMatPropInt(fp, ils, ile, "softID", &ID, &found);

	// error checking
	if(!found)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No softening law ID specified! ");
	}
	if(ID > numSoft - 1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect softening law numbering!");
	}

	// get pointer to specified softening law
	s = matSoft + ID;

	// check ID
	if(s->ID != -1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect softening law numbering!");
	}

	// set ID
	s->ID = ID;

	// read and store softening law parameters
	getMatPropScalar(fp, ils, ile, "A",    &s->A,    &found);
	getMatPropScalar(fp, ils, ile, "APS1", &s->APS1, &found);
	getMatPropScalar(fp, ils, ile, "APS2", &s->APS2, &found);

	if(!s->A || !s->APS1 || !s->APS2)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "All parameters must be specified simultaneously for softening law %lld", (LLD)ID);
	}

	PetscPrintf(PETSC_COMM_WORLD,"SoftLaw [%lld]: A = %g, APS1 = %g, APS2 = %g \n", (LLD)(s->ID), s->A, s->APS1, s->APS2);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void getLineStruct(
	FILE *fp, PetscInt *ls, PetscInt *le, PetscInt max_num,
	PetscInt *count_starts, PetscInt *count_ends,
	const char key[], const char key_end[])
{
	// get the positions in the file for the material structures

	char line[MAX_LINE_LEN];
	PetscInt comment, start_match, end_match;
	PetscInt i = 0, ii = 0;

	// reset to start of file
	rewind(fp);

	while( !feof(fp) )
	{
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
		if( start_match == _TRUE )
		{
			// check bounds
			if(i+1 > max_num) { i++; return; }

			// mark file position and scan until end_match is found
			ls[i++] = (PetscInt)ftell( fp );
		}

		// find end of material structure
		end_match = material_key_matches( key_end, line );
		if( end_match == _TRUE )
		{
			// check bounds
			if(ii+1 > max_num) { ii++; return; }

			// mark file position
			le[ii++] = (PetscInt)ftell( fp );
		}

		// if no match continue
		if( (start_match == _FALSE) && (end_match == _FALSE) ) { continue; }
	}

	(*count_starts)  = i;
	(*count_ends)    = ii;
}
//---------------------------------------------------------------------------
void getMatPropInt(FILE *fp, PetscInt ils, PetscInt ile,
	const char key[], PetscInt *value, PetscInt *found)
{
	// get integer within specified positions of the file

	char line[MAX_LINE_LEN];
	PetscInt comment, pos;
	PetscInt match, int_val;

	// init flag
	(*found) = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) )
	{
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
void getMatPropScalar(FILE *fp, PetscInt ils, PetscInt ile,
	const char key[], PetscScalar *value, PetscInt *found)
{
	// get scalar within specified positions of the file

	char          line[MAX_LINE_LEN];
	PetscInt      comment, pos;
	PetscInt      match;
	PetscScalar   double_val;

	// init flag
	(*found) = _FALSE;

	// reset to start of file
	rewind( fp );

	while( !feof(fp) )
	{
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
