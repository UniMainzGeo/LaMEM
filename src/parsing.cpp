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
//...................   Input file parsing routines   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "tools.h"
//---------------------------------------------------------------------------
PetscErrorCode FBLoad(FB **pfb)
{
	FB        *fb;
	FILE      *fp;
	size_t    sz;
	PetscBool found;
	char      filename[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PetscMalloc(sizeof(FB), &fb); CHKERRQ(ierr);
	ierr = PetscMemzero(fb, sizeof(FB)); CHKERRQ(ierr);

	if(ISRankZero(PETSC_COMM_WORLD))
	{
		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

		// check whether input file is specified
		ierr = PetscOptionsGetCheckString("-ParamFile", filename, &found); CHKERRQ(ierr);

		if(found != PETSC_TRUE)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Input file name is not specified. You must add the -ParamFile option to specify a LaMEM input file as in:  ./LaMEM -ParamFile your_input_file.dat \n");
		}

		// open input file
		fp = fopen(filename, "rb");

		// read additional PETSc options from input file
		if(fp == NULL)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open input file %s\n", filename);
		}

		PetscPrintf(PETSC_COMM_WORLD, "Parsing input file : %s \n", filename);

		// get file size
		fseek(fp, 0L, SEEK_END);

		sz = (size_t)ftell(fp);

		rewind(fp);

		// read entire file into buffer
		ierr = PetscMalloc((sz + 1)*sizeof(char), &fb->fbuf); CHKERRQ(ierr);

		fread(fb->fbuf, sz*sizeof(char), 1, fp);

		fclose(fp);

		// pad with string terminator
		fb->fbuf[sz] = '\0';

		// set number of characters
		fb->nchar = (PetscInt)sz + 1;
	}

	// broadcast
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Bcast(&fb->nchar, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}

	if(!ISRankZero(PETSC_COMM_WORLD))
	{
		ierr = PetscMalloc((size_t)fb->nchar*sizeof(char), &fb->fbuf); CHKERRQ(ierr);
	}

	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Bcast(fb->fbuf, (PetscMPIInt)fb->nchar, MPI_CHAR, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	
	// parse buffer
	ierr = FBParseBuffer(fb); CHKERRQ(ierr);

	// print message
	PetscPrintf(PETSC_COMM_WORLD, "Finished parsing input file \n");

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	// return pointer
	(*pfb) = fb;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FBDestroy(FB **pfb)
{
	FB *fb;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// get pointer
	fb = (*pfb);

	if(!fb) PetscFunctionReturn(0);

	ierr = PetscFree(fb->fbuf);    CHKERRQ(ierr);
	ierr = PetscFree(fb->lbuf);    CHKERRQ(ierr);
	ierr = PetscFree(fb->pfLines); CHKERRQ(ierr);
	ierr = PetscFree(fb->pbLines); CHKERRQ(ierr);
	ierr = FBFreeBlocks(fb);       CHKERRQ(ierr);
	ierr = PetscFree(fb);          CHKERRQ(ierr);

	// clear pointer
	(*pfb) = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FBParseBuffer(FB *fb)
{
	char      *line, *b, p;
	size_t    len, maxlen;
	PetscInt  i, nchar, nlines, comment, cnt, block, *fblock;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// process buffer
	b     = fb->fbuf;
	nchar = fb->nchar;

	// purge line delimiters, replace tabs with spaces
	for(i = 0; i < nchar; i++)
	{
		if(b[i] == '\r') b[i] = '\0';
		if(b[i] == '\n') b[i] = '\0';
		if(b[i] == '\t') b[i] = ' ';
	}

	// purge comments
	for(i = 0, comment = 0; i < nchar; i++)
	{
		if(comment) { if(b[i] == '\0') { comment = 0; } b[i] = '\0';   }
		else        { if(b[i] == '#')  { comment = 1;   b[i] = '\0'; } }
	}

	// check equal signs
	for(i = 0; i < nchar; i++)
	{
		if(b[i] == '=')
		{
			if(!i)
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Input file cannot start with equal sign");
			}

			if(b[i-1] != ' ' || b[i+1] != ' ')
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Equal signs must be surrounded by spaces or tabs");
			}
		}
	}

	// purge empty lines, count actual number of lines
	for(i = 0, cnt = 0, p = '\0', nlines = 0; i < nchar; i++)
	{
		if(b[i] == '\0' && p == '\0') continue;
		p        = b[i];
		b[cnt++] = p;
		if(p == '\0') nlines++;
	}

	// collect garbage
	ierr = PetscMemzero(b + cnt, (size_t)(nchar - cnt)*sizeof(char)); CHKERRQ(ierr);

	// store actual number of characters
	fb->nchar = cnt;

	// count lines form flat and block access spaces, get line buffer size
	fb->nbLines = 0;
	fb->nfLines = 0;
	maxlen      = 0;

	ierr = makeIntArray(&fblock, NULL, nlines); CHKERRQ(ierr);

	for(i = 0, line = b, block = 0; i < nlines; i++)
	{
		if(block) { if(strstr(line, "<") && strstr(line, ">")) { block = 0; } fblock[i] = 1;   }
		else      { if(strstr(line, "<") && strstr(line, ">")) { block = 1;   fblock[i] = 1; } }

		if(fblock[i]) fb->nbLines++;
		else          fb->nfLines++;

		len = strlen(line);

		if(len > maxlen) maxlen = len;

		line += len + 1;
	}

	// allocate line buffer
	ierr = PetscMalloc((maxlen + 1)*sizeof(char), &fb->lbuf);         CHKERRQ(ierr);
	ierr = PetscMemzero(fb->lbuf, (size_t)(maxlen + 1)*sizeof(char)); CHKERRQ(ierr);

	// setup line pointers
	ierr = PetscMalloc((size_t)fb->nbLines*sizeof(char*), &fb->pbLines); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)fb->nfLines*sizeof(char*), &fb->pfLines); CHKERRQ(ierr);

	fb->nbLines = 0;
	fb->nfLines = 0;

	for(i = 0, line = b; i < nlines; i++)
	{
		if(fblock[i]) fb->pbLines[fb->nbLines++] = line;
		else          fb->pfLines[fb->nfLines++] = line;

		line += strlen(line) + 1;
	}

	ierr = PetscFree(fblock); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FBFindBlocks(FB *fb, ParamType ptype, const char *keybeg, const char *keyend)
{
	// find line ranges of data blocks

	PetscInt i, nbeg, nend;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	nbeg = 0;
	nend = 0;

	// count number of blocks
	for(i = 0; i < fb->nbLines; i++)
	{
		if(strstr(fb->pbLines[i], keybeg)) nbeg++;
		if(strstr(fb->pbLines[i], keyend)) nend++;
	}

	if(nbeg != nend)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "%s - %s identifiers don't match\n", keybeg, keyend);
	}

	fb->nblocks = nbeg;

	// check whether blocks are specified
	if(!fb->nblocks)
	{
		if     (ptype == _REQUIRED_) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "%s - %s blocks must be defined\n", keybeg, keyend);
		else if(ptype == _OPTIONAL_) PetscFunctionReturn(0);
	}

	// find & store block line ranges
	ierr = makeIntArray(&fb->blBeg, NULL, fb->nblocks); CHKERRQ(ierr);
	ierr = makeIntArray(&fb->blEnd, NULL, fb->nblocks); CHKERRQ(ierr);

	nbeg = 0;
	nend = 0;

	for(i = 0; i < fb->nbLines; i++)
	{
		if(strstr(fb->pbLines[i], keybeg)) fb->blBeg[nbeg++] = i+1;
		if(strstr(fb->pbLines[i], keyend)) fb->blEnd[nend++] = i;
	}

	// check block line ranges
	for(i = 0; i < fb->nblocks; i++)
	{
		if(fb->blBeg[i] >= fb->blEnd[i])
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect order of %s - %s identifiers\n", keybeg, keyend);
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FBFreeBlocks(FB *fb)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fb->nblocks = 0;
	fb->blockID = 0;

	ierr = PetscFree(fb->blBeg); CHKERRQ(ierr);
	ierr = PetscFree(fb->blEnd); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
char ** FBGetLineRanges(FB *fb, PetscInt *lnbeg, PetscInt *lnend)
{
	// return input file line ranges and pointers for parsing depending on access mode

	if(fb->nblocks)
	{
		// block access mode
		(*lnbeg) = fb->blBeg[fb->blockID];
		(*lnend) = fb->blEnd[fb->blockID];

		return fb->pbLines;
	}
	else
	{
		// flat access mode
		(*lnbeg) = 0;
		(*lnend) = fb->nfLines;

		return fb->pfLines;
	}
}
//-----------------------------------------------------------------------------
PetscErrorCode FBGetIntArray(
		FB         *fb,
		const char *key,
		PetscInt   *nvalues,
		PetscInt   *values,
		PetscInt    num,
		PetscBool  *found)
{
	PetscFunctionBeginUser;

	char     *ptr, *line, **lines;
	PetscInt  i, lnbeg, lnend, count;

	// initialize
	(*nvalues) = 0;
	(*found)   = PETSC_FALSE;

	// get line buffer & pointers
	line  = fb->lbuf;
	lines = FBGetLineRanges(fb, &lnbeg, &lnend);

	for(i = lnbeg; i < lnend; i++)
	{
		// copy line for parsing
		strcpy(line, lines[i]);

		// check for key match
		ptr = strtok(line, " ");

		if(!ptr || strcmp(ptr, key)) continue;

		// check equal sign
		ptr = strtok(NULL, " ");

		if(!ptr || strcmp(ptr, "="))
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);
		}

		// retrieve values after equal sign
		count = 0;
		ptr   = strtok(NULL, " ");

		while(ptr != NULL && count < num)
		{
			values[count++] = (PetscInt)strtol(ptr, NULL, 0);

			ptr = strtok(NULL, " ");
		}

		if(!count) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);

		(*nvalues) = count;
		(*found)   = PETSC_TRUE;

		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FBGetScalarArray(
		FB          *fb,
		const char  *key,
		PetscInt    *nvalues,
		PetscScalar *values,
		PetscInt     num,
		PetscBool   *found)
{
	PetscFunctionBeginUser;

	char     *ptr, *line, **lines;
	PetscInt  i, lnbeg, lnend, count;

	// initialize
	(*nvalues) = 0;
	(*found)   = PETSC_FALSE;

	// get line buffer & pointers
	line  = fb->lbuf;
	lines = FBGetLineRanges(fb, &lnbeg, &lnend);

	for(i = lnbeg; i < lnend; i++)
	{
		// copy line for parsing
		strcpy(line, lines[i]);

		// check for key match
		ptr = strtok(line, " ");

		if(!ptr || strcmp(ptr, key)) continue;

		// check equal sign
		ptr = strtok(NULL, " ");

		if(!ptr || strcmp(ptr, "="))
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);
		}

		// retrieve values after equal sign
		count = 0;
		ptr   = strtok(NULL, " ");

		while(ptr != NULL && count < num)
		{
			values[count++] = (PetscScalar)strtod(ptr, NULL);

			ptr = strtok(NULL, " ");
		}

		if(!count) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);

		(*nvalues) = count;
		(*found)   = PETSC_TRUE;

		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FBGetString(
		FB         *fb,
		const char *key,
		char       *str,    // output string
		PetscBool  *found)
{
	PetscFunctionBeginUser;

	char     *ptr, *line, **lines;
	PetscInt  i, lnbeg, lnend;

	// initialize
	(*found) = PETSC_FALSE;

	// get line buffer & pointers
	line  = fb->lbuf;
	lines = FBGetLineRanges(fb, &lnbeg, &lnend);

	for(i = lnbeg; i < lnend; i++)
	{
		// copy line for parsing
		strcpy(line, lines[i]);

		// check for key match
		ptr = strtok(line, " ");

		if(!ptr || strcmp(ptr, key)) continue;

		// check equal sign
		ptr = strtok(NULL, " ");

		if(!ptr || strcmp(ptr, "="))
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);
		}

		// retrieve values after equal sign
		ptr = strtok(NULL, " ");

		if(!ptr) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);

		// make sure string fits & is null terminated (two null characters are reserved in the end)
		if(strlen(ptr) > (_str_len_ - 2) )
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "String %s is more than %lld symbols long, (_str_len_ in parsing.h) \n", key, (LLD)(_str_len_ - 2));
		}

		// copy & pad the rest of the string with zeros
		strncpy(str, ptr, _str_len_);

		(*found) = PETSC_TRUE;

		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
// Wrappers
//-----------------------------------------------------------------------------
PetscErrorCode getIntParam(
		FB         *fb,
		ParamType   ptype,
		const char *key,
		PetscInt   *val,
		PetscInt    num,
		PetscInt    maxval)
{
	PetscInt  i, nval;
	PetscBool found;
	char     *dbkey;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(num < 1) PetscFunctionReturn(0);

	found = PETSC_FALSE;

	if(!fb->nblocks)
	{
		asprintf(&dbkey, "-%s", key);
	}
	else
	{
		asprintf(&dbkey, "-%s[%i]", key, (int) fb->ID);
	}

	nval = num;

	ierr = PetscOptionsGetIntArray(NULL, NULL, dbkey, val, &nval, &found); CHKERRQ(ierr);

	free(dbkey);

	if(found != PETSC_TRUE && fb)
	{
		ierr = FBGetIntArray(fb, key, &nval, val, num, &found); CHKERRQ(ierr);
	}

	// check whether parameter is set
	if(found != PETSC_TRUE)
	{
		if     (ptype == _REQUIRED_) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define parameter \"[-]%s\"\n", key);
		else if(ptype == _OPTIONAL_) PetscFunctionReturn(0);
	}

	// check number of entries
	if(nval < num) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "%lld entry(ies) are missing in parameter \"[-]%s\" \n",
		(LLD)(num-nval), key);

	// check for out-of-bound entries
	if(maxval > 0)
	{
		for(i = 0; i < num; i++)
		{
			if(val[i] > maxval)
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Entry %lld in parameter \"[-]%s\" is larger than allowed : val=%lld, max=%lld\n",
					(LLD)(i+1), key, (LLD)val[i], (LLD)maxval);
			}
		}
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode getScalarParam(
		FB          *fb,
		ParamType    ptype,
		const char  *key,
		PetscScalar *val,
		PetscInt     num,
		PetscScalar  scal)
{
	PetscInt  i, nval;
	PetscBool found;
	char     *dbkey;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	if(num < 1) PetscFunctionReturn(0);

	found = PETSC_FALSE;

	if(!fb->nblocks)
	{
		asprintf(&dbkey, "-%s", key);
	}
	else
	{
		asprintf(&dbkey, "-%s[%i]", key, (int) fb->ID);
	}
	
	nval = num;

	ierr = PetscOptionsGetScalarArray(NULL, NULL, dbkey, val, &nval, &found); CHKERRQ(ierr);
	
	free(dbkey);

	
	if(found != PETSC_TRUE && fb)
	{
		ierr = FBGetScalarArray(fb, key, &nval, val, num, &found); CHKERRQ(ierr);
	}

	// check data item exists
	if(found != PETSC_TRUE)
	{
		if     (ptype == _REQUIRED_) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define parameter \"[-]%s\"\n", key);
		else if(ptype == _OPTIONAL_) PetscFunctionReturn(0);
	}

	// check number of entries
	if(nval < num) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "%lld entry(ies) are missing in parameter \"[-]%s\" \n", (LLD)(num-nval), key);

	// nondimensionalize
	for(i = 0; i < num; i++) val[i] /= scal;

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode getStringParam(
		FB          *fb,
		ParamType    ptype,
		const char  *key,
		char        *str,        // output string
		const char  *_default_)  // default value (optional)
{
	// default ->  NULL             str -> cleared
	// default -> "_none_"          str -> not cleared
	// default ->  any other value  str -> default

	PetscBool found;
	char     *dbkey;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	found = PETSC_FALSE;

	// set defaults
	if(_default_) { if(strcmp(_default_, "_none_")) { PetscCall(PetscStrncpy(str, _default_, _str_len_)); } }
	else          {                                   PetscCall(PetscMemzero(str,            _str_len_)); }

	if(!fb->nblocks)
	{
		asprintf(&dbkey, "-%s", key);
	}
	else
	{
		asprintf(&dbkey, "-%s[%i]", key, (int) fb->ID);
	}
	
	ierr = PetscOptionsGetCheckString(dbkey, str, &found); CHKERRQ(ierr);

	free(dbkey);

	if(found != PETSC_TRUE && fb)
	{
		ierr = FBGetString(fb, key, str, &found);  CHKERRQ(ierr);
	}

	// check data item exists
	if(!strlen(str))
	{
		if     (ptype == _REQUIRED_) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define parameter \"[-]%s\"\n", key);
		else if(ptype == _OPTIONAL_) PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode  PetscOptionsGetCheckString(
	const char   key[],
	char         str[],
	PetscBool   *set)
{
	// prohibit empty parameters & check for overruns (two null characters are reserved in the end)

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL, NULL, key, str, _str_len_, set); CHKERRQ(ierr);

	if(*set && !strlen(str))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);
	}

	if(*set && strlen(str) > (_str_len_ - 2))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "String %s is more than %lld symbols long, (_str_len_ in parsing.h) \n", key, (LLD)(_str_len_ - 2));
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------

