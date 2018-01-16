/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   parsing.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//...................   Input file parsing routines   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "tools.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FBLoad"
PetscErrorCode FBLoad(FB **pfb)
{
	FB        *fb;
	FILE      *fp;
	size_t    sz;
	PetscBool found;
	char      filename[_STR_LEN_], *all_options;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscMalloc(sizeof(FB), &fb); CHKERRQ(ierr);
	ierr = PetscMemzero(fb, sizeof(FB)); CHKERRQ(ierr);

	if(ISRankZero(PETSC_COMM_WORLD))
	{
		// check whether input file is specified
		ierr = PetscOptionsGetCheckString("-ParamFile", filename, &found); CHKERRQ(ierr);

		// read additional PETSc options from input file
		if(found != PETSC_TRUE)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Input file name is not specified\n");
		}

		// open input file
		fp = fopen(filename, "r");

		if(fp == NULL)
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open input file %s\n", filename);
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

	// copy all command line and previously specified options to buffer
	ierr = PetscOptionsGetAll(NULL, &all_options);  CHKERRQ(ierr);

	// remove command line options from database
	ierr = PetscOptionsClear(NULL); CHKERRQ(ierr);

	// Set default solver options if defined in file
	ierr = StokesSetDefaultSolverOptions(fb); CHKERRQ(ierr);

	// load additional options from file
	ierr = PetscOptionsReadFromFile(fb); CHKERRQ(ierr);

	// push command line options to the end of database (priority)
	ierr = PetscOptionsInsertString(NULL, all_options); CHKERRQ(ierr);
	
	// clean
	ierr = PetscFree(all_options); CHKERRQ(ierr);

	// return pointer
	(*pfb) = fb;

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FBDestroy"
PetscErrorCode FBDestroy(FB **pfb)
{
	FB *fb;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
#undef __FUNCT__
#define __FUNCT__ "FBParseBuffer"
PetscErrorCode FBParseBuffer(FB *fb)
{
	char      *line, *b, p;
	size_t    len, maxlen;
	PetscInt  i, nchar, nlines, comment, cnt, block, *fblock;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// process buffer
	b     = fb->fbuf;
	nchar = fb->nchar;

	// purge line delimiters
	for(i = 0; i < nchar; i++)
	{
		if(b[i] == '\r') b[i] = '\0';
		if(b[i] == '\n') b[i] = '\0';
	}

	// purge comments
	for(i = 0, comment = 0; i < nchar; i++)
	{
		if(comment) { if(b[i] == '\0') { comment = 0; } b[i] = '\0';   }
		else        { if(b[i] == '#')  { comment = 1;   b[i] = '\0'; } }
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
#undef __FUNCT__
#define __FUNCT__ "FBFindBlocks"
PetscErrorCode FBFindBlocks(FB *fb, ParamType ptype, const char *keybeg, const char *keyend)
{
	// find line ranges of data blocks

	PetscInt i, nbeg, nend;

	PetscErrorCode ierr;
	PetscFunctionBegin;

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
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "%s - %s identifiers don't match\n", keybeg, keyend);
	}

	fb->nblocks = nbeg;

	// check whether blocks are specified
	if(!fb->nblocks)
	{
		if     (ptype == _REQUIRED_) SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "%s - %s blocks must be defined\n", keybeg, keyend);
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
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect order of %s - %s identifiers\n", keybeg, keyend);
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FBFreeBlocks"
PetscErrorCode FBFreeBlocks(FB *fb)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

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
#undef __FUNCT__
#define __FUNCT__ "FBGetIntArray"
PetscErrorCode FBGetIntArray(
		FB         *fb,
		const char *key,
		PetscInt   *nvalues,
		PetscInt   *values,
		PetscInt    num,
		PetscBool  *found)
{
	PetscFunctionBegin;

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
		ptr = strtok(line, " \t");

		if(!ptr || strcmp(ptr, key)) continue;

		// check equal sign
		ptr = strtok(NULL, " \t");

		if(!ptr || strcmp(ptr, "="))
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);
		}

		// retrieve values after equal sign
		count = 0;
		ptr   = strtok(NULL, " \t");

		while(ptr != NULL && count < num)
		{
			values[count++] = (PetscInt)strtol(ptr, NULL, 0);

			ptr = strtok(NULL, " \t");
		}

		if(!count) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);

		(*nvalues) = count;
		(*found)   = PETSC_TRUE;

		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FBGetScalarArray"
PetscErrorCode FBGetScalarArray(
		FB          *fb,
		const char  *key,
		PetscInt    *nvalues,
		PetscScalar *values,
		PetscInt     num,
		PetscBool   *found)
{
	PetscFunctionBegin;

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
		ptr = strtok(line, " \t");

		if(!ptr || strcmp(ptr, key)) continue;

		// check equal sign
		ptr = strtok(NULL, " \t");

		if(!ptr || strcmp(ptr, "="))
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);
		}

		// retrieve values after equal sign
		count = 0;
		ptr   = strtok(NULL, " \t");

		while(ptr != NULL && count < num)
		{
			values[count++] = (PetscScalar)strtod(ptr, NULL);

			ptr = strtok(NULL, " \t");
		}

		if(!count) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);

		(*nvalues) = count;
		(*found)   = PETSC_TRUE;

		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FBGetString"
PetscErrorCode FBGetString(
		FB         *fb,
		const char *key,
		char       *str,    // output string
		PetscBool  *found)
{
	PetscFunctionBegin;

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
		ptr = strtok(line, " \t");

		if(!ptr || strcmp(ptr, key)) continue;

		// check equal sign
		ptr = strtok(NULL, " \t");

		if(!ptr || strcmp(ptr, "="))
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);
		}

		// retrieve values after equal sign
		ptr = strtok(NULL, " \t");

		if(!ptr) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);

		// make sure string fits & is null terminated (two null characters are reserved in the end)
		if(strlen(ptr) > _STR_LEN_-2)
		{
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "String %s is more than %lld symbols long, (_STR_LEN_ in parsing.h) \"%s\" \n", key, _STR_LEN_-2);
		}

		// copy & pad the rest of the string with zeros
		strncpy(str, ptr, _STR_LEN_);

		(*found) = PETSC_TRUE;

		PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
// Wrappers
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getIntParam"
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
	PetscFunctionBegin;

	if(num < 1) PetscFunctionReturn(0);

	found = PETSC_FALSE;

	// PETSc options are not checked in block access mode
	if(!fb->nblocks)
	{
		asprintf(&dbkey, "-%s", key);

		nval = num;

		ierr = PetscOptionsGetIntArray(NULL, NULL, dbkey, val, &nval, &found); CHKERRQ(ierr);

		free(dbkey);
	}

	if(found != PETSC_TRUE && fb)
	{
		ierr = FBGetIntArray(fb, key, &nval, val, num, &found); CHKERRQ(ierr);
	}

	// check whether parameter is set
	if(found != PETSC_TRUE)
	{
		if     (ptype == _REQUIRED_) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define parameter \"[-]%s\"\n", key);
		else if(ptype == _OPTIONAL_) PetscFunctionReturn(0);
	}

	// check number of entries
	if(nval < num) SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "%lld entry(ies) are missing in parameter \"[-]%s\" \n",
		(LLD)(num-nval), key);

	// check for out-of-bound entries
	if(maxval > 0)
	{
		for(i = 0; i < num; i++)
		{
			if(val[i] > maxval)
			{
				SETERRQ4(PETSC_COMM_WORLD, PETSC_ERR_USER, "Entry %lld in parameter \"[-]%s\" is larger than allowed : val=%lld, max=%lld\n",
					(LLD)(i+1), key, (LLD)val[i], (LLD)maxval);
			}
		}
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getScalarParam"
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
	PetscFunctionBegin;

	if(num < 1) PetscFunctionReturn(0);

	found = PETSC_FALSE;

	if(!fb->nblocks)
	{
		asprintf(&dbkey, "-%s", key);

		nval = num;

		ierr = PetscOptionsGetScalarArray(NULL, NULL, dbkey, val, &nval, &found); CHKERRQ(ierr);

		free(dbkey);
	}

	if(found != PETSC_TRUE && fb)
	{
		ierr = FBGetScalarArray(fb, key, &nval, val, num, &found); CHKERRQ(ierr);
	}

	// check data item exists
	if(found != PETSC_TRUE)
	{
		if     (ptype == _REQUIRED_) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define parameter \"[-]%s\"\n", key);
		else if(ptype == _OPTIONAL_) PetscFunctionReturn(0);
	}

	// check number of entries
	if(nval < num) SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "%lld entry(ies) are missing in parameter \"[-]%s\" \n", (LLD)(num-nval), key);

	// nondimensionalize
	for(i = 0; i < num; i++) val[i] /= scal;

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "getStringParam"
PetscErrorCode getStringParam(
		FB          *fb,
		ParamType    ptype,
		const char  *key,
		char        *str,        // output string
		const char  *_default_)  // default value (optional)
{
	PetscBool found;
	char     *dbkey;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	found = PETSC_FALSE;

	// set defaults
	if(_default_) { ierr = PetscStrncpy(str, _default_, _STR_LEN_); CHKERRQ(ierr); }
	else          { ierr = PetscMemzero(str,            _STR_LEN_); CHKERRQ(ierr); }

	if(!fb->nblocks)
	{
		asprintf(&dbkey, "-%s", key);

		ierr = PetscOptionsGetCheckString(dbkey, str, &found); CHKERRQ(ierr);

		free(dbkey);
	}

	if(found != PETSC_TRUE && fb)
	{
		ierr = FBGetString(fb, key, str, &found);  CHKERRQ(ierr);
	}

	// check data item exists
	if(!strlen(str))
	{
		if     (ptype == _REQUIRED_) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define parameter \"[-]%s\"\n", key);
		else if(ptype == _OPTIONAL_) PetscFunctionReturn(0);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
// PETSc options parsing functions
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PetscOptionsReadFromFile"
PetscErrorCode PetscOptionsReadFromFile(FB *fb)
{
	// * load additional options from input file
	// * push command line options to the end of database
	// (PETSc prioritizes options appearing LAST)

	PetscBool found;
	PetscInt  jj, i, lnbeg, lnend;
	char     *line, **lines, *key, *val, *option, filename[_STR_LEN_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!fb) PetscFunctionReturn(0);


	// setup block access mode
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<PetscOptionsStart>", "<PetscOptionsEnd>"); CHKERRQ(ierr);

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
			key = strtok(line, " \t");

			if(!key) continue;

			// get value
			val = strtok(NULL, " \t");

			if(!val) option = key;
			else     asprintf(&option, "%s %s", key, val);

			// add to PETSc options
			PetscPrintf(PETSC_COMM_WORLD, "   Adding PETSc option: %s\n", option);

			ierr = PetscOptionsInsertString(NULL, option); CHKERRQ(ierr);

			if(val) free(option);
		}

		fb->blockID++;
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	// print message
	ierr = PetscOptionsGetCheckString("-ParamFile", filename, &found); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Finished parsing input file : %s \n", filename);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PetscOptionsReadRestart"
PetscErrorCode PetscOptionsReadRestart(FILE *fp)
{
	// load options from restart file, replace existing

	size_t len;
	char   *all_options;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsClear(NULL); CHKERRQ(ierr);

	// length already includes terminating null character
	fread(&len, sizeof(size_t), 1, fp);

	ierr = PetscMalloc(sizeof(char)*len, &all_options); CHKERRQ(ierr);

	fread(all_options, sizeof(char)*len, 1, fp); CHKERRQ(ierr);

	ierr = PetscOptionsInsertString(NULL, all_options); CHKERRQ(ierr);

	ierr = PetscFree(all_options); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PetscOptionsWriteRestart"
PetscErrorCode PetscOptionsWriteRestart(FILE *fp)
{
	// save all existing options to restart file

	size_t len;
	char   *all_options;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsGetAll(NULL, &all_options);  CHKERRQ(ierr);

	// include terminating null character
	len = strlen(all_options) + 1;

	fwrite(&len, sizeof(size_t), 1, fp);

	fwrite(all_options, sizeof(char)*len, 1, fp);

	ierr = PetscFree(all_options); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PetscOptionsGetCheckString"
PetscErrorCode  PetscOptionsGetCheckString(
	const char   key[],
	char         str[],
	PetscBool   *set)
{
	// prohibit empty parameters & check for overruns (two null characters are reserved in the end)

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsGetString(NULL, NULL, key, str, _STR_LEN_, set); CHKERRQ(ierr);

	if((*set) == PETSC_TRUE && !strlen(str))
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);
	}

	if(strlen(str) > _STR_LEN_-2)
	{
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "String %s is more than %lld symbols long, (_STR_LEN_ in parsing.h) \"%s\" \n", key, _STR_LEN_-2);
	}

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "StokesSetDefaultSolverOptions"
PetscErrorCode StokesSetDefaultSolverOptions(FB *fb)
{
	PetscErrorCode ierr;
 	char     		SolverType[_STR_LEN_], DirectSolver[_STR_LEN_], str[_STR_LEN_], SmootherType[_STR_LEN_];
	PetscScalar 	scalar;
	PetscInt 		integer;
	
	PetscFunctionBegin;
	
	// Set some 'best-guess' default solver paramaters to help the average user
	// All options can be overridden by the usual PETSC options
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting default solver options ****** \n"); 

	// Set default parameters for the outer iterations
	ierr = PetscOptionsInsertString(NULL, "-js_ksp_monitor"); 			CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-js_ksp_converged_reason"); 	CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-js_ksp_min_it 1"); 			CHKERRQ(ierr);
	
	// Read input file to see if we set solver options 
	ierr = getStringParam(fb, _OPTIONAL_, "SolverType",          SolverType,         NULL);          CHKERRQ(ierr);

	// depending on what was chosen, set the command-line parameters accordingly
	// These parameters can be overruled by parameters given in the PetscOptionsStart/PetscOptionsEnd block or on the command-line
	if     		(!strcmp(SolverType, "direct")){
		// Direct solver
		ierr = PetscOptionsInsertString(NULL, "-pcmat_type mono"); 	CHKERRQ(ierr);
		ierr = PetscOptionsInsertString(NULL, "-jp_type user"); 	CHKERRQ(ierr);
		ierr = PetscOptionsInsertString(NULL, "-jp_pc_type lu"); 	CHKERRQ(ierr);
		
		// Set penalty parameter if specified
		scalar 	= NULL;
		ierr 	= getScalarParam(fb, _OPTIONAL_, "DirectPenalty",       &scalar,        1, 1.0);          CHKERRQ(ierr);
		if (scalar){
 			sprintf(str, "-pcmat_pgamma %e", scalar);	ierr = PetscOptionsInsertString(NULL, str); 	CHKERRQ(ierr);
		}

		// if the type of direct solver is specified, use that
		ierr = getStringParam(fb, _OPTIONAL_, "DirectSolver",        DirectSolver,       NULL );          CHKERRQ(ierr);
		if     		(!strcmp(DirectSolver, "mumps")){        ierr = PetscOptionsInsertString(NULL, "-jp_pc_factor_mat_solver_package mumps"); 			CHKERRQ(ierr); }
		else if     (!strcmp(DirectSolver, "superlu_dist")){ ierr = PetscOptionsInsertString(NULL, "-jp_pc_factor_mat_solver_package superlu_dist"); 	CHKERRQ(ierr); }
		else if     (!strcmp(DirectSolver, "pastix"))	   { ierr = PetscOptionsInsertString(NULL, "-jp_pc_factor_mat_solver_package pastix"); 			CHKERRQ(ierr); }
		else {
			// solver is not specified; set a default one
			if (ISParallel(PETSC_COMM_WORLD))
			{
				// We need to set one of the parallel solvers. Determine if we have one of them installed in the current PETSC version
				ierr = PetscOptionsInsertString(NULL, "-jp_pc_factor_mat_solver_package superlu_dist"); 	CHKERRQ(ierr);
			}
		}

	
		// Report
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Solver Type 		: 	%s \n",SolverType); 

	}
	else if 	(!strcmp(SolverType, "multigrid")){
		// Multigrid solver

		ierr = PetscOptionsInsertString(NULL, "-pcmat_type mono"); 					CHKERRQ(ierr);
		ierr = PetscOptionsInsertString(NULL, "-jp_type mg"); 						CHKERRQ(ierr);
		ierr = PetscOptionsInsertString(NULL, "-gmg_pc_type mg"); 					CHKERRQ(ierr);
		ierr = PetscOptionsInsertString(NULL, "-gmg_pc_mg_galerkin"); 				CHKERRQ(ierr);
		ierr = PetscOptionsInsertString(NULL, "-gmg_pc_mg_type multiplicative"); 	CHKERRQ(ierr);
		ierr = PetscOptionsInsertString(NULL, "-gmg_pc_mg_cycle_type v"); 			CHKERRQ(ierr);
		ierr = PetscOptionsInsertString(NULL, "-gmg_pc_mg_log"); 					CHKERRQ(ierr);

		integer 	= 3;
		ierr 	= getIntParam(fb, _OPTIONAL_, "MGLevels",       &integer,        1, 100);          CHKERRQ(ierr);
		if (integer){
 			sprintf(str, "-gmg_pc_mg_levels %i", integer);	ierr = PetscOptionsInsertString(NULL, str); 	CHKERRQ(ierr);
		}
		
		integer 	= 10;
		ierr 	= getIntParam(fb, _OPTIONAL_, "MGSweeps",       &integer,        1, 100);          CHKERRQ(ierr);
		if (integer){
 			sprintf(str, "-gmg_mg_levels_ksp_max_it %i", integer);	ierr = PetscOptionsInsertString(NULL, str); 	CHKERRQ(ierr);
		}

		/* Specify smoother type options */
		ierr = getStringParam(fb, _OPTIONAL_, "MGSmoother",          SmootherType,         "chebyshev");          CHKERRQ(ierr);
		if 	(!strcmp(SmootherType, "jacobi")){

			ierr = PetscOptionsInsertString(NULL, "-gmg_mg_levels_ksp_type richardson"); 		CHKERRQ(ierr);
			ierr = PetscOptionsInsertString(NULL, "-gmg_mg_levels_pc_type jacobi"); 			CHKERRQ(ierr);

			scalar 	= 0.6;
			ierr 	= getScalarParam(fb, _OPTIONAL_, "MGJacobiDamp",       &scalar,        1, 1.0);          CHKERRQ(ierr);
			if (scalar){
 				sprintf(str, "-gmg_mg_levels_ksp_richardson_scale %f", scalar);	ierr = PetscOptionsInsertString(NULL, str); 	CHKERRQ(ierr);
			}
		}
		else if (!strcmp(SmootherType, "chebyshev")){
			ierr = PetscOptionsInsertString(NULL, "-gmg_mg_levels_ksp_type chebyshev"); 		CHKERRQ(ierr);

		}

		/* Specify coarse grid direct solver options */
		ierr = getStringParam(fb, _OPTIONAL_, "MGCoarseSolver",          SolverType,         "direct");          CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"MGCoarseSolver=%s \n",SolverType);
		if 	(!strcmp(SolverType, "direct") | !strcmp(SolverType, "mumps") | !strcmp(SolverType, "superlu_dist")){
			PetscPrintf(PETSC_COMM_WORLD,"Setting direct solver options for crs \n",SolverType);

			ierr = PetscOptionsInsertString(NULL, "-crs_ksp_type preonly"); 		CHKERRQ(ierr);
			ierr = PetscOptionsInsertString(NULL, "-crs_pc_type lu"); 		CHKERRQ(ierr);
			if (ISParallel(PETSC_COMM_WORLD)){
				if (!strcmp(SolverType, "direct") | !strcmp(SolverType, "superlu_dist")){
					ierr = PetscOptionsInsertString(NULL, "-crs_pc_factor_mat_solver_package superlu_dist"); 		CHKERRQ(ierr);
				}
				else if (!strcmp(SolverType, "mumps")){
					ierr = PetscOptionsInsertString(NULL, "-crs_pc_factor_mat_solver_package mumps"); 		CHKERRQ(ierr);
				}
			}
		}
		else if (!strcmp(SolverType, "redundant")){
			ierr = PetscOptionsInsertString(NULL, "-crs_ksp_type preonly"); 		CHKERRQ(ierr);
			ierr = PetscOptionsInsertString(NULL, "-crs_pc_type redundant"); 		CHKERRQ(ierr);
			
			// define number of redundant solves
			integer 	= 	4;
			ierr 		= 	getIntParam(fb, _OPTIONAL_, "MGRedundantNum",       &integer,        1, 100);          CHKERRQ(ierr);
			sprintf(str, "-crs_pc_redundant_number %i", integer);	ierr = PetscOptionsInsertString(NULL, str); 	CHKERRQ(ierr);

			ierr = getStringParam(fb, _OPTIONAL_, "MGRedundantSolver",          SolverType,         "superlu_dist");          CHKERRQ(ierr);
			sprintf(str, "-crs_redundant_pc_factor_mat_solver_package %s", SolverType);	ierr = PetscOptionsInsertString(NULL, str); 	CHKERRQ(ierr);
			

		}
		// More can be added here later, such as telescope etc. (once we have a bit more experience with those solvers)

	} 
	else{
		ierr = PetscPrintf(PETSC_COMM_WORLD,"NOT SETTING DEFAULT OPTIONS ****** \n"); 
	} 


	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
