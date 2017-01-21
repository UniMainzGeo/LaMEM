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
	PetscInt  i;
	size_t    sz, len;
	PetscBool found;
	char      *ptr, *comment;
	char      filename[MAX_STR_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscMalloc(sizeof(FB), &fb); CHKERRQ(ierr);
	ierr = PetscMemzero(fb, sizeof(FB)); CHKERRQ(ierr);

	if(ISRankZero(PETSC_COMM_WORLD))
	{
		// check whether input file is specified
		ierr = PetscOptionsGetString(NULL, NULL, "-ParamFile", filename, MAX_STR_LEN, &found); CHKERRQ(ierr);

		// read additional PETSc options from input file
		if(found != PETSC_TRUE)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Input file name is not specified%s\n");
		}

		// open input file
		fp = fopen(filename, "r");

		if(fp == NULL)
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open input file %s\n", filename);
		}

		PetscPrintf(PETSC_COMM_WORLD," Parsing input file : %s \n", filename);

		// get file size
		fseek(fp, 0L, SEEK_END);

		sz = (size_t)ftell(fp);

		rewind(fp);

		// read entire file into buffer
		ierr = PetscMalloc((sz + 1)*sizeof(char), &fb->buff); CHKERRQ(ierr);

		fread(fb->buff, sz*sizeof(char), 1, fp);

		fclose(fp);

		// pad with line terminator
		fb->buff[sz] = '\0';

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
		ierr = PetscMalloc((size_t)fb->nchar*sizeof(char), &fb->buff); CHKERRQ(ierr);
	}

	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Bcast(fb->buff, (PetscMPIInt)fb->nchar, MPI_CHAR, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}

	// count lines, replace end-of-line characters with null terminators
	fb->nlines = 0;

	for(i = 0; i < fb->nchar; i++)
	{
		if(fb->buff[i] == '\n'
		|| fb->buff[i] == '\r') fb->buff[i] = '\0';
		if(fb->buff[i] == '\0') fb->nlines++;
	}

	// setup line pointers, trim comments
	ierr = PetscMalloc((size_t)fb->nlines*sizeof(char*), &fb->line); CHKERRQ(ierr);

	ptr = fb->buff;

	for(i = 0; i < fb->nlines; i++)
	{
		fb->line[i] = ptr;

		len = strlen(ptr);

		comment = strchr(ptr, '#');

		if(comment) strncpy(comment, "\0", strlen(comment));

		// update pointer
		ptr += len + 1;
	}

	// load additional options from file
	ierr = PetscOptionsReadFromFile(fb); CHKERRQ(ierr);

	// return pointer
	(*pfb) = fb;

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

	ierr = PetscFree(fb->buff); CHKERRQ(ierr);
	ierr = PetscFree(fb->line); CHKERRQ(ierr);
	ierr = FBFreeBlocks(fb);    CHKERRQ(ierr);
	ierr = PetscFree(fb);       CHKERRQ(ierr);

	// clear pointer
	(*pfb) = NULL;

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
	for(i = 0; i < fb->nlines; i++)
	{
		if(strstr(fb->line[i], keybeg)) nbeg++;
		if(strstr(fb->line[i], keyend)) nend++;
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

	for(i = 0; i < fb->nlines; i++)
	{
		if(strstr(fb->line[i], keybeg)) fb->blBeg[nbeg++] = i+1;
		if(strstr(fb->line[i], keyend)) fb->blEnd[nend++] = i;
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
#undef __FUNCT__
#define __FUNCT__ "FBGetRanges"
PetscErrorCode FBGetRanges(FB *fb, PetscInt *lnbeg, PetscInt *lnend)
{
	// return input file line ranges for parsing depending on access mode

	PetscFunctionBegin;

	if(fb->nblocks)
	{
		// block access mode
		(*lnbeg) = fb->blBeg[fb->blockID];
		(*lnend) = fb->blEnd[fb->blockID];
	}
	else
	{
		// flat access mode
		(*lnbeg) = 0;
		(*lnend) = fb->nlines;
	}

	PetscFunctionReturn(0);
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
	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscInt  val;
	char     *ptr, *pEnd;
	PetscInt  i, lnbeg, lnend, count;

	// initialize
	(*nvalues) = 0;
	(*found)   = PETSC_FALSE;

	ierr = FBGetRanges(fb, &lnbeg, &lnend); CHKERRQ(ierr);

	for(i = lnbeg; i < lnend; i++)
	{
		// check for key match
		if(!strstr(fb->line[i], key)) continue;

		// find equal sign
		ptr = strstr(fb->line[i], "=");

		if(!ptr) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);

		// retrieve values after equal sign
		ptr++;
		count = 0;

		while(count < num)
		{
			val = (PetscInt)strtol(ptr, &pEnd, 0);

			if(ptr == pEnd) break;

			ptr = pEnd;

			values[count++] = val;
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
	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscScalar  val;
	char        *ptr, *pEnd;
	PetscInt     i, lnbeg, lnend, count;

	// initialize
	(*nvalues) = 0;
	(*found)   = PETSC_FALSE;

	ierr = FBGetRanges(fb, &lnbeg, &lnend); CHKERRQ(ierr);

	for(i = lnbeg; i < lnend; i++)
	{
		// check for key match
		if(!strstr(fb->line[i], key)) continue;

		// find equal sign
		ptr = strstr(fb->line[i], "=");

		if(!ptr) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);

		// retrieve values after equal sign
		ptr++;
		count = 0;

		while(count < num)
		{
			val = (PetscScalar)strtod(ptr, &pEnd);

			if(ptr == pEnd) break;

			ptr = pEnd;

			values[count++] = val;
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
		size_t      fsz,    // full size of output string
		PetscBool  *found)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	char     *ptr, *tmp;
	PetscInt  i, lnbeg, lnend;

	// initialize
	(*found) = PETSC_FALSE;

	ierr = FBGetRanges(fb, &lnbeg, &lnend); CHKERRQ(ierr);

	for(i = lnbeg; i < lnend; i++)
	{
		// check for key match
		if(!strstr(fb->line[i], key)) continue;

		// find equal sign
		ptr = strstr(fb->line[i], "=");

		if(!ptr) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No equal sign specified for parameter \"%s\"\n", key);

		// retrieve first token after equal sign
		ptr++;

		// copy line, since strtok modifies it
		asprintf(&tmp, "%s", ptr);

		ptr = strtok(tmp, " \t");

		if(!ptr) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "No value specified for parameter \"%s\"\n", key);

		// make sure string fits & is null terminated
		if(strlen(ptr) + 1 > fsz)
		{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Input string is not large enough to hold result for parameter \"%s\" \n", key);
		}

		// copy & pad the rest of the string with zeros
		strncpy(str, ptr, fsz);

		free(tmp);

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
		if     (ptype == _REQUIRED_) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define parameter \"[-]%s\"\n", key);
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
		char        *str,  // output string
		size_t       fsz)  // full size of output string
{
	PetscBool found;
	char     *dbkey;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	found = PETSC_FALSE;

	if(!fb->nblocks)
	{
		asprintf(&dbkey, "-%s", key);

		ierr = PetscOptionsGetString(NULL, NULL, dbkey, str, fsz, &found); CHKERRQ(ierr);

		free(dbkey);
	}

	if(found != PETSC_TRUE && fb)
	{
		ierr = FBGetString(fb, key, str, fsz, &found);  CHKERRQ(ierr);
	}

	// check data item exists
	if(found != PETSC_TRUE)
	{
		if     (ptype == _REQUIRED_) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Define parameter \"[-]%s\"\n", key);
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
	// (PETSc prioritises options appearing LAST)

	char     *all_options, *tmp, *key, *val, *option;
	PetscInt  jj, i, lnbeg, lnend;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!fb) PetscFunctionReturn(0);

	// copy all command line options to buffer
	ierr = PetscOptionsGetAll(NULL, &all_options);  CHKERRQ(ierr);

	// remove command line options from database
	ierr = PetscOptionsClear(NULL); CHKERRQ(ierr);

	// setup block access mode
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<PetscOptionsStart>", "<PetscOptionsEnd>"); CHKERRQ(ierr);

	for(jj = 0; jj < fb->nblocks; jj++)
	{
		ierr = FBGetRanges(fb, &lnbeg, &lnend); CHKERRQ(ierr);

		for(i = lnbeg; i < lnend; i++)
		{
			// skip empty line
			if(!strlen(fb->line[i])) continue;

			// copy line, since strtok modifies it
			asprintf(&tmp, "%s", fb->line[i]);

			// get key
			key = strtok(tmp,  " \t");

			if(!key) { free(tmp); continue; }

			// get value
			val = strtok(NULL, " \t");

			if(!val) option = key;
			else     asprintf(&option, "%s %s", key, val);

			// add to PETSc options
			PetscPrintf(PETSC_COMM_WORLD, "Adding PETSc option: %s\n", option);

			ierr = PetscOptionsInsertString(NULL, option); CHKERRQ(ierr);

			if(val) free(option);
			free(tmp);
		}

		fb->blockID++;
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	// push command line options to the end of database (priority)
	ierr = PetscOptionsInsertString(NULL, all_options); CHKERRQ(ierr);

	ierr = PetscFree(all_options); CHKERRQ(ierr);

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
