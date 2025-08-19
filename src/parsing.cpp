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
PetscErrorCode FBLoad(FB **pfb, PetscBool DisplayOutput, char *restartFileName)
{
	FB        *fb;
	FILE      *fp;
	size_t    sz;
	PetscBool found;
	char      buffer[_str_len_], *filename, *all_options;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PetscMalloc(sizeof(FB), &fb); CHKERRQ(ierr);
	ierr = PetscMemzero(fb, sizeof(FB)); CHKERRQ(ierr);

	if(ISRankZero(PETSC_COMM_WORLD))
	{
		if(!restartFileName)
		{
			// check whether input file is specified
			ierr = PetscOptionsGetCheckString("-ParamFile", buffer, &found); CHKERRQ(ierr);

			if(found != PETSC_TRUE)
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Input file name is not specified. You must add the -ParamFile option to specify a LaMEM input file as in:  ./LaMEM -ParamFile your_input_file.dat \n");
			}

			filename = buffer;
		}
		else
		{
			// set restart input file
			filename = restartFileName;
		}

		// open input file
		fp = fopen(filename, "rb");

		// read additional PETSc options from input file
		if(fp == NULL)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open input file %s\n", filename);
		}

		if(DisplayOutput)
		{
			PetscPrintf(PETSC_COMM_WORLD, "Parsing input file : %s \n", filename);
		}

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

	// set simplified solver options from the input file
	ierr = solverOptionsReadFromFile(fb); CHKERRQ(ierr);

	// load additional options from file
	ierr = PetscOptionsReadFromFile(fb, DisplayOutput); CHKERRQ(ierr);

	// push command line options to the end of database (priority)
	ierr = PetscOptionsInsertString(NULL, all_options); CHKERRQ(ierr);
	
	// set required options (priority)
	ierr = solverOptionsSetRequired(); CHKERRQ(ierr);

	// print message
	if(DisplayOutput)
	{
		PetscPrintf(PETSC_COMM_WORLD, "Finished parsing input file \n");
	}

	// clean
	ierr = PetscFree(all_options); CHKERRQ(ierr);

	// return pointer
	(*pfb) = fb;

	if (DisplayOutput &&  ISRankZero(PETSC_COMM_WORLD)){
		PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
	}

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
	PetscBool found;
	char     *dbkey;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	found = PETSC_FALSE;

	// set defaults
	if(_default_) { ierr = PetscStrncpy(str, _default_, _str_len_); CHKERRQ(ierr); }
	else          { ierr = PetscMemzero(str,            _str_len_); CHKERRQ(ierr); }

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
// PETSc options parsing functions
//-----------------------------------------------------------------------------
PetscErrorCode PetscOptionsReadFromFile(FB *fb, PetscBool DisplayOutput)
{
	// * load additional options from input file
	// * push command line options to the end of database
	// (PETSc prioritizes options appearing LAST)

	PetscInt  jj, i, lnbeg, lnend;
	char     *line, **lines, *key, *val, *option;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
			key = strtok(line, " ");

			if(!key) continue;

			// get value
			val = strtok(NULL, " ");

			if(!val) option = key;
			else     asprintf(&option, "%s %s", key, val);

			// add to PETSc options
			if (DisplayOutput){
				PetscPrintf(PETSC_COMM_WORLD, "   Adding PETSc option: %s\n", option);
			}
			ierr = PetscOptionsInsertString(NULL, option); CHKERRQ(ierr);

			if(val) free(option);
		}

		fb->blockID++;
	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode PetscOptionsReadRestart(FILE *fp)
{
	// load options from restart file, replace existing

	size_t len;
	char   *all_options;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

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
PetscErrorCode PetscOptionsWriteRestart(FILE *fp)
{
	// save all existing options to restart file

	size_t len;
	char   *all_options;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	ierr = PetscOptionsGetAll(NULL, &all_options);  CHKERRQ(ierr);

	// include terminating null character
	len = strlen(all_options) + 1;

	fwrite(&len, sizeof(size_t), 1, fp);

	fwrite(all_options, sizeof(char)*len, 1, fp);

	ierr = PetscFree(all_options); CHKERRQ(ierr);

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
//---------------------------------------------------------------------------
PetscErrorCode set_integer_option(const char *key, const PetscInt val, const char *prefix = NULL)
{
	PetscFunctionBeginUser;
	char *opt;
	if(prefix) asprintf(&opt,"-%s_%s %lld", prefix, key, (LLD)val);
	else       asprintf(&opt,"-%s %lld",            key, (LLD)val);
	PetscCall(PetscOptionsInsertString(NULL, opt));
	free(opt);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode set_scalar_option(const char *key, const PetscScalar val, const char *prefix = NULL)
{
	PetscFunctionBeginUser;
	char *opt;
	if(prefix) asprintf(&opt,"-%s_%s %g", prefix, key, val);
	else       asprintf(&opt,"-%s %g",            key, val);
	PetscCall(PetscOptionsInsertString(NULL, opt));
	free(opt);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode set_string_option(const char *key, const char *val, const char *prefix = NULL)
{
	PetscFunctionBeginUser;
	char *opt;
	if(prefix) asprintf(&opt,"-%s_%s %s", prefix, key, val);
	else       asprintf(&opt,"-%s %s",            key, val);
	PetscCall(PetscOptionsInsertString(NULL, opt));
	free(opt);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode set_empty_option(const char *key, const char *prefix = NULL)
{
	PetscFunctionBeginUser;
	char *opt;
	if(prefix) asprintf(&opt,"-%s_%s", prefix, key);
	else       asprintf(&opt,"-%s",            key);
	PetscCall(PetscOptionsInsertString(NULL, opt));
	free(opt);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode set_mg_options(const char *prefix, PetscInt nlevels, PetscInt nsweeps, PetscScalar damping)
{
	PetscFunctionBeginUser;

	PetscCall(set_string_option ("pc_type", "mg", prefix));
	PetscCall(set_integer_option("pc_mg_levels", nlevels, prefix));
	PetscCall(set_empty_option  ("pc_mg_galerkin", prefix));
	PetscCall(set_string_option ("pc_mg_type", "multiplicative", prefix));
	PetscCall(set_string_option ("pc_mg_cycle_type", "v", prefix));
	PetscCall(set_string_option ("mg_levels_ksp_type", "richardson", prefix));
	PetscCall(set_scalar_option ("mg_levels_ksp_richardson_scale", damping, prefix));
	PetscCall(set_integer_option("mg_levels_ksp_max_it", nsweeps, prefix));
	PetscCall(set_string_option ("mg_levels_pc_type", "jacobi", prefix));

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
PetscErrorCode solverOptionsReadFromFile(FB *fb)
{
	// set "best-guess" solver options to help an inexperienced user
	// all options can be overridden by the usual PETSC options


	PetscMPIInt      size;

	
	PetscErrorCode ierr;
	PetscFunctionBeginUser;
	
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<SolverOptionsStart>", "<SolverOptionsEnd>"); CHKERRQ(ierr);

	// do not set options if not requested explicitly
	if(!fb->nblocks) PetscFunctionReturn(0);

	if(fb->nblocks > 1)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many solver options blocks. Only one is allowed");
	}

	// get number of ranks
	MPI_Comm_size(PETSC_COMM_WORLD, &size);


	// set defaults
	PetscScalar snes_param[ ]                  =  { 1e-5, -1.0, 50.0  }; // rtol, atol, maxit (-1 = automatic setting)
	PetscScalar ksp_param [ ]                  =  { 1e-6, -1.0, 200.0 }; // rtol, atol, maxit (-1 = automatic setting)
	PetscScalar picard_to_newton[ ]            =  { 1e-2,  1.2, 20.0  }; // picard-newton-rtol, newton-picard-rtol, picard-newton-maxit
	PetscInt    use_line_search                =  1;
	PetscInt    use_eisenstat_walker           =  0;
	PetscInt    use_mat_free_jac               =  0;
	PetscInt    mat_free_levels                =  0;
	char        stokes_solver[_str_len_];      // powell_hestenes        // [powell_hestenes, coupled_mg, block_mg, wbfbt];
	char        direct_solver_type[_str_len_]; // superlu_dist           // [mumps, superlu_dist, lu]
	PetscScalar penalty                        =  1e3;                   // only for powell_hestenes
	PetscInt    num_levels                     = -1;                     // -1 = automatic setting
	char        smoother_ksp[_str_len_];       // richardson             // [richardson, chebyshev, gmres];
	char        smoother_pc[_str_len_];        // jacobi                 // [jacobi, bjacobi, asm]
	PetscScalar smoother_damping               =  0.5;                   // only for richardson
	PetscInt    smoother_num_sweeps            =  10;                    // maxit
	PetscInt    coarse_num_cpu                 = -1;                     // -1 = automatic setting
	PetscInt    coarse_cells_per_cpu           =  2048;                  // required to set automatic value for coarse_num_cpu
	char        coarse_ksp[_str_len_];         // preonly                // [preonly, gmres];
	char        coarse_pc[_str_len_];          // direct                 // [direct, hypre, bjacobi, asm];
	PetscScalar coarse_param[]                 = { 1e-3, 100 } ;         // rtol, maxit
	PetscInt    subdomain_overlap              =  1  ;                   // only for asm
	PetscInt    subdomain_ilu_levels           =  0 ;                    // only for bjacobi and asm
	PetscInt    subdomain_num_cells            = -1 ;                    // -1 = all cells of cpu are assigned to one subdomain
	PetscScalar steady_thermal_param[ ]        = { 1e-8, 500.0 };        // tol, maxit
	char        steady_thermal_pc[_str_len_]; // mg                      // [mg, default]


	// read simplified solver options
	ierr = getScalarParam(fb, _OPTIONAL_, "snes_param",            snes_param,            3, 1.0);             CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "ksp_param",             ksp_param,             3, 1.0);             CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "picard_to_newton",      picard_to_newton,      3, 1.0);             CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "use_line_search",      &use_line_search,       1, 1);               CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "use_eisenstat_walker", &use_eisenstat_walker,  1, 1);               CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "use_mat_free_jac",     &use_mat_free_jac,      1, 1);               CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "mat_free_levels",      &mat_free_levels,       1, 16);              CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "stokes_solver",         stokes_solver,         "powell_hestenes" ); CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "direct_solver_type",    direct_solver_type,    "superlu_dist" );    CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "penalty",              &penalty,               1, 1.0);             CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "num_levels",           &num_levels,            1, 32);              CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "smoother_ksp",          smoother_ksp,          "richardson");       CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "smoother_pc",           smoother_pc,           "jacobi");           CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "smoother_damping",     &smoother_damping,      1, 1.0);             CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "smoother_num_sweeps",  &smoother_num_sweeps,   1, 1000);            CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "coarse_num_cpu",       &coarse_num_cpu,        1, 2096);            CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "coarse_cells_per_cpu", &coarse_cells_per_cpu,  1, 32768);           CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "coarse_ksp",            coarse_ksp,            "preonly");          CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "coarse_pc",             coarse_pc,             "direct");           CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "coarse_param",          coarse_param,          2, 1.0);             CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "subdomain_overlap",     &subdomain_overlap,    1, 10);              CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "subdomain_ilu_levels",  &subdomain_ilu_levels, 1, 8);               CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "subdomain_num_cells",   &subdomain_num_cells,  1, 65536);           CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "steady_thermal_param",  steady_thermal_param,  2, 1.0);             CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "steady_thermal_pc",     steady_thermal_pc,     "mg");               CHKERRQ(ierr);
/*

	if     (!strcmp(SolverType, "direct"))    solType = _DIRECT_STOKES_;
	else if(!strcmp(SolverType, "multigrid")) solType = _MULTIGRID_STOKES_;
	else if(!strcmp(SolverType, "block"))     solType = _BLOCK_STOKES_;
	else if(!strcmp(SolverType, "wbfbt"))     solType = _wBFBT_STOKES_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect solver type (SolverType): %s", SolverType);

	if(!(!strcmp(DirectSolver, "mumps") || !strcmp(DirectSolver, "superlu_dist")))
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect direct solver package (DirectSolver): %s", DirectSolver);
	}

	if     (!strcmp(MGCoarseSolver, "direct"))  crsType = _DIRECT_COARSE_;
	else if(!strcmp(MGCoarseSolver, "hypre"))   crsType = _HYPRE_COARSE_;
	else if(!strcmp(MGCoarseSolver, "asm"))     crsType = _ASM_COARSE_;
	else if(!strcmp(MGCoarseSolver, "bjacobi")) crsType = _BJACOBI_COARSE_;
	else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect coarse solver type (MGCoarseSolver): %s", MGCoarseSolver);

	// SNES
	ierr = PetscOptionsInsertString(NULL, "-snes_monitor");                        CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_rtol 1e-3");                      CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_atol 1e-7");                      CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_stol 1e-8");                      CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_max_funcs 1000000");              CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_max_it 50");                      CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_linesearch_type l2");             CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_linesearch_max_it 5");            CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_linesearch_maxstep 1.0");         CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_linesearch_minlambda 0.05");      CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_PicardSwitchToNewton_rtol 1e-2"); CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_NewtonSwitchToPicard_it 20");     CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-snes_NewtonSwitchToPicard_rtol 1.2");  CHKERRQ(ierr);

	// KSP
	ierr = PetscOptionsInsertString(NULL, "-js_ksp_type fgmres");      CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-js_ksp_rtol  1e-6");       CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-js_ksp_max_it 500");       CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-js_ksp_monitor");          CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(NULL, "-js_ksp_converged_reason"); CHKERRQ(ierr);

	if(solType == _DIRECT_STOKES_)
	{
		PetscCall(set_string_option("jp_type", "bf"));
		PetscCall(set_scalar_option("jp_pgamma", pgamma));
		PetscCall(set_string_option("bf_vs_type", "user"));
		PetscCall(set_string_option("vs_ksp_type", "preonly"));
		PetscCall(set_string_option("vs_pc_type","lu"));
		PetscCall(set_string_option("vs_pc_factor_mat_solver_type", DirectSolver));
	}
	else if(solType == _MULTIGRID_STOKES_)
	{
		PetscCall(set_string_option ("jp_type", "mg"));

		// muligrid defaults
		PetscCall(set_mg_options("gmg", nlevels, nsweeps, damping));
	}
	else if(solType == _BLOCK_STOKES_)
	{
		PetscCall(set_string_option("jp_type", "bf"));
		PetscCall(set_string_option("bf_vs_type", "mg"));
		PetscCall(set_string_option("vs_ksp_type", "preonly"));

		// muligrid defaults
		PetscCall(set_mg_options("gmg", nlevels, nsweeps, damping));
	}
	else if(solType == _wBFBT_STOKES_)
	{
		PetscCall(set_string_option("jp_type", "bf"));
		PetscCall(set_string_option("bf_vs_type", "mg"));
		PetscCall(set_empty_option ("bf_schur_wbfbt"));
		PetscCall(set_string_option("vs_ksp_type", "preonly"));
		PetscCall(set_string_option("ks_ksp_type", "preonly"));

		// muligrid defaults
		PetscCall(set_mg_options("gmg", nlevels, nsweeps, damping));
		PetscCall(set_mg_options("ks",  nlevels, nsweeps, damping));
	}
*/
	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

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
//---------------------------------------------------------------------------
