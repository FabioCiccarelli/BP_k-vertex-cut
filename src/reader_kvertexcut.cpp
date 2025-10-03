/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_kvertexcut.c
 * @brief  K-vertex cut problem reader for DIMACS format
 * @author Fabio Ciccarelli
 *
 * This file implements the reader/parser used to read k-vertex cut input data in DIMACS format.
 * The DIMACS format for graphs has the following structure:
 * - Comment lines (when present) start with 'c'
 * - Problem line: "p edge nnodes nedges"
 * - Edge lines: "e u v" where u and v are node indices (1-indexed)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_setppc.h"

#include "probdata_kvertexcut.h"
#include "reader_kvertexcut.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "kvcreader"
#define READER_DESC             "file reader for DIMACS graph data format"
#define READER_EXTENSION        "dimacs"

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadKvc)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   int* tail;
   int* head;
   SCIP_Bool error;

   char name[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   char keyword[SCIP_MAXSTRLEN];
   int nnodes;
   int nedges;
   int k;
   int nread;
   int u, v;
   int nedges_read;
   int lineno;

   *result = SCIP_DIDNOTRUN;

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   lineno = 0;
   nnodes = 0;
   nedges = 0;
   k = 1; /* default value, will be overridden by command line */

   /* extract problem name from filename */
   {
      const char* basename = strrchr(filename, '/');
      if( basename != NULL )
         basename++;
      else
         basename = filename;
      
      /* copy and remove extension */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", basename);
      char* dot = strrchr(name, '.');
      if( dot != NULL )
         *dot = '\0';
   }

   /* read header information */
   error = FALSE;
   while( !SCIPfeof(file) && !error )
   {
      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* skip comment lines */
      if( buffer[0] != 'e' && buffer[0] != 'p' )
         continue;

      /* parse problem line */
      if( buffer[0] == 'p' )
      {
         nread = sscanf(buffer, "%s %s %d %d", keyword, keyword, &nnodes, &nedges);
         if( nread != 4 || strcmp(keyword, "edge") != 0 )
         {
            SCIPwarningMessage(scip, "invalid problem line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            error = TRUE;
            break;
         }

         SCIPdebugMsg(scip, "problem: %d nodes, %d edges\n", nnodes, nedges);
         break;
      }
   }

   if( nnodes == 0 || nedges == 0 )
   {
      SCIPwarningMessage(scip, "no valid problem line found in file <%s>\n", filename);
      (void)SCIPfclose(file);
      return SCIP_READERROR;
   }

   /* allocate memory for edges */
   SCIP_CALL( SCIPallocBufferArray(scip, &tail, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &head, nedges) );

   /* read edges */
   nedges_read = 0;

   while( !SCIPfeof(file) && !error && nedges_read < nedges )
   {
      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* skip comment lines */
      if( buffer[0] != 'e' && buffer[0] != 'p' )
         continue;

      /* parse edge line */
      if( buffer[0] == 'e' )
      {
         nread = sscanf(buffer, "%s %d %d", keyword, &u, &v);
         if( nread != 3 )
         {
            SCIPwarningMessage(scip, "invalid edge line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            error = TRUE;
            break;
         }

         /* check if nodes are in valid range */
         if( u < 1 || u > nnodes || v < 1 || v > nnodes )
         {
            SCIPwarningMessage(scip, "invalid edge (%d,%d) in line %d: nodes must be in range [1,%d]\n", 
                              u, v, lineno, nnodes);
            error = TRUE;
            break;
         }

         tail[nedges_read] = u;
         head[nedges_read] = v;
         nedges_read++;

         SCIPdebugMsg(scip, "edge %d: (%d, %d)\n", nedges_read, u, v);
      }
   }

   if( nedges_read < nedges && !error )
   {
      SCIPwarningMessage(scip, "found only %d edges, expected %d in file <%s>\n", nedges_read, nedges, filename);
      error = TRUE;
   }

   if( !error )
   {
      /* TODO: get k value from command line parameters */
      /* For now, use a default value */
      SCIP_CALL( SCIPgetIntParam(scip, "k", &k) );;
      SCIPinfoMessage(scip, NULL, "Using k=%d\n", k);

      /* create k-vertex cut problem */
      SCIP_CALL( SCIPprobdataCreate(scip, name, nnodes, nedges, tail, head, k) );
   }

   /* cleanup */
   (void)SCIPfclose(file);
   SCIPfreeBufferArray(scip, &head);
   SCIPfreeBufferArray(scip, &tail);

   if( error )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the kvc file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderKvc(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create k-vertex cut reader data */
   readerdata = NULL;

   /* include k-vertex cut reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadKvc) );

   return SCIP_OKAY;
}

/**@} */
