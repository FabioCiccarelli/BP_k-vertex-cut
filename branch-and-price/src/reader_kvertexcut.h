/**@file   reader_kvertexcut.h
 * @brief  k-vertex cut problem reader file
 * @author Fabio Ciccarelli
 *
 * This file implements the reader/parser used to read the k-vertex cut input data.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_KVC_H__
#define __SCIP_READER_KVC_H__


#include "scip/scip.h"


/** includes the k-vertex cut file reader into SCIP */
SCIP_RETCODE SCIPincludeReaderKvc(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
