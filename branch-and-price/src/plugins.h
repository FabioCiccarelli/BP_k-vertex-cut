/**@file   plugins.h
 * @brief  include SCIP plugins and custom plugins
 * @author Christopher Hojny
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PLUGINS_H__
#define __SCIP_PLUGINS_H__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** includes default SCIP plugins and custom plugins into SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePlugins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
