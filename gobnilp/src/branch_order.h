
/**@file   branch_order.h
 * @ingroup BRANCHINGRULES
 * @brief  order branching rule
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_ORDER_H__
#define __SCIP_BRANCH_ORDER_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the order branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleOrder(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup BRANCHINGRULES
 *
 * @{
 */

/* TODO place other public methods in this group to facilitate navigation through the documentation */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
