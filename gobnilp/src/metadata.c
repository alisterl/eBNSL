/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2017 James Cussens, Mark Bartlett        *
 *                                                                       *
 *   This program is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU General Public License as      *
 *   published by the Free Software Foundation; either version 3 of the  *
 *   License, or (at your option) any later version.                     *
 *                                                                       *
 *   This program is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *
 *   General Public License for more details.                            *
 *                                                                       *
 *   You should have received a copy of the GNU General Public License   *
 *   along with this program; if not, see                                *
 *   <http://www.gnu.org/licenses>.                                      *
 *                                                                       *
 *   Additional permission under GNU GPL version 3 section 7             *
 *                                                                       *
 *   If you modify this Program, or any covered work, by linking or      *
 *   combining it with SCIP (or a modified version of that library),     *
 *   containing parts covered by the terms of the ZIB Academic License,  *
 *   the licensors of this Program grant you additional permission to    *
 *   convey the resulting work.                                          *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/** @file
 *  Contains metadata relating to the problem that is not necessarily part of the LP that
 *  will be solved.
 *
 *  At present, the information stored is that of what family sets exist, which variable
 *  is associated with each and the names of each node.
 *
 *  This information is stored as a constraint.  However, it puts no limitations on the
 *  problem and in fact should never be called.  The reason for storing this in a
 *  constraint is so that the information contained can be written to and read from .cip
 *  format files in order to preserve the information needed if the problem is saved in
 *  this way.
 */

#include <string.h>
#include "metadata.h"
#include "utils.h"

/** The name of this 'constraint handler'. */
#define CONSHDLR_NAME "metadata"
/** The name of the 'constraint' referring to parent set data. */
#define PARENT_SET_NAME "parent set data"
/** The name of the 'constraint' referring to pedigree data. */
#define PEDIGREE_NAME   "pedigree data"
/** The name of the 'constraint' referring to property data. */
#define PROPERTY_NAME   "property data"

struct SCIP_ConshdlrData
{
   /** The data on parent sets and variables that will be stored by the metadata. */
   ParentSetData* psd;
   /** The data specific to pedigrees that is stored as metadata. */
   PedigreeData* pd;
   /** The property data that is stored as metadata. */
   PropertyData* prop;
};

/** Frees the metadata held when the program exits. */
static SCIP_DECL_CONSFREE(consFreeMetadata)
{
   SCIP_CONSHDLRDATA* data = SCIPconshdlrGetData(conshdlr);
   if( data != NULL )
   {
      if( data->psd != NULL )
         PS_deallocateParentSetData(scip, &(data->psd), FALSE);
      if( data->pd != NULL )
         PE_deallocatePedigreeData(scip, &(data->pd));
      if( data->prop != NULL )
         PR_deallocatePropertyData(scip, &(data->prop));

      SCIPfreeBlockMemory(scip, &data);
   }

   SCIPconshdlrSetData(conshdlr, NULL);
   return SCIP_OKAY;
}

/* Reading and writing to files */
/** Prints the metadata to a .cip file as a fake constraint entry. */
static SCIP_DECL_CONSPRINT(consPrintMetadata)
{
   SCIP_CONSHDLRDATA* data = SCIPconshdlrGetData(conshdlr);
   const char* name = SCIPconsGetName(cons);
   if( strcmp(name, PARENT_SET_NAME) == 0 )
      PS_writeToFile(scip, file, data->psd);
   else if( strcmp(name, PEDIGREE_NAME) == 0 )
      PE_writeToFile(scip, file, data->pd);
   else if( strcmp(name, PROPERTY_NAME) == 0 )
      PR_writeToFile(scip, file, data->prop);
   return SCIP_OKAY;
}
/** Parses the metadata returned from a .cip file and stores it. */
static SCIP_DECL_CONSPARSE(consParseMetadata)
{
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   *success = FALSE;

   if( strcmp(name, PARENT_SET_NAME) == 0 )
   {
      SCIP_CALL( PS_parse(scip, (char*)str, &(conshdlrdata->psd)) );
      *success = TRUE;
   }
   else if( strcmp(name, PEDIGREE_NAME) == 0 )
   {
      SCIP_CALL( PE_parse(scip, (char*)str, &(conshdlrdata->pd)) );
      *success = TRUE;
   }
   else if( strcmp(name, PROPERTY_NAME) == 0 )
   {
      SCIP_CALL( PR_parse(scip, (char*)str, &(conshdlrdata->prop)) );
      *success = TRUE;
   }

   if( *success == TRUE )
      SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, NULL, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/* Fundamental callback methods not actually needed and will always succeed. */
/** Unneeded fundamental callback method, which will always succeed. */
static SCIP_DECL_CONSENFOLP(consEnfolpMetadata)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}
/** Unneeded fundamental callback method, which will always succeed. */
static SCIP_DECL_CONSENFOPS(consEnfopsMetadata)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}
/** Unneeded fundamental callback method, which will always succeed. */
static SCIP_DECL_CONSCHECK(consCheckMetadata)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}
/** Unneeded fundamental callback method, which will always succeed. */
static SCIP_DECL_CONSLOCK(consLockMetadata)
{
   return SCIP_OKAY;
}

/* Utility functions to make adding new meta data more stright forward. */
/** Creates a constraint with unusual settings that should mean it is never called.
 *
 *  @param scip The SCIP instance the constaint features in.
 *  @param name The name to give the constraint.
 *  @return SCIP_OKAY if the creation succeeds, or an appropriate error code otherwise.
 */
static SCIP_RETCODE createBlankConstraint(SCIP* scip, const char* name)
{
   SCIP_CONS* cons = NULL;
   SCIP_CALL(SCIPcreateCons(scip, &cons, name, SCIPfindConshdlr(scip, CONSHDLR_NAME), NULL,
                            FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   return SCIP_OKAY;
}
/** Gets the data associated with the constraint if it has already been set.
 *
 *  @param scip The SCIP instance the constraint features in.
 *  @param conshdlrdata A pointer to return the data.
 *  @return SCIP_OKAY if the data could be retrieved, or an appropriate error message otherwise.
 */
static SCIP_RETCODE getConshdlrDataIfExists(SCIP* scip, SCIP_CONSHDLRDATA** conshdlrdata)
{
   SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Need to initialise metadata before trying to access its data.\n");
      return SCIP_ERROR;
   }
   else
   {
      (*conshdlrdata) = SCIPconshdlrGetData(conshdlr);
      return SCIP_OKAY;
   }
}

/* Externally visible interface */
/** Initialised the metadata, so that information can be stored in it.
 *
 *  This function must be called before attempting to add data using MD_setParentSetData().
 *
 *  @param scip The SCIP instance to which the metadata will belong.
 *  @return SCIP_OKAY if initialisation succeeded or an appropriate error code otherwise.
 */
SCIP_RETCODE MD_initialiseMetadata(SCIP* scip)
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->psd  = NULL;
   conshdlrdata->pd   = NULL;
   conshdlrdata->prop = NULL;

   /* Include the handler in a nice clear way */
   SCIP_CALL(SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME,
                                      "pseudo-constraint handler, just used to read and UT_write to files", -1, -1, -1, TRUE,
                                      consEnfolpMetadata, consEnfopsMetadata, consCheckMetadata, consLockMetadata, conshdlrdata));
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeMetadata) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseMetadata) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintMetadata) );
   return SCIP_OKAY;
}
/** Sets the parent set data associated with the problem.
 *
 *  This method should only be called once, when the metadata first becomes available
 *  for storage.  If the information needs changing at a later date, call
 *  MD_getParentSetData() and modify the data structure directly.
 *
 *  @param scip The SCIP instance the metdata relates to.
 *  @param psd The parent set data for the problem
 *  @return SCIP_OKAY if storage worked, or an appropriate error message otherwise.
 */
SCIP_RETCODE MD_setParentSetData(SCIP* scip, ParentSetData* psd)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CALL( getConshdlrDataIfExists(scip, &conshdlrdata) );
   if( conshdlrdata->psd != NULL )
   {
      SCIPerrorMessage("Parent set data is already set.\n");
      return SCIP_ERROR;
   }
   else
   {
      SCIP_CALL( createBlankConstraint(scip, PARENT_SET_NAME) );
      SCIP_CALL( PS_copyParentSetData(scip, psd, &(conshdlrdata->psd)) );
      return SCIP_OKAY;
   }
}
/** Gets the parent set information associated with the problem.
 *
 *  Note that this will return the actual stored data structure, not a copy.
 *
 *  @param scip The SCIP instance the metadata belongs to.
 *  @return The parent set data stored for the problem.
 */
ParentSetData* MD_getParentSetData(SCIP* scip)
{
   SCIP_CONSHDLR* handler = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( handler == NULL )
      return NULL;
   else
      return (SCIPconshdlrGetData(handler))->psd;
}
/** Sets the pedigree data associated with the problem.
 *
 *  This method should only be called once, when the metadata first becomes available
 *  for storage.  If the information needs changing at a later date, call
 *  MD_getPedigreeData() and modify the data structure directly.
 *
 *  @param scip The SCIP instance the metdata relates to.
 *  @param pd The pedigree data for the problem
 *  @return SCIP_OKAY if storage worked, or an appropriate error message otherwise.
 */
SCIP_RETCODE MD_setPedigreeData(SCIP* scip, PedigreeData* pd)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CALL( getConshdlrDataIfExists(scip, &conshdlrdata) );
   if( conshdlrdata->pd != NULL )
   {
      SCIPerrorMessage("Pedigree data is already set.\n");
      return SCIP_ERROR;
   }
   else
   {
      SCIP_CALL( createBlankConstraint(scip, PEDIGREE_NAME) );
      SCIP_CALL( PE_copyPedigreeData(scip, pd, &(conshdlrdata->pd)) );
      return SCIP_OKAY;
   }
}
/** Gets the pedigree information associated with the problem.
 *
 *  Note that this will return the actual stored data structure, not a copy.
 *
 *  @param scip The SCIP instance the metadata belongs to.
 *  @return The pedigree data stored for the problem.
 */
PedigreeData* MD_getPedigreeData(SCIP* scip)
{
   SCIP_CONSHDLR* handler = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( handler == NULL )
      return NULL;
   else
      return (SCIPconshdlrGetData(handler))->pd;
}
/** Sets the property data associated with the problem.
 *
 *  This method should only be called once, when the metadata first becomes available
 *  for storage.  If the information needs changing at a later date, call
 *  MD_getPropertyData() and modify the data structure directly.
 *
 *  @param scip The SCIP instance the metdata relates to.
 *  @param prop The property data for the problem
 *  @return SCIP_OKAY if storage worked, or an appropriate error message otherwise.
 */
SCIP_RETCODE MD_setPropertyData(SCIP* scip, PropertyData* prop)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CALL( getConshdlrDataIfExists(scip, &conshdlrdata) );
   if( conshdlrdata->prop != NULL )
   {
      SCIPerrorMessage("Property data is already set.\n");
      return SCIP_ERROR;
   }
   else
   {
      SCIP_CALL( createBlankConstraint(scip, PROPERTY_NAME) );
      SCIP_CALL( PR_copyPropertyData(scip, prop, &(conshdlrdata->prop)) );
      return SCIP_OKAY;
   }
}
/** Gets the property information associated with the problem.
 *
 *  Note that this will return the actual stored data structure, not a copy.
 *
 *  @param scip The SCIP instance the metadata belongs to.
 *  @return The property data stored for the problem.
 */
PropertyData* MD_getPropertyData(SCIP* scip)
{
   SCIP_CONSHDLR* handler = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( handler == NULL )
      return NULL;
   else
      return (SCIPconshdlrGetData(handler))->prop;
}
