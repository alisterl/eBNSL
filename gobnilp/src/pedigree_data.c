/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2017 James Cussens, Mark Barlett         *
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
 *  Contains functions related to managing PedigreeData.
 */

#include "pedigree_data.h"
#include "utils.h"

/** Deallocates the memory associated with a PedigreeData structure.
 *
 *  @param scip The SCIP instance which the data refers to.
 *  @param pd A pointer to the data structure to deallocate.
 *  @return SCIP_OKAY if the allocation succeeded.
 */
SCIP_RETCODE PE_deallocatePedigreeData(SCIP* scip, PedigreeData** pd)
{
   if( (*pd) != NULL )
   {
      if( (*pd)->SexVars != NULL )
         SCIPfreeMemoryArray(scip, &((*pd)->SexVars));
      if( (*pd)->ages != NULL )
         SCIPfreeMemoryArray(scip, &((*pd)->ages));
      if( (*pd)->sexes != NULL )
         SCIPfreeMemoryArray(scip, &((*pd)->sexes));
      SCIPfreeMemory(scip, pd);
   }
   return SCIP_OKAY;
}
/** Makes a deep copy of a PedigreeData structure.
 *
 *  @param scip The SCIP instance to which the data belongs.
 *  @param original The original data structure.
 *  @param duplicate A pointer to the duplicated data structure.
 *  @return SCIP_OKAY if copying suceeded.
 */
SCIP_RETCODE PE_copyPedigreeData(SCIP* scip, PedigreeData* original, PedigreeData** duplicate)
{
   int i;

   /* Allocate the memory */
   SCIP_CALL( SCIPallocMemory(scip, duplicate) );
   (*duplicate)->SexVars = NULL;
   (*duplicate)->ages = NULL;
   (*duplicate)->sexes = NULL;
   if( original->SexVars != NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->SexVars), original->n) );
   if( original->ages != NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->ages),    original->n) );
   if( original->sexes != NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->sexes),   original->n) );

   /* Copy the data across */
   (*duplicate)->n = original->n;
   for( i = 0; i < original->n; i++ )
   {
      if( original->SexVars != NULL )
         (*duplicate)->SexVars[i] = original->SexVars[i];
      if( original->ages != NULL )
         (*duplicate)->ages[i] = original->ages[i];
      if( original->sexes != NULL )
         (*duplicate)->sexes[i] = original->sexes[i];
   }

   return SCIP_OKAY;
}

/** Writes a PedigreeData structure to file.
 *
 *  @param scip The SCIP instance the data belongs to.
 *  @param file The file to write to.
 *  @param pd The data to write.
 *  @return SCIP_OKAY if the writing succeeded.
 */
SCIP_RETCODE PE_writeToFile(SCIP* scip, FILE* file, PedigreeData* pd)
{
   if( pd != NULL )
   {
      SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
      SCIPinfoMessage(scip, file, "%d", pd->n);
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      if( pd->SexVars != NULL )
         SCIP_CALL( UT_writeVarArray(scip, file, pd->SexVars, pd->n) );
      else
      {
         int i;
         SCIP_VAR** tmp_array;
         SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, pd->n) );
         for( i = 0; i < pd-> n; i++ )
            tmp_array[i] = NULL;
         SCIP_CALL( UT_writeVarArray(scip, file, tmp_array, pd->n) );
         SCIPfreeMemoryArray(scip, &tmp_array);
      }
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      if( pd->ages != NULL )
         SCIP_CALL( UT_writeIntArray(scip, file, pd->ages, pd->n) );
      else
      {
         int i;
         int* tmp_array;
         SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, pd->n) );
         for( i = 0; i < pd-> n; i++ )
            tmp_array[i] = -1;
         SCIP_CALL( UT_writeIntArray(scip, file, tmp_array, pd->n) );
         SCIPfreeMemoryArray(scip, &tmp_array);
      }
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      if( pd->sexes != NULL )
         SCIP_CALL( UT_writeCharArray(scip, file, pd->sexes, pd->n) );
      else
      {
         int i;
         char* tmp_array;
         SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, pd->n) );
         for( i = 0; i < pd-> n; i++ )
            tmp_array[i] = 'U';
         SCIP_CALL( UT_writeCharArray(scip, file, tmp_array, pd->n) );
         SCIPfreeMemoryArray(scip, &tmp_array);
      }
      SCIPinfoMessage(scip, file, "%c", UT_LIST_END);
   }
   return SCIP_OKAY;
}
/** Parses a PedigreeData structure from a sting.
 *
 *  @param scip The SCIP instance the data will belong to.
 *  @param str The string to parse.
 *  @param pd A pointer to the data structure resulting from parsing.
 *  @return SCIP_OKAY if parsing succeeded.
 */
SCIP_RETCODE PE_parse(SCIP* scip, char* str, PedigreeData** pd)
{
   int i;
   char** all_data = NULL;
   SCIP_Bool empty_sex = TRUE;
   SCIP_Bool empty_age = TRUE;
   SCIP_Bool unknown_sex = TRUE;

   int n;
   SCIP_VAR** SexVars = NULL;
   int* ages = NULL;
   char* sexes = NULL;

   /* Get the data as string arrays */
   SCIP_CALL( SCIPallocMemoryArray(scip, &all_data, 4) );
   for( i = 0; i < 4; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(all_data[i]), 1000000) );
   SCIP_CALL( UT_parseArray((char*)str, &all_data) );

   /* Find the number of nodes */
   sscanf(all_data[0], "%d", &n);

   /* Get the sex variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &SexVars, n) );
   SCIP_CALL( UT_parseVarArray(scip, all_data[1], &SexVars, n) );
   for( i = 0; i < n; i++ )
      if( SexVars[i] != NULL )
         empty_sex = FALSE;

   /* Get the ages of each individual */
   SCIP_CALL( SCIPallocMemoryArray(scip, &ages, n) );
   SCIP_CALL( UT_parseIntArray(all_data[2], &ages, n) );
   for( i = 0; i < n; i++ )
      if( ages[i] != -1 )
         empty_age = FALSE;

   /* Get the sexes of each individual */
   SCIP_CALL( SCIPallocMemoryArray(scip, &sexes, n) );
   SCIP_CALL( UT_parseCharArray(all_data[3], &sexes, n) );
   for( i = 0; i < n; i++ )
      if( sexes[i] != 'U' )
         unknown_sex = FALSE;

   /* Set the constraint data */
   SCIP_CALL( SCIPallocMemory(scip, pd) );
   (*pd)->n = n;
   if( !empty_sex )
      (*pd)->SexVars = SexVars;
   else
   {
      (*pd)->SexVars = NULL;
      SCIPfreeMemoryArray(scip, &SexVars);
   }
   if( !empty_age )
      (*pd)->ages = ages;
   else
   {
      (*pd)->ages = NULL;
      SCIPfreeMemoryArray(scip, &ages);
   }
   if( !unknown_sex )
      (*pd)->sexes = sexes;
   else
   {
      (*pd)->sexes = NULL;
      SCIPfreeMemoryArray(scip, &sexes);
   }

   /* Clean up */
   for( i = 0; i < 4; i++ )
      SCIPfreeMemoryArray(scip, &(all_data[i]));
   SCIPfreeMemoryArray(scip, &all_data);

   return SCIP_OKAY;
}

