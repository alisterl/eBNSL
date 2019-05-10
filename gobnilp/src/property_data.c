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
 *  Contains functions related to managing PropertyData.
 */

#include "property_data.h"
#include "parent_set_data.h"
#include "utils.h"
#include <stdint.h>
#include <string.h>

/** Deallocates the memory associated with a PropertyData structure.
 *
 *  @param scip The SCIP instance which the data refers to.
 *  @param pd A pointer to the data structure to deallocate.
 *  @return SCIP_OKAY if the allocation succeeded.
 */
SCIP_RETCODE PR_deallocatePropertyData(
   SCIP* scip, 
   PropertyData** pd    
   )
{
   int i;
   int j;

   if( (*pd) != NULL )
   {
      if( (*pd)->global_property_names != NULL )
      {
         for( i = 0; i < (*pd)->num_global; i++ )
            if( (*pd)->global_property_names[i] != NULL )
               SCIPfreeMemoryArray(scip, &((*pd)->global_property_names[i]));
         SCIPfreeMemoryArray(scip, &((*pd)->global_property_names));
      }
      if( (*pd)->global_property_values != NULL )
      {
         for( i = 0; i < (*pd)->num_global; i++ )
            if( (*pd)->global_property_values[i] != NULL )
               SCIPfreeMemoryArray(scip, &((*pd)->global_property_values[i]));
         SCIPfreeMemoryArray(scip, &((*pd)->global_property_values));
      }
      if( (*pd)->num_properties != NULL )
      {
         if( (*pd)->property_names != NULL )
         {
            for( i = 0; i < (*pd)->n; i++ )
               if( (*pd)->num_properties[i] > 0 )
               {
                  for( j = 0; j < (*pd)->num_properties[i]; j++ )
                     if( (*pd)->property_names[i][j] != NULL )
                        SCIPfreeMemoryArray(scip, &((*pd)->property_names[i][j]));
                  SCIPfreeMemoryArray(scip, &((*pd)->property_names[i]));
               }
            SCIPfreeMemoryArray(scip, &((*pd)->property_names));
         }
         if( (*pd)->property_values != NULL )
         {
            for( i = 0; i < (*pd)->n; i++ )
               if( (*pd)->num_properties[i] > 0 )
               {
                  for( j = 0; j < (*pd)->num_properties[i]; j++ )
                     if( (*pd)->property_values[i][j] != NULL )
                        SCIPfreeMemoryArray(scip, &((*pd)->property_values[i][j]));
                  SCIPfreeMemoryArray(scip, &((*pd)->property_values[i]));
               }
            SCIPfreeMemoryArray(scip, &((*pd)->property_values));
         }
         SCIPfreeMemoryArray(scip, &((*pd)->num_properties));
      }
      SCIPfreeMemory(scip, pd);
   }
   return SCIP_OKAY;
}

/** Makes a deep copy of a PropertyData structure.
 *
 *  @param scip The SCIP instance to which the data belongs.
 *  @param original The original data structure.
 *  @param duplicate A pointer to the duplicated data structure.
 *  @return SCIP_OKAY if copying suceeded.
 */
SCIP_RETCODE PR_copyPropertyData(
   SCIP* scip, 
   PropertyData* original, 
   PropertyData** duplicate
   )
{
   int i;
   int j;

   /* Allocate the memory */
   SCIP_CALL( SCIPallocMemory(scip, duplicate) );
   (*duplicate)->global_property_names = NULL;
   (*duplicate)->global_property_values = NULL;
   (*duplicate)->num_properties = NULL;
   (*duplicate)->property_names = NULL;
   (*duplicate)->property_values = NULL;
   if( original->global_property_names != NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->global_property_names), original->num_global) );
   if( original->global_property_values != NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->global_property_values), original->num_global) );
   if( original->num_properties != NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->num_properties), original->n) );
   if( original->property_names != NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->property_names), original->n) );
   if( original->property_values != NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->property_values), original->n) );
   for( i = 0; i < original->n; i++ )
      if( original->num_properties[i] > 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->property_names[i]), original->num_properties[i]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*duplicate)->property_values[i]), original->num_properties[i]) );
      }
      else
      {
         (*duplicate)->property_names[i] = NULL;
         (*duplicate)->property_values[i] = NULL;
      }

   /* Copy the data across */
   (*duplicate)->num_global = original->num_global;
   for( i = 0; i < original->num_global; i++ )
   {
      if( original->global_property_names != NULL )
         (*duplicate)->global_property_names[i] = SCIPstrdup(scip, original->global_property_names[i]);
      if( original->global_property_values != NULL )
         (*duplicate)->global_property_values[i] = SCIPstrdup(scip, original->global_property_values[i]);
   }
   (*duplicate)->n = original->n;
   for( i = 0; i < original->n; i++ )
   {
      (*duplicate)->num_properties[i] = original->num_properties[i];
      for( j = 0; j < original->num_properties[i]; j++ )
      {
         if( original->property_names[i][j] != NULL )
            (*duplicate)->property_names[i][j] = SCIPstrdup(scip, original->property_names[i][j]);
         if( original->property_values[i][j] != NULL )
            (*duplicate)->property_values[i][j] = SCIPstrdup(scip, original->property_values[i][j]);
      }
   }

   return SCIP_OKAY;
}

/** Writes a PropertyData structure to file.
 *
 *  @param scip The SCIP instance the data belongs to.
 *  @param file The file to write to.
 *  @param pd The data to write.
 *  @return SCIP_OKAY if the writing succeeded.
 */
SCIP_RETCODE PR_writeToFile(
   SCIP* scip, 
   FILE* file, 
   PropertyData* pd
   )
{
   int i;

   if( pd != NULL )
   {
      SCIPinfoMessage(scip, file, "%c", UT_LIST_START);

      /* Global */
      SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
      SCIPinfoMessage(scip, file, "%d", pd->num_global);
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      if( pd->global_property_names != NULL )
         SCIP_CALL( UT_writeStringArray(scip, file, pd->global_property_names, pd->num_global) );
      else
         SCIPinfoMessage(scip, file, "%c%c", UT_LIST_START, UT_LIST_END);
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      if( pd->global_property_values != NULL )
         SCIP_CALL( UT_writeStringArray(scip, file, pd->global_property_values, pd->num_global) );
      else
         SCIPinfoMessage(scip, file, "%c%c", UT_LIST_START, UT_LIST_END);
      SCIPinfoMessage(scip, file, "%c", UT_LIST_END);

      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);

      /* Individual */
      SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
      SCIPinfoMessage(scip, file, "%d", pd->n);
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
      for( i = 0; i < pd->n; i++ )
      {
         if( i != 0 )
            SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
         SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
         SCIPinfoMessage(scip, file, "%d", pd->num_properties[i]);
         SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
         if( pd->num_properties[i] > 0 )
            SCIP_CALL( UT_writeStringArray(scip, file, pd->property_names[i], pd->num_properties[i]) );
         else
            SCIPinfoMessage(scip, file, "%c%c", UT_LIST_START, UT_LIST_END);
         SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
         if( pd->num_properties[i] > 0 )
            SCIP_CALL( UT_writeStringArray(scip, file, pd->property_values[i], pd->num_properties[i]) );
         else
            SCIPinfoMessage(scip, file, "%c%c", UT_LIST_START, UT_LIST_END);
         SCIPinfoMessage(scip, file, "%c", UT_LIST_END);
      }
      SCIPinfoMessage(scip, file, "%c", UT_LIST_END);
      SCIPinfoMessage(scip, file, "%c", UT_LIST_END);

      SCIPinfoMessage(scip, file, "%c", UT_LIST_END);
   }
   return SCIP_OKAY;
}
/** Parses a PropertyData structure from a sting.
 *
 *  @param scip The SCIP instance the data will belong to.
 *  @param str The string to parse.
 *  @param pd A pointer to the data structure resulting from parsing.
 *  @return SCIP_OKAY if parsing succeeded.
 */
SCIP_RETCODE PR_parse(
   SCIP* scip, 
   char* str, 
   PropertyData** pd
   )
{
   int i;
   int j;
   char** all_data = NULL;
   char** global_data = NULL;
   char** individual_data = NULL;
   char** item_data = NULL;
   char** item2_data = NULL;

   int num_global;
   char** global_property_names = NULL;
   char** global_property_values = NULL;

   int n;
   int* num_properties = NULL;
   char*** property_names = NULL;
   char*** property_values = NULL;

   /* Get the data as string arrays */
   SCIP_CALL( SCIPallocMemoryArray(scip, &all_data, 2) );
   for( i = 0; i < 2; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(all_data[i]), 1000000) );
   SCIP_CALL( UT_parseArray((char*)str, &all_data) );

   /* Get the global data as string arrays */
   SCIP_CALL( SCIPallocMemoryArray(scip, &global_data, 3) );
   for( i = 0; i < 3; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(global_data[i]), 1000000) );
   SCIP_CALL( UT_parseArray(all_data[0], &global_data) );

   /* Extract the global properties */
   num_global = atoi(global_data[0]);
   if( num_global > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &global_property_names, num_global) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &global_property_values, num_global) );
      for( i = 0; i < num_global; i++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(global_property_names[i]), SCIP_MAXSTRLEN) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(global_property_values[i]), SCIP_MAXSTRLEN) );
      }
      fprintf(stderr, "%s %d\n", global_data[1], num_global);
      SCIP_CALL( UT_parseStringArray(global_data[1], &global_property_names, num_global) );
      SCIP_CALL( UT_parseStringArray(global_data[2], &global_property_values, num_global) );
   }

   /* Get the individual data as string arrays */
   SCIP_CALL( SCIPallocMemoryArray(scip, &individual_data, 2) );
   for( i = 0; i < 2; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(individual_data[i]), 1000000) );
   SCIP_CALL( UT_parseArray(all_data[1], &individual_data) );

   /* Extract the individual properties */
   n = atoi(individual_data[0]);
   if( n > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &num_properties, n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &property_names, n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &property_values, n) );

      SCIP_CALL( SCIPallocMemoryArray(scip, &item_data, n) );
      for( i = 0; i < n; i++ )
         SCIP_CALL( SCIPallocMemoryArray(scip, &(item_data[i]), 1000000) );
      SCIP_CALL( UT_parseArray(individual_data[1], &item_data) );

      /* Extract each item's data */
      for( i = 0; i < n; i++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &item2_data, 3) );
         for( j = 0; j < 3; j++ )
            SCIP_CALL( SCIPallocMemoryArray(scip, &(item2_data[j]), 1000000) );
         SCIP_CALL( UT_parseArray(item_data[i], &item2_data) );

         num_properties[i] = atoi(item2_data[0]);
         if( num_properties[i] > 0 )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(property_names[i]), num_properties[i]) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &(property_values[i]), num_properties[i]) );
            for( j = 0; j < num_properties[i]; j++ )
            {
               SCIP_CALL( SCIPallocMemoryArray(scip, &(property_names[i][j]), SCIP_MAXSTRLEN) );
               SCIP_CALL( SCIPallocMemoryArray(scip, &(property_values[i][j]), SCIP_MAXSTRLEN) );
            }
            SCIP_CALL( UT_parseStringArray(item2_data[1], &(property_names[i]), num_properties[i]) );
            SCIP_CALL( UT_parseStringArray(item2_data[2], &(property_values[i]), num_properties[i]) );
         }

         for( j = 0; j < 3; j++ )
            SCIPfreeMemoryArray(scip, &(item2_data[j]));
         SCIPfreeMemoryArray(scip, &item2_data);
      }
   }

   /* Store the results */
   SCIP_CALL( SCIPallocMemory(scip, pd) );
   (*pd)->num_global = num_global;
   (*pd)->global_property_names = global_property_names;
   (*pd)->global_property_values = global_property_values;
   (*pd)->n = n;
   (*pd)->num_properties = num_properties;
   (*pd)->property_names = property_names;
   (*pd)->property_values = property_values;

   /* Clean up */
   for( i = 0; i < 2; i++ )
      SCIPfreeMemoryArray(scip, &(all_data[i]));
   SCIPfreeMemoryArray(scip, &(all_data));
   for( i = 0; i < 3; i++ )
      SCIPfreeMemoryArray(scip, &(global_data[i]));
   SCIPfreeMemoryArray(scip, &(global_data));
   for( i = 0; i < 2; i++ )
      SCIPfreeMemoryArray(scip, &(individual_data[i]));
   SCIPfreeMemoryArray(scip, &(individual_data));
   for( i = 0; i < n; i++ )
      SCIPfreeMemoryArray(scip, &(item_data[i]));
   SCIPfreeMemoryArray(scip, &item_data);

   return SCIP_OKAY;
}

/** Initialises a data structure with empty values.
 *  @param pd The data structure to initialise.
 */
void PR_initialise(
   PropertyData* pd
   )
{
   pd->num_global = 0;
   pd->global_property_names = NULL;
   pd->global_property_values = NULL;
   pd->n = 0;
   pd->num_properties = NULL;
   pd->property_names = NULL;
   pd->property_values = NULL;
}

/** Find out whether a string is a global property 
* @return TRUE if the string is a global property, else FALSE */
SCIP_Bool PR_hasGlobalProperty(
   SCIP* scip,        /**< SCIP data structure (not used) */
   PropertyData* pd,  /**< Property data structure */
   const char* name   /**< String to be queried */
   )
{
   int i;
   for( i = 0; i < pd->num_global; i++ )
      if( strcmp(name, pd->global_property_names[i]) == 0 )
         return TRUE;
   return FALSE;
}

/** Get the value of a global property 
* @return If the global property exists then its value else NULL */
char* PR_getGlobalProperty(
   SCIP* scip,         /**< SCIP data structure (not used) */
   PropertyData* pd,   /**< Property data structure */
   const char* name    /**< Global property whose value is sought */
   )
{
   int i;
   for( i = 0; i < pd->num_global; i++ )
      if( strcmp(name, pd->global_property_names[i]) == 0 )
         return pd->global_property_values[i];
   return NULL;
}

/** Set a global property to a given value 
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setGlobalProperty(
   SCIP* scip,           /**< SCIP data structure */ 
   PropertyData* pd,     /**< Property data structure */
   const char* name,     /**< Global property to be set */
   const char* value     /**< Value of global property */
   )
{
   int i;
   for( i = 0; i < pd->num_global; i++ )
      if( strcmp(name, pd->global_property_names[i]) == 0 )
      {
         SCIPfreeMemoryArray(scip, &(pd->global_property_values[i]));
         pd->global_property_values[i] = SCIPstrdup(scip, value);
         return SCIP_OKAY;
      }

   if( pd->num_global == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pd->global_property_names), 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pd->global_property_values), 1) );
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(pd->global_property_names), pd->num_global + 1) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(pd->global_property_values), pd->num_global + 1) );
   }

   pd->global_property_names[pd->num_global] = SCIPstrdup(scip, name);
   pd->global_property_values[pd->num_global] = SCIPstrdup(scip, value);
   pd->num_global += 1;

   return SCIP_OKAY;
}

/** Set a global property to a given value. 
 * Value is converted from an int to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setGlobalPropertyFromInt(
   SCIP* scip,         /**< SCIP data structure */ 
   PropertyData* pd,   /**< Property data structure */
   const char* name,   /**< Global property to be set */
   int value           /**< Value for the global property */
   )
{
   char tmp_string[SCIP_MAXSTRLEN];
   sprintf(tmp_string, "%d", value);
   SCIP_CALL( PR_setGlobalProperty(scip, pd, name, tmp_string) );
   return SCIP_OKAY;
}

/** Set a global property to a given value 
 * Value is converted from a SCIP_Real to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setGlobalPropertyFromReal(
   SCIP* scip,         /**< SCIP data structure */ 
   PropertyData* pd,   /**< Property data structure */
   const char* name,   /**< Global property to be set */
   SCIP_Real value     /**< Value for the global property */
   )
{
   char tmp_string[SCIP_MAXSTRLEN];
   sprintf(tmp_string, "%f", value);
   SCIP_CALL( PR_setGlobalProperty(scip, pd, name, tmp_string) );
   return SCIP_OKAY;
}

/** Set a global property to a given value 
 * Value is converted from a SCIP_Bool to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setGlobalPropertyFromBool(
   SCIP* scip,          /**< SCIP data structure */
   PropertyData* pd,    /**< Property data structure */
   const char* name,    /**< Global property to be set */
   SCIP_Bool value      /**< Value for the global property */
   )
{
   if( value )
      SCIP_CALL( PR_setGlobalProperty(scip, pd, name, "TRUE") );
   else
      SCIP_CALL( PR_setGlobalProperty(scip, pd, name, "FALSE") );
   return SCIP_OKAY;
}

/** Set a global property to a given value 
 * Value is converted from an array of strings to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setGlobalPropertyFromArray(
   SCIP* scip,           /**< SCIP data structure */
   PropertyData* pd,     /**< Property data structure */
   const char* name,     /**< Global property to be set */
   const char** value,   /**< Value for the global property */
   int length            /**< Length of the array */
   )
{
   int i;
   char tmp_value[SCIP_MAXSTRLEN];
   tmp_value[0] = '\0';
   for( i = 0; i < length; i++ )
   {
      if( i > 0 )
         strcat(tmp_value, " ");
      strcat(tmp_value, value[i]);
   }
   SCIP_CALL( PR_setGlobalProperty(scip, pd, name, tmp_value) );
   return SCIP_OKAY;
}

/** Set a global property to a given value 
 * Value is converted from an array of reals to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setGlobalPropertyFromRealArray(
   SCIP* scip,         /**< SCIP data structure */ 
   PropertyData* pd,   /**< Property data structure */
   const char* name,   /**< Global property to be set */
   SCIP_Real* value,   /**< Value for the global property */
   int length          /**< Length of the array */
   )
{
   int i;
   char tmp_value[SCIP_MAXSTRLEN];
   tmp_value[0] = '\0';
   for( i = 0; i < length; i++ )
   {
      char next_value[SCIP_MAXSTRLEN];
      if( i > 0 )
         strcat(tmp_value, " ");
      sprintf(next_value, "%f", value[i]);
      strcat(tmp_value, next_value);
   }
   SCIP_CALL( PR_setGlobalProperty(scip, pd, name, tmp_value) );
   return SCIP_OKAY;
}

/** Does a property data structure have a given property for a given individual?
* @return TRUE if it does, FALSE otherwise
*/
SCIP_Bool PR_hasProperty(
   SCIP* scip,          /**< SCIP data structure */ 
   PropertyData* pd,    /**< Property data structure */
   int individual,      /**< Individual */
   const char* name     /**< Property name */
   )
{
   int i;
   for( i = 0; i < pd->num_properties[individual]; i++ )
      if( strcmp(name, pd->property_names[individual][i]) == 0 )
         return TRUE;
   return FALSE;
}

/** Returns the value of a given property for a given individual
 * @return The value of the property for the individual or NULL if missing
*/
char* PR_getProperty(
   SCIP* scip,         /**< SCIP data structure */ 
   PropertyData* pd,   /**< Property data structure */
   int individual,     /**< Individual */
   const char* name    /**< Property name */
   )
{
   int i;
   for( i = 0; i < pd->num_properties[individual]; i++ )
      if( strcmp(name, pd->property_names[individual][i]) == 0 )
         return pd->property_values[individual][i];
   return NULL;
}
/** Set a property for an individual to a given value 
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setProperty(
   SCIP* scip,         /**< SCIP data structure */ 
   PropertyData* pd,   /**< Property data structure */
   int individual,     /**< Individual */
   const char* name,   /**< Property name */
   const char* value   /**< Property value */
   )
{
   int i;

   if( pd->num_properties == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pd->num_properties), pd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pd->property_names), pd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pd->property_values), pd->n) );
      for( i = 0; i < pd->n; i++ )
      {
         pd->num_properties[i] = 0;
         pd->property_names[i] = NULL;
         pd->property_values[i] = NULL;
      }
   }

   for( i = 0; i < pd->num_properties[individual]; i++ )
      if( strcmp(name, pd->property_names[individual][i]) == 0 )
      {
         SCIPfreeMemoryArray(scip, &(pd->property_values[individual][i]));
         pd->property_values[individual][i] = SCIPstrdup(scip, value);
         return SCIP_OKAY;
      }

   if( pd->num_properties[individual] == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pd->property_names[individual]), 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pd->property_values[individual]), 1) );
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(pd->property_names[individual]), pd->num_properties[individual] + 1) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(pd->property_values[individual]), pd->num_properties[individual] + 1) );
   }

   pd->property_names[individual][pd->num_properties[individual]] = SCIPstrdup(scip, name);
   pd->property_values[individual][pd->num_properties[individual]] = SCIPstrdup(scip, value);
   pd->num_properties[individual] += 1 ;

   return SCIP_OKAY;
}

/** Set a property for an individual to a given value 
 * Value is converted from an int to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setPropertyFromInt(
   SCIP* scip,         /**< SCIP data structure */
   PropertyData* pd,   /**< Property data structure */
   int individual,     /**< Individual */             
   const char* name,   /**< Property name */          
   int value           /**< Property value */          
   )
{
   char tmp_string[SCIP_MAXSTRLEN];
   sprintf(tmp_string, "%d", value);
   SCIP_CALL( PR_setProperty(scip, pd, individual, name, tmp_string) );
   return SCIP_OKAY;
}

/** Set a property for an individual to a given value 
 * Value is converted from a SCIP_Real to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setPropertyFromReal(
   SCIP* scip,     /**< SCIP data structure */ 
   PropertyData* pd,  /**< Property data structure */
   int individual,    /**< Individual */             
   const char* name,  /**< Property name */          
   SCIP_Real value    /**< Property value */         
   )
{
   char tmp_string[SCIP_MAXSTRLEN];
   sprintf(tmp_string, "%f", value);
   SCIP_CALL( PR_setProperty(scip, pd, individual, name, tmp_string) );
   return SCIP_OKAY;
}

/** Set a property for an individual to a given value 
 * Value is converted from a SCIP_Bool to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setPropertyFromBool(
   SCIP* scip,     /**< SCIP data structure */ 
   PropertyData* pd, /**< Property data structure */
   int individual,   /**< Individual */             
   const char* name, /**< Property name */          
   SCIP_Bool value   /**< Property value */         
   )
{
   if( value )
      SCIP_CALL( PR_setProperty(scip, pd, individual, name, "TRUE") );
   else
      SCIP_CALL( PR_setProperty(scip, pd, individual, name, "FALSE") );
   return SCIP_OKAY;
}

/** Set a property for an individual to a given value 
 * Value is converted from an array of strings to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setPropertyFromArray(
   SCIP* scip,   /**< SCIP data structure */ 
   PropertyData* pd, /**< Property data structure */
   int individual,   /**< Individual */             
   const char* name, /**< Property name */          
   const char** value/**< Property value */         , 
   int length        /**< Length of the array */
   )
{
   int i;
   char tmp_value[SCIP_MAXSTRLEN];
   tmp_value[0] = '\0';
   for( i = 0; i < length; i++ )
   {
      if( i > 0 )
         strcat(tmp_value, " ");
      strcat(tmp_value, value[i]);
   }
   SCIP_CALL( PR_setProperty(scip, pd, individual, name, tmp_value) );
   return SCIP_OKAY;
}
/** Set a property for an individual to a given value 
 * Value is converted from an array of reals to a string
 * @return SCIP_OKAY if all is well */
SCIP_RETCODE PR_setPropertyFromRealArray(
   SCIP* scip,        /**< SCIP data structure */ 
   PropertyData* pd,  /**< Property data structure */
   int individual,    /**< Individual */             
   const char* name,  /**< Property name */          
   SCIP_Real* value,  /**< Property value */         
   int length         /**< Length of the array */
   )
{
   int i;
   char tmp_value[SCIP_MAXSTRLEN];
   tmp_value[0] = '\0';
   for( i = 0; i < length; i++ )
   {
      char next_value[SCIP_MAXSTRLEN];
      if( i > 0 )
         strcat(tmp_value, " ");
      sprintf(next_value, "%f", value[i]);
      strcat(tmp_value, next_value);
   }
   SCIP_CALL( PR_setProperty(scip, pd, individual, name, tmp_value) );
   return SCIP_OKAY;
}
