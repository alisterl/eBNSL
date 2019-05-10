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
 *  Contains some useful functions that may be needed in several other files.
 *
 *  Currently, there are four sorts of function included in this file.
 *  - Functions to simplify creating parameters.
 *    - UT_addBoolParam()
 *    - UT_addIntParam()
 *    - UT_addRealParam()
 *    - UT_addStringParam()
 *  - Functions to simplify creating linear constraints.
 *    - UT_createEmptyLinearConstraint()
 *    - UT_createEmptyLTEConstraint()
 *    - UT_createEmptyGTEConstraint()
 *  - Functions to help read and write data to files.
 *    - UT_writeStringArray()
 *    - UT_writeIntArray()
 *    - UT_writeVarArray()
 *    - UT_writeIntArrayArray()
 *    - UT_writeVarArrayArray()
 *    - UT_writeIntArrayArrayArray()
 *    - UT_parseArray()
 *    - UT_parseStringArray()
 *    - UT_parseIntArray()
 *    - UT_parseVarArray()
 *    - UT_parseIntArrayArray()
 *    - UT_parseVarArrayArray()
 *    - UT_parseIntArrayArrayArray()
 *  - Functions to create and use a hash table for storing 'arrow' and 'edge' variables
 *  - get_index
 */

#include <assert.h>
#include "utils.h"
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#define HASHMAP_SIZE 20

/** Gets the index of a named variable.
 *  Adapted from function of same name in probdata_bn  
 *  @return the index of the given variable in the given parent sets data structure
 */
int get_index(
   char* nodeName,     /**< (the name of ) the variable */
   ParentSetData* psd  /**< the parent sets data structure, only nodeNames is set */
   )
{
   return get_index_names(nodeName, psd->n, psd->nodeNames);
}



/** Gets the index of a named variable.
 *  Adapted from function of same name in probdata_bn  
 *  @return the index of the given variable in the given parent sets data structure
 */
int get_index_names(
   char* nodeName,     /**< (the name of ) the variable */
   int n,
   char** nodeNames
   )
{
   int i;
   char* end;

   /* from stackoverflow 122616 */

   /* Trim leading space */
   while( isspace(*nodeName) ) 
      nodeName++;

   /* Trim trailing space */
   end = nodeName + strlen(nodeName) - 1;
   while( end > nodeName && isspace(*end) ) 
      end--;

   /* Write new null terminator */
   *(end+1) = 0;

   for( i = 0; i < n; ++i )
      if( strcmp(nodeName, nodeNames[i]) == 0 )
         return i;
   SCIPerrorMessage("Not recognised as a variable name: %s\n", nodeName);
   return -1;
}


SCIP_RETCODE free_allsubsets(
   SCIP* scip,             /**< SCIP instance */
   int*** subsets_ptr      /**< data to free */
)
{
   int i;
   const int n = (*subsets_ptr)[0][0];
   
   assert(subsets_ptr != NULL);
   assert(*subsets_ptr != NULL);
   assert((*subsets_ptr)[0] != NULL);
   
   for( i = 0; i <= n; ++i )
      SCIPfreeMemoryArray(scip,&((*subsets_ptr)[i]));
   SCIPfreeMemoryArray(scip, subsets_ptr);
   return SCIP_OKAY;
}

SCIP_RETCODE allsubsets(
   SCIP* scip,             /**< SCIP instance */
   const int* set,         /**< set whose subsets we want */
   const int n,            /**< number of elements in set */
   int*** subsets_ptr      /**< (pointer to ) output 
                            subsets[0][0] is the number of subsets 
                            subsets[j][0] is the size of the jth subset
                            subsets[j][i] is the ith element of the jth subset */
   )
{

   int old_nsubsets;
   int new_nsubsets;
   int new_one;
   int i;
   int j;
   int k;
   int j_size;

   assert(scip != NULL);
   assert(set != NULL);
   assert(subsets_ptr != NULL);

   if( n == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, subsets_ptr, 2) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*subsets_ptr)[0]), 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*subsets_ptr)[1]), 1) );
      (*subsets_ptr)[0][0] = 1;         /* number of subsets */
      (*subsets_ptr)[1][0] = 0;         /* size of empty set */
      return SCIP_OKAY;
   }

   if( n == 1 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, subsets_ptr, 3) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*subsets_ptr)[0]), 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*subsets_ptr)[1]), 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*subsets_ptr)[2]), 2) );
      (*subsets_ptr)[0][0] = 2;         /* number of subsets */
      (*subsets_ptr)[1][0] = 0;         /* size of empty set */
      (*subsets_ptr)[2][0] = 1;         /* size of singleton set */
      (*subsets_ptr)[2][1] = set[0];    /* singleton set */
      return SCIP_OKAY;
   }

   SCIP_CALL( allsubsets(scip, set+1, n-1, subsets_ptr) );

   new_one = set[0];
   old_nsubsets = (*subsets_ptr)[0][0];
   new_nsubsets = 2*old_nsubsets;
   /* make room for subsets which include set[0] */
   SCIP_CALL( SCIPreallocMemoryArray(scip, subsets_ptr, new_nsubsets+1) );

   for( i = 1; i <= old_nsubsets; i++ )
   {
      j = i + old_nsubsets;
      j_size = (*subsets_ptr)[i][0] + 1;
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*subsets_ptr)[j]), j_size+1) );
      (*subsets_ptr)[j][0] = j_size;
      for( k = 1; k < j_size; k++ )
         (*subsets_ptr)[j][k] = (*subsets_ptr)[i][k];
      (*subsets_ptr)[j][k] = new_one; 
   }
   (*subsets_ptr)[0][0] = new_nsubsets;
   return SCIP_OKAY;
}

/** Frees memory (previously) used by a hash table */
SCIP_RETCODE hashtablefreeArrow(
   SCIP* scip,             /**< SCIP instance */
   ParentSetData* psd      /**< Parent set data structure containing hashtable */
   )
{

   assert(scip != NULL);
   assert(psd != NULL);

   SCIPhashmapFree(&(psd->arrow));
   SCIPhashmapFree(&(psd->edge));

   return SCIP_OKAY;
}
   
/** Creates a hash table for arrow variables and a hash table for edge variables */
SCIP_RETCODE hashtableCreateArrow(
   SCIP* scip,          /**< SCIP instance */
   ParentSetData* psd   /**< Parent set data structure containing hashtable */
   )
{

   assert(scip != NULL);
   assert(psd != NULL);

   SCIP_CALL( SCIPhashmapCreate(&(psd->arrow), SCIPblkmem(scip), HASHMAP_SIZE) );
   SCIP_CALL( SCIPhashmapCreate(&(psd->edge), SCIPblkmem(scip), HASHMAP_SIZE) );
   
   return SCIP_OKAY;
}

/** Retrieve an arrow or edge variable for specified child and parent or NULL if none exists 
 * @return an arrow or edge variable for specified child and parent or NULL if none exists 
*/
SCIP_VAR* get_arrowedge(
   const ParentSetData* psd,  /**< the parent sets data structure */
   int i,                     /**< child variable */
   int j,                     /**< parent variable */
   SCIP_Bool arrow            /**< If TRUE return arrow variable for \a i and \a j, else return edge variable for \a i and \a j */ 
   )
{
   SCIP_HASHMAP* hashmap;
   int l;
   int idx;

   if( j == i )
      return NULL;

   if( arrow )
   {
      hashmap = psd->arrow;
   }
   else
   {
      hashmap = psd->edge;
      if( i > j )
      {
         l = i;
         i = j;
         j = l;
      }
   }

   assert(hashmap != NULL);
   idx = (psd->n)*i + j;
   return (SCIP_VAR*) SCIPhashmapGetImage(hashmap, (void*)(size_t)idx);
}

/** Retrieve an arrow variable for specified child and parent or NULL if none exists 
 * @return an arrow variable for specified child and parent or NULL if none exists
 */ 
SCIP_VAR* get_arrow(
   const ParentSetData* psd,  /**< the parent sets data structure */
   const int i,               /**< child variable */
   const int j                /**< parent variable */
   )
{
   return get_arrowedge(psd, i, j, TRUE);
}

/** Retrieve an edge variable for specified child and parent or NULL if none exists 
 * @return an edge variable for specified child and parent or NULL if none exists
*/
SCIP_VAR* get_edge(
   const ParentSetData* psd,  /**< the parent sets data structure */
   const int i,               /**< child variable */
   const int j                /**< parent variable */
   )
{
   return get_arrowedge(psd, i, j, FALSE);
}


/** Store an arrow or edge variable for specified child and parent */
static
SCIP_RETCODE put_arrowedge(
   SCIP* scip,          /**< SCIP instance */
   ParentSetData* psd,  /**< the parent sets data structure */
   int i,               /**< child variable */
   int j,               /**< parent variable */
   SCIP_VAR* arrow_i_j, /**< variable to store */
   SCIP_Bool arrow      /**< If TRUE store arrow variable for \a i and \a j, else store edge variable for \a i and \a j */ 
   )
{

   SCIP_HASHMAP* hashmap;
   int l;
   int idx;

   assert(scip != NULL);
   assert(psd != NULL);
   
   if( i == j )
      return SCIP_ERROR;

   if( arrow )
   {
      hashmap = psd->arrow;
   }
   else
   {
      hashmap = psd->edge;
      if( i > j )
      {
         l = i;
         i = j;
         j = l;
      }
   }
   assert(hashmap != NULL);
   assert(i < (psd->n));
   assert(j < (psd->n));
   idx = (psd->n) * i + j;
   SCIP_CALL( SCIPhashmapInsert(hashmap, (void*)(size_t)idx, arrow_i_j) );

   return SCIP_OKAY;
}

/** Store an arrow variable for specified child and parent */
SCIP_RETCODE put_arrow(
   SCIP* scip,          /**< SCIP instance */
   ParentSetData* psd,  /**< the parent sets data structure */
   const int i,         /**< child variable */
   const int j,         /**< parent variable */
   SCIP_VAR* arrow_i_j  /**< variable to store */
   )
{
   return put_arrowedge(scip, psd, i, j, arrow_i_j, TRUE);
}

/** Store an edge variable for specified child and parent */
SCIP_RETCODE put_edge(
   SCIP* scip,
   ParentSetData* psd,  /**< the parent sets data structure */
   const int i,         /**< child variable */
   const int j,         /**< parent variable */
   SCIP_VAR* arrow_i_j  /**< variable to store */
   )
{
   return put_arrowedge(scip, psd, i, j, arrow_i_j, FALSE);
}

char* SCIPstrdup(SCIP* scip, const char* str)
{
   char* answer;
   SCIPallocMemoryArray(SCIP, &answer, strlen(str) + 1);
   strcpy(answer, str);
   return answer;
}


/* Some convenient wrappers for creating new parameters that set many values to sensible defaults */
/** Adds a boolean parameter to those recognised by SCIP.
 *
 *  This is just a shortcut for SCIPaddBoolParam() with various options set to their most common values.
 *  Use the full function if you need any of the more advanced options.
 *
 *  @param scip The SCIP instance to add the parameter to.
 *  @param name The parameter's name.
 *  @param desc A description of the parameter.
 *  @param value The parameter's default value.
 *  @return SCIP_OKAY if the operation suceeded.  Otherwise an appropriate error message.
 */
SCIP_RETCODE UT_addBoolParam(
   SCIP* scip,
   const char* name,
   const char* desc,
   SCIP_Bool value
   )
{
   return SCIPaddBoolParam(scip, name, desc, NULL, FALSE, value, NULL, NULL);
}

/** Adds an integer parameter to those recognised by SCIP.
 *
 *  This is just a shortcut for SCIPaddIntParam() with various options set to their most common values.
 *  Use the full function if you need any of the more advanced options. 
 *
 *  @param scip The SCIP instance to add the parameter to.
 *  @param name The parameter's name.
 *  @param desc A description of the parameter.
 *  @param value The parameter's default value.
 *  @param min The parameter's minimum value.
 *  @param max The parameter's maximum value.
 *  @return SCIP_OKAY if the operation suceeded.  Otherwise an appropriate error message.
 */
SCIP_RETCODE UT_addIntParam(
   SCIP* scip,
   const char* name,
   const char* desc,
   int value,
   int min,
   int max
   )
{
   return SCIPaddIntParam(scip, name, desc, NULL, FALSE, value, min, max, NULL, NULL);
}

/** Adds a long integer parameter to those recognised by SCIP.
 *
 *  This is just a shortcut for SCIPaddLongintParam() with various options set to their most common values.
 *  Use the full function if you need any of the more advanced options.  
 *
 *  @param scip The SCIP instance to add the parameter to.
 *  @param name The parameter's name.
 *  @param desc A description of the parameter.
 *  @param value The parameter's default value.
 *  @param min The parameter's minimum value.
 *  @param max The parameter's maximum value.
 *  @return SCIP_OKAY if the operation suceeded.  Otherwise an appropriate error message.
 */
SCIP_RETCODE UT_addLongintParam(
   SCIP* scip,
   const char* name,
   const char* desc,
   SCIP_Longint value,
   SCIP_Longint min,
   SCIP_Longint max
   )
{
   return SCIPaddLongintParam(scip, name, desc, NULL, FALSE, value, min, max, NULL, NULL);
}

/** Adds a real valued parameter to those recognised by SCIP.
 *
 *  This is just a shortcut for SCIPaddRealParam() with various options set to their most common values.
 *  Use the full function if you need any of the more advanced options.
 *
 *  @param scip The SCIP instance to add the parameter to.
 *  @param name The parameter's name.
 *  @param desc A description of the parameter.
 *  @param min The parameter's minimum value.
 *  @param max The parameter's maximum value.
 *  @param value The parameter's default value.
 *  @return SCIP_OKAY if the operation suceeded.  Otherwise an appropriate error message.
 */
SCIP_RETCODE UT_addRealParam(
   SCIP* scip,
   const char* name,
   const char* desc,
   SCIP_Real value,
   SCIP_Real min,
   SCIP_Real max
   )
{
   return SCIPaddRealParam(scip, name, desc, NULL, FALSE, value, min, max, NULL, NULL);
}

/** Adds a string parameter to those recognised by SCIP.
 *
 *  This is just a shortcut for SCIPaddStringParam() with various options set to their most common values.
 *  Use the full function if you need any of the more advanced options.
 *
 *  @param scip The SCIP instance to add the parameter to.
 *  @param name The parameter's name.
 *  @param desc A description of the parameter.
 *  @param value The parameter's default value.
 *  @return SCIP_OKAY if the operation suceeded.  Otherwise an appropriate error message.
 */
SCIP_RETCODE UT_addStringParam(
   SCIP* scip,
   const char* name,
   const char* desc,
   const char* value
   )
{
   return SCIPaddStringParam(scip, name, desc, NULL, FALSE, value, NULL, NULL);
}

/* Some convenient wrappers for creating empty linear constraints that set many values to sensible defaults */
/** Creates an initially empty linear constraint with most options set to sensible defaults.
 *
 *  @param scip SCIP data structure
 *  @param cons pointer to hold the created constraint
 *  @param name name of constraint
 *  @param lhs left hand side of constraint
 *  @param rhs right hand side of constraint
 *  @return SCIP_OKAY if the operation succeeded or an error code otherwise.
 */
SCIP_RETCODE UT_createEmptyLinearConstraint(
   SCIP* scip,
   SCIP_CONS** cons,
   const char* name,
   SCIP_Real lhs,
   SCIP_Real rhs
   )
{
   SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, 0, NULL, NULL, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   return SCIP_OKAY;
}

/** Creates an initially empty less than or equal to linear constraint with most options set to sensible defaults.
 *
 *  @param scip SCIP data structure
 *  @param cons pointer to hold the created constraint
 *  @param name name of constraint
 *  @param rhs The value the constraint should be less than or equal to
 *  @return SCIP_OKAY if the operation succeeded or an error code otherwise.
 */
SCIP_RETCODE UT_createEmptyLTEConstraint(
   SCIP* scip,
   SCIP_CONS** cons,
   const char* name,
   SCIP_Real rhs
   )
{
   SCIP_CALL( UT_createEmptyLinearConstraint(scip, cons, name, -SCIPinfinity(scip), rhs) );
   return SCIP_OKAY;
}
/** Creates an initially empty greater than or equal to linear constraint with most options set to sensible defaults.
 *
 *  @param scip SCIP data structure
 *  @param cons pointer to hold the created constraint
 *  @param name name of constraint
 *  @param lhs The value the constraint should be greater than or equal to
 *  @return SCIP_OKAY if the operation succeeded or an error code otherwise.
 */
SCIP_RETCODE UT_createEmptyGTEConstraint(
   SCIP* scip,
   SCIP_CONS** cons,
   const char* name,
   SCIP_Real lhs
   )
{
   SCIP_CALL( UT_createEmptyLinearConstraint(scip, cons, name, lhs, SCIPinfinity(scip)) );
   return SCIP_OKAY;
}

/* Functions for writing and parsing items */

/* Escaping functions */
/** The character used to escape special characters when writing to a file. */
#define ESC_CHAR '\\'
/** Whether a character needs an escape character adding before it can be written to fiel.
 *
 *  @param item The character to consider
 *  @return True if the character needs to be escaped.  False otherwise.
 */
static
SCIP_Bool needsEscaping(
   char item
   )
{
   switch(item)
   {
   case UT_LIST_START:
   case UT_LIST_END:
   case UT_LIST_SEP:
   case ESC_CHAR:
   case ' ':
      return TRUE;
   default:
      return FALSE;
   }
}
/** Adds an necessary escape characters to a string so it will be re-read correctly.
 *
 *  @param input_string The initial string.
 *  @param output_string A copy of the string after any escape characters have been added.
 *  @return SCIP_OKAY if the operation succeeded.
 */
static
SCIP_RETCODE escapeString(
   char* input_string,
   char** output_string
   )
{
   int input_pos;
   int output_pos;
   for( input_pos = 0, output_pos = 0; input_string[input_pos] != '\0'; input_pos++, output_pos++ )
   {
      if( needsEscaping(input_string[input_pos]) )
      {
         (*output_string)[output_pos++] = ESC_CHAR;
      }
      (*output_string)[output_pos] = input_string[input_pos];
   }
   (*output_string)[output_pos] = '\0';
   return SCIP_OKAY;
}

/** Removes an escape charcters from a string to return the string before it was escaped.
 *
 *  @param input_string The initial string.
 *  @param output_string A copy of the string after any escape characters have been removed.
 *  @return SCIP_OKAY if the operation succeeded.
 */
static
SCIP_RETCODE unescapeString(
   char* input_string,
   char** output_string
   )
{
   int input_pos;
   int output_pos = 0;
   for( input_pos = 0; input_string[input_pos] != '\0'; input_pos++ )
      if( input_string[input_pos] != ESC_CHAR )
         (*output_string)[output_pos++] = input_string[input_pos];
   (*output_string)[output_pos] = '\0';
   return SCIP_OKAY;
}

/* Writing functions */
/** Writes an array of strings to a file.
 *
 *  @param scip The SCIP instance that the data belongs to.
 *  @param file The file to write to.
 *  @param input_array The string array to write to file.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_writeStringArray(
   SCIP* scip,
   FILE* file,
   char** input_array,
   int length
   )
{
   int i;
   char* next_item;
   SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
   for( i = 0; i < length; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &next_item, 2*strlen(input_array[i])+1) );
      if( i > 0 )
         SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIP_CALL( escapeString(input_array[i], &next_item) );
      SCIPinfoMessage(scip, file, "%s", next_item);
      SCIPfreeMemoryArray(scip, &next_item);
   }
   SCIPinfoMessage(scip, file, "%c", UT_LIST_END);

   return SCIP_OKAY;
}

/** Writes an array of integers to a file.
 *
 *  @param scip The SCIP instance that the data belongs to.
 *  @param file The file to write to.
 *  @param input_array The integer array to write to file.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_writeIntArray(
   SCIP* scip,
   FILE* file,
   int* input_array,
   int length
   )
{
   int i;
   char** tmp_array;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array), length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), SCIP_MAXSTRLEN) );

   for( i = 0; i < length; i++ )
      sprintf(tmp_array[i], "%d", input_array[i]);
   SCIP_CALL( UT_writeStringArray(scip, file, tmp_array, length) );

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &(tmp_array));

   return SCIP_OKAY;
}

/** Writes an array of characters to a file.
 *
 *  @param scip The SCIP instance that the data belongs to.
 *  @param file The file to write to.
 *  @param input_array The character array to write to file.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_writeCharArray(
   SCIP* scip,
   FILE* file,
   char* input_array,
   int length
   )
{
   int i;
   char** tmp_array;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array), length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), SCIP_MAXSTRLEN) );

   for( i = 0; i < length; i++ )
      sprintf(tmp_array[i], "%c", input_array[i]);
   SCIP_CALL( UT_writeStringArray(scip, file, tmp_array, length) );

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &(tmp_array));

   return SCIP_OKAY;
}

/** Writes the names of an array of variables to a file.
 *
 *  @param scip The SCIP instance that the data belongs to.
 *  @param file The file to write to.
 *  @param input_array The variable array to write to file.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_writeVarArray(
   SCIP* scip,
   FILE* file,
   SCIP_VAR** input_array,
   int length
   )
{
   int i;
   char** tmp_array;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array), length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), SCIP_MAXSTRLEN) );

   for( i = 0; i < length; i++ )
      if( input_array[i] != NULL )
         sprintf(tmp_array[i], "%s", SCIPvarGetName(input_array[i]));
      else
         sprintf(tmp_array[i], "%s", "NULL");
   SCIP_CALL( UT_writeStringArray(scip, file, tmp_array, length) );

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &(tmp_array));
   return SCIP_OKAY;
}

/** Writes an array of arrays of integers to a file.
 *
 *  @param scip The SCIP instance that the data belongs to.
 *  @param file The file to write to.
 *  @param input_array The integer array array to write to file.
 *  @param lengths The number of items in each sub array.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_writeIntArrayArray(
   SCIP* scip,
   FILE* file,
   int** input_array,
   int* lengths,
   int length
   )
{
   int i;

   SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
   for( i = 0; i < length; i++ )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIP_CALL( UT_writeIntArray(scip, file, input_array[i], lengths[i]) );
   }
   SCIPinfoMessage(scip, file, "%c", UT_LIST_END);
   return SCIP_OKAY;
}

/** Writes an array of arrays of variables to a file.
 *
 *  @param scip The SCIP instance that the data belongs to.
 *  @param file The file to write to.
 *  @param input_array The variable array array to write to file.
 *  @param lengths The number of items in each sub array.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_writeVarArrayArray(
   SCIP* scip,
   FILE* file,
   SCIP_VAR*** input_array,
   int* lengths,
   int length
   )
{
   int i;

   SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
   for( i = 0; i < length; i++ )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIP_CALL( UT_writeVarArray(scip, file, input_array[i], lengths[i]) );
   }
   SCIPinfoMessage(scip, file, "%c", UT_LIST_END);
   return SCIP_OKAY;
}

/** Writes an array of arrays of arrays of integers to a file.
 *
 *  @param scip The SCIP instance that the data belongs to.
 *  @param file The file to write to.
 *  @param input_array The integer array array array to write to file.
 *  @param lengths The number of items in each sub sub array.
 *  @param lengths_of_lengths The number of items in each sub array.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_writeIntArrayArrayArray(
   SCIP* scip,
   FILE* file,
   int*** input_array,
   int** lengths,
   int* lengths_of_lengths,
   int length
   )
{
   int i;

   SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
   for( i = 0; i < length; i++ )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIP_CALL( UT_writeIntArrayArray(scip, file, input_array[i], lengths[i], lengths_of_lengths[i]) );
   }
   SCIPinfoMessage(scip, file, "%c", UT_LIST_END);
   return SCIP_OKAY;
}

/* Parsing functions */
/** Parses an array from a string.
 *
 *  Use this function only if the array contains multiple types of items, or items that
 *  none of the other `UT_parseXXXArray()` functions deals with.
 *
 *  @param input_string The string representation of the array to parse.
 *  @param output_array The array represented by the string.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_parseArray(
   char* input_string,
   char*** output_array
   )
{
   int i;
   int pos = 0;
   int begin_pos = 0;
   int array_pos = 0;
   int sublist_depth = 0;
   int state = 0;

   while( input_string[pos] != '\0' )
   {
      /*printf("%d[%d] - %d %c %s\n", state, sublist_depth, pos, input_string[pos], input_string); */
      switch(state)
      {
      case 0:
         /* At beginning of string */
         switch(input_string[pos])
         {
         case ' ':
            /* Skip white space here */
            pos++;
            begin_pos = pos;
            break;
         case UT_LIST_START:
            pos++;
            begin_pos = pos;
            state = 1;
            break;
         default:
            SCIPerrorMessage("Unable to UT_parse \"%s\".  Does not begin with %c.\n", input_string, UT_LIST_START);
            return SCIP_ERROR;
            break;
         }
         break;
      case 1:
         /* Expecting beginning of new item or end of list */
         switch(input_string[pos])
         {
         case ' ':
            /* Skip white space here */
            pos++;
            begin_pos = pos;
            break;
         case UT_LIST_START:
            /* Entering new item which is a sub list */
            pos++;
            state = 2;
            sublist_depth = 1;
            break;
         case UT_LIST_END:
            /* Reached end of list (hopefully) */
            pos++;
            begin_pos = pos;
            state = 3;
            break;
         case UT_LIST_SEP:
            /* Shouldn't happen */
            SCIPerrorMessage("Unable to UT_parse \"%s\".  Error at position %d.\n", input_string, pos);
            return SCIP_ERROR;
            break;
         case ESC_CHAR:
            /* Entering new item that happens to start with ESC_CHAR */
            pos += 2;
            state = 4;
            break;
         default:
            /* Entering new item */
            pos++;
            state = 4;
            break;
         }
         break;
      case 2:
         /* In the midst of an item which is a sublist */
         switch(input_string[pos])
         {
         case UT_LIST_START:
            /* Another sublist */
            pos++;
            sublist_depth += 1;
            break;
         case UT_LIST_END:
            /* End of a sublist */
            sublist_depth -= 1;
            if( sublist_depth == 0 )
            {
               pos++;
               /* Hit the end of the item */
               for( i = 0; i < pos - begin_pos; i++ )
               {
                  (*output_array)[array_pos][i] = input_string[i + begin_pos];
               }
               (*output_array)[array_pos][i] = '\0';
               begin_pos = pos;
               array_pos++;
               state = 5;
            }
            else
               pos++;
            break;
         case ESC_CHAR:
            pos += 2;
            break;
         default:
            pos++;
            break;
         }
         break;
      case 3:
         /* Expecting to be at the end of the list */
         switch(input_string[pos])
         {
         case ' ':
            /* Just skip over any white space */
            begin_pos = pos;
            pos++;
            break;
         case UT_LIST_START:
         case UT_LIST_END:
         case UT_LIST_SEP:
         case ESC_CHAR:
         default:
            /* Shouldn't happen */
            SCIPerrorMessage("Unable to UT_parse \"%s\".  Unexpected %c at position %d.\n", input_string, input_string[pos], pos);
            return SCIP_ERROR;
            break;
         }
         break;
      case 4:
         /* In the midst of an item */
         switch(input_string[pos])
         {
         case ' ':
            /* Hit the end of the item */
            for( i = 0; i < pos - begin_pos; i++ )
               (*output_array)[array_pos][i] = input_string[i + begin_pos];
            (*output_array)[array_pos][i] = '\0';
            array_pos++;
            pos++;
            begin_pos = pos;
            state = 5;
            break;
         case UT_LIST_START:
            /* Shouldn't happen */
            SCIPerrorMessage("Unable to UT_parse \"%s\".  Unexpected %c at position %d.\n", input_string, UT_LIST_START, pos);
            return SCIP_ERROR;
            break;
         case UT_LIST_END:
            /* End of the item, and of the list */
            for( i = 0; i < pos - begin_pos; i++ )
               (*output_array)[array_pos][i] = input_string[i + begin_pos];
            (*output_array)[array_pos][i] = '\0';
            array_pos++;
            pos++;
            begin_pos = pos;
            state = 3;
            break;
         case UT_LIST_SEP:
            /* Hit the end of the item */
            for( i = 0; i < pos - begin_pos; i++ )
               (*output_array)[array_pos][i] = input_string[i + begin_pos];
            (*output_array)[array_pos][i] = '\0';
            array_pos++;
            pos++;
            begin_pos = pos;
            state = 6;
            break;
         case ESC_CHAR:
            /* Escaped char in item */
            pos += 2;
            break;
         default:
            /* Normal char in item */
            pos ++;
            break;
         }
         break;
      case 5:
         /* Expecting either a UT_LIST_SEP or a UT_LIST_END */
         switch(input_string[pos])
         {
         case ' ':
            /* Skip white space here */
            pos++;
            begin_pos = pos;
            break;
         case UT_LIST_END:
            pos++;
            begin_pos = pos;
            state = 3;
            break;
         case UT_LIST_SEP:
            pos++;
            begin_pos = pos;
            state = 6;
            break;
         case UT_LIST_START:
         case ESC_CHAR:
         default:
            SCIPerrorMessage("Unable to UT_parse \"%s\".  Unexpected %c at position %d.\n", input_string, input_string[pos], pos);
            return SCIP_ERROR;
            break;
         }
         break;
      case 6:
         /* Expecting beginning of new item */
         switch(input_string[pos])
         {
         case ' ':
            /* Skip white space here */
            pos++;
            begin_pos = pos;
            break;
         case UT_LIST_START:
            /* Entering new item which is a sub list */
            pos++;
            state = 2;
            sublist_depth = 1;
            break;
         case UT_LIST_END:
         case UT_LIST_SEP:
            /* Shouldn't happen */
            SCIPerrorMessage("Unable to UT_parse \"%s\".  Error at position %d.\n", input_string, pos);
            return SCIP_ERROR;
            break;
         case ESC_CHAR:
            /* Entering new item that happens to start with ESC_CHAR */
            pos += 2;
            state = 4;
            break;
         default:
            /* Entering new item */
            pos++;
            state = 4;
            break;
         }
         break;
      }
   }
   return SCIP_OKAY;
}
/** Parses a string array from a string.
 *
 *  @param input_string The string representation of the array to parse.
 *  @param output_array The array represented by the string.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_parseStringArray(
   char* input_string,
   char*** output_array,
   int length
   )
{
   int i;
   char** tmp_array;
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), strlen(input_string)) );

   SCIP_CALL( UT_parseArray(input_string, &tmp_array) );
   for( i = 0; i < length; i++ )
      unescapeString(tmp_array[i], &((*output_array)[i]));

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &tmp_array);
   return SCIP_OKAY;
}

/** Parses an integer array from a string.
 *
 *  @param input_string The string representation of the array to parse.
 *  @param output_array The array represented by the string.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_parseIntArray(
   char* input_string,
   int** output_array,
   int length
   )
{
   int i;
   char** tmp_array;
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), strlen(input_string)) );

   SCIP_CALL( UT_parseStringArray(input_string, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      sscanf(tmp_array[i], "%d", &((*output_array)[i]));

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &tmp_array);
   return SCIP_OKAY;
}
/** Parses a character array from a string.
 *
 *  @param input_string The string representation of the array to parse.
 *  @param output_array The array represented by the string.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_parseCharArray(
   char* input_string,
   char** output_array,
   int length
   )
{
   int i;
   char** tmp_array;

   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), strlen(input_string)) );

   SCIP_CALL( UT_parseStringArray(input_string, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      sscanf(tmp_array[i], "%c", &((*output_array)[i]));

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &tmp_array);
   return SCIP_OKAY;
}

/** Parses a variable array from a string.
 *
 *  The variables are found from their names in the input.  It is therefore important
 *  that the SCIP instance knows the variables so that they can be found from these
 *  names.
 *
 *  @param scip The SCIP instance in which the variables appear.
 *  @param input_string The string representation of the array to parse.
 *  @param output_array The array represented by the string.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_parseVarArray(
   SCIP* scip,
   char* input_string,
   SCIP_VAR*** output_array,
   int length
   )
{
   int i;
   char** tmp_array;
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), strlen(input_string)) );

   SCIP_CALL( UT_parseStringArray(input_string, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      (*output_array)[i] = SCIPfindVar(scip, tmp_array[i]);

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &tmp_array);
   return SCIP_OKAY;
}

/** Parses an array of integer arrays from a string.
 *
 *  @param input_string The string representation of the array to parse.
 *  @param output_array The array represented by the string.
 *  @param lengths The number of items in each sub array.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_parseIntArrayArray(
   char* input_string,
   int*** output_array,
   int* lengths,
   int length
   )
{
   int i;
   char** tmp_array;
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, length) );

   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), strlen(input_string)) );

   SCIP_CALL( UT_parseArray(input_string, &tmp_array) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( UT_parseIntArray(tmp_array[i], &((*output_array)[i]), lengths[i]) );

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &tmp_array);
   return SCIP_OKAY;
}

/** Parses an array of variable arrays from a string.
 *
 *  The variables are found from their names in the input.  It is therefore important
 *  that the SCIP instance knows the variables so that they can be found from these
 *  names.
 *
 *  @param scip The SCIP instance in which the variables appear.
 *  @param input_string The string representation of the array to parse.
 *  @param output_array The array represented by the string.
 *  @param lengths The number of items in each sub array.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_parseVarArrayArray(
   SCIP* scip,
   char* input_string,
   SCIP_VAR**** output_array,
   int* lengths,
   int length
   )
{
   int i;
   char** tmp_array;
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), strlen(input_string)) );

   SCIP_CALL( UT_parseArray(input_string, &tmp_array) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( UT_parseVarArray(scip, tmp_array[i], &((*output_array)[i]), lengths[i]) );

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &tmp_array);
   return SCIP_OKAY;
}

/** Parses an array of arrays of integer arrays from a string.
 *
 *  @param input_string The string representation of the array to parse.
 *  @param output_array The array represented by the string.
 *  @param lengths The number of items in each sub sub array.
 *  @param lengths_of_lengths The number of items in each sub array.
 *  @param length The number of items in the array.
 *  @return SCIP_OKAY if the operation succeeded.
 */
SCIP_RETCODE UT_parseIntArrayArrayArray(
   char* input_string,
   int**** output_array,
   int** lengths,
   int* lengths_of_lengths,
   int length
   )
{
   int i;
   char** tmp_array;
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_array, length) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_array[i]), strlen(input_string)) );

   SCIP_CALL( UT_parseArray(input_string, &tmp_array) );
   for( i = 0; i < length; i++ )
      SCIP_CALL( UT_parseIntArrayArray(tmp_array[i], &((*output_array)[i]), lengths[i], lengths_of_lengths[i]) );

   for( i = 0; i < length; i++ )
      SCIPfreeMemoryArray(scip, &(tmp_array[i]));
   SCIPfreeMemoryArray(scip, &tmp_array);
   return SCIP_OKAY;
}

/* Functions for reading from a file */
static
SCIP_Bool isSplitter(
   char item,
   char* splitters,
   int num_splitters
   )
{
   int s;

   for( s = 0; s < num_splitters; s++ )
      if( splitters[s] == item )
         return TRUE;
   return FALSE;
}

static
int convertToSimpleForm(
   char* line,
   char* splitters,
   int num_splitters,
   SCIP_Bool multiple_splitters_as_one,
   char** output_line
   )
{
   char basic_splitter = splitters[0];
   int current_read_pos = 0;
   int current_write_pos = 0;
   SCIP_Bool last_char_was_splitter = FALSE;

   for( current_read_pos = 0; !(line[current_read_pos] == '\0' || line[current_read_pos] == '\n' || line[current_read_pos] == '\r'); current_read_pos++ )
      if( isSplitter(line[current_read_pos], splitters, num_splitters) )
      {
         if( !last_char_was_splitter || !multiple_splitters_as_one )
         {
            (*output_line)[current_write_pos] = basic_splitter;
            current_write_pos++;
         }
         last_char_was_splitter = TRUE;
      }
      else
      {
         last_char_was_splitter = FALSE;
         (*output_line)[current_write_pos] = line[current_read_pos];
         current_write_pos++;
      }
   (*output_line)[current_write_pos] = '\0';
   return current_write_pos;
}

static
int countNumberOfSegments(
   char* line,
   int length,
   char splitter
   )
{
   int num = 1;
   int i;
   for( i = 0; i < length; i++ )
      if( line[i] == splitter )
         num += 1;
   return num;
}

static
SCIP_RETCODE readNextLine(
   SCIP* scip,
   FILE* file,
   int* length,
   char** buffer,
   char** status
   )
{
   int initial_length = (*length);

   (*status) = fgets(*buffer, *length, file);

   /* Repeat until we have a new line in the buffer (or at EOF, to cope with missing newlines) */
   while( strchr(*buffer, '\n') == NULL && !feof(file) )
   {
      int pos = (*length);

      /* Reallocate memory */
      (*length) = SCIPcalcMemGrowSize(scip, (*length) + initial_length);
      SCIP_CALL( SCIPreallocBufferArray(scip, buffer, (*length)) );

      /* read next line */
      (*status) = fgets(&(*buffer[pos]), initial_length, file);
   }

   return SCIP_OKAY;
}

/** Reads a file and splits it into lines, which are themselves split into white space delimited tokens.
 *
 *  @param scip The SCIP instance being used.
 *  @param file The file to read.
 *  @param lines A pointer to the lines of text that are read from the file.
 *  @param splitters The field delimiters to use.
 *  @param num_splitters The number of field delimiters to use.
 *  @param multiple_splitters_as_one Whether multiple field delimiters should be merged in to one.
 *  @param num_lines The number of lines of text read.
 *  @param line_lengths The number of tokens in each line.
 *  @return SCIP_OKAY if the file was correctly read or a NULL pointer otherwise.
 */
SCIP_RETCODE UT_readFileAndSplit(
   SCIP* scip,
   FILE* file,
   char* splitters,
   int num_splitters,
   SCIP_Bool multiple_splitters_as_one,
   char**** lines,
   int* num_lines,
   int** line_lengths
   )
{
   char*** tmp_lines = NULL;
   int* tmp_line_lengths = NULL;
   int i, j;
   char* line_status;
   char*  line;
   int   line_length = 131071;
   int   current_line_number = 0;
   int   num_lines_increment = 10000;

   /* For each line in the file... */
   SCIP_CALL( SCIPallocBufferArray(scip, &line, line_length) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_lines, num_lines_increment) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_line_lengths, num_lines_increment) );
   SCIP_CALL( readNextLine(scip, file, &line_length, &line, &line_status) );
   while( line_status != NULL )
   {
      char* basic_line;
      int num_chars;
      int num_segments;
      int item_start_pos = -1;
      int current_pos = 0;
      int current_item_number = 0;

      if( (current_line_number  > 0) && (current_line_number % num_lines_increment == 0) )
      {
         int new_length = SCIPcalcMemGrowSize(scip, current_line_number + num_lines_increment);
         SCIP_CALL( SCIPreallocMemoryArray(scip, &tmp_lines, new_length) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &tmp_line_lengths, new_length) );
      }

      SCIP_CALL( SCIPallocMemoryArray(scip, &basic_line, line_length) );
      num_chars = convertToSimpleForm(line, splitters, num_splitters, multiple_splitters_as_one, &basic_line);
      num_segments = countNumberOfSegments(basic_line, num_chars, splitters[0]);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_lines[current_line_number]), num_segments) );

      /* printf("\"%s\"\n\"%s\"\n", line, basic_line); */

      for( current_pos = 0; current_pos < num_chars; current_pos++ )
      {
         if( basic_line[current_pos] == splitters[0] )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_lines[current_line_number][current_item_number]), current_pos - item_start_pos) );
            for( i = item_start_pos + 1; i < current_pos; i++ )
               tmp_lines[current_line_number][current_item_number][i - item_start_pos - 1] = basic_line[i];
            tmp_lines[current_line_number][current_item_number][current_pos - item_start_pos - 1] = '\0';
            current_item_number++;
            tmp_line_lengths[current_line_number] = current_item_number;
            item_start_pos = current_pos;
         }
      }

      if( item_start_pos != current_pos )
      {
         /* Deal with anything after the final separator */
         SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_lines[current_line_number][current_item_number]), current_pos - item_start_pos) );
         for( i = item_start_pos + 1; i < current_pos; i++ )
            tmp_lines[current_line_number][current_item_number][i - item_start_pos - 1] = basic_line[i];
         tmp_lines[current_line_number][current_item_number][current_pos - item_start_pos - 1] = '\0';
         current_item_number++;
         tmp_line_lengths[current_line_number] = current_item_number;
      }

      SCIPfreeMemoryArray(scip, &basic_line);

      current_line_number++;
      SCIP_CALL( readNextLine(scip, file, &line_length, &line, &line_status) );
   }

   /* Copy it to the real data structures */
   (*num_lines) = current_line_number;
   SCIP_CALL( SCIPallocMemoryArray(scip, lines, current_line_number) );
   SCIP_CALL( SCIPallocMemoryArray(scip, line_lengths, current_line_number) );
   for( i = 0; i < current_line_number; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*lines)[i]), tmp_line_lengths[i]) );
      (*line_lengths)[i] = tmp_line_lengths[i];
      for( j = 0; j < tmp_line_lengths[i]; j++ )
         (*lines)[i][j] = tmp_lines[i][j];
      SCIPfreeMemoryArray(scip, &(tmp_lines[i]));
   }
   SCIPfreeMemoryArray(scip, &tmp_lines);
   SCIPfreeMemoryArray(scip, &tmp_line_lengths);
   SCIPfreeBufferArray(scip, &line);

   /*
   for (i = 0; i < (*num_lines); i++) {
      for (j = 0; j < (*line_lengths)[i]; j++)
         printf("\"%s\" ", (*lines)[i][j]);
      printf("\n");
   }
   */


   return SCIP_OKAY;
}

/**< either check that a DAG is the distinguished representative of its Markov equivalence class 
 *  or construct the distinguished representative
 *  @return If just checking, return whether the input DAG is the distinguished representative of its Markov equivalence class. If constructing, return whether it is possible to construct the  distinguished representative. 
 */
SCIP_Bool is_dr_feasible(
   SCIP* scip,             /**< SCIP data structure */
   const ParentSetData* psd,     /**< Parent set data structure */
   SCIP_SOL* sol,          /**< DAG to check/use as input */
   SCIP_Bool constructing, /**< whether we are constructing the distinguished representative (if any) */
   SCIP_SOL* dr_sol        /**< dr_sol is the constructed distinguished representative (if any) */
   )
{

   int i;
   int j;
   int k;
   int kk;
   int l;
   int h;
   
   int n = psd->n;

   int n_todo;
   int old_n_todo;

   SCIP_Bool* todo;
   int** ch;
   int* n_ch;
   int* pa_idx;
   int* children;

   int n_children;
   int child;
   int childpa;

   SCIP_VAR* var;

   SCIP_Bool ok;
   SCIP_Bool all_found;
   SCIP_Bool result;

   int* newparents;
   int n_newparents;

   SCIP_Real val;
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &todo, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &ch, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &n_ch, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &children, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pa_idx, n) );
   if( constructing )
      SCIP_CALL( SCIPallocMemoryArray(scip, &newparents, n) );
   for( i = 0; i < n; ++i )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(ch[i]), n) );

   /* collect children and parents in this solution for each node */
   SCIPdebugMessage("Checking this solution\n");
   for( i = 0; i < n; ++i )
   {
      n_ch[i] = 0;
      for( j = 0; j < n; ++j )
      {
         var = get_arrow(psd,j,i);
         if( var == NULL )
            continue;
         
         val = SCIPgetSolVal(scip, sol, var);
         /* non-integral solutions are not distinguished representatives */
         if( !SCIPisIntegral(scip, val) )
         {
            result = FALSE;
            goto ABORT;
         }
         if( SCIPisGT(scip, val, 0.5) )
            ch[i][n_ch[i]++] = j;
      }
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(ch[i]), n_ch[i]) );

      pa_idx[i] = -1;
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         if( SCIPisGT(scip, SCIPgetSolVal(scip, sol, psd->PaVars[i][k]), 0.5) )
         {
            SCIPdebugMessage("%s \n", SCIPvarGetName(psd->PaVars[i][k]));
            pa_idx[i] =  k;
            break;
         }
      }
      
      if( pa_idx[i] == -1 )
      {
         /* somehow couldn't find a set family variable */
         result = FALSE;
         goto ABORT;
      }
   }

   /* all nodes yet to be deleted */
   for( i = 0; i < n; ++i )
      todo[i] = TRUE;

   n_todo = n;
   result = TRUE;
   while( result == TRUE && n_todo > 0 )
   {
      old_n_todo = n_todo;
      for( i = n-1; i > -1; --i )
      {
         if( !todo[i] )
            continue;

         /* store index of current parent set for i */
         k = pa_idx[i];

         /* collect undeleted children */
         n_children = 0;
         for( j = 0; j < n_ch[i]; ++j )
         {
            child = ch[i][j];
            if( todo[child] )
               children[n_children++] = child;
         }

         if( n_children == 0 && !constructing )
         {
            /* i is a sink in subgraph of undeleted vertices
               so clearly can be 'made a sink' without altering MEC! 
               If not constructing, just delete it and continue */
            SCIPdebugMessage("Sink %d\n", i);
            todo[i] = FALSE;
            n_todo--;
            break;
         }

         /* i is not a sink, so unless there is something to stop it 'becoming one'
            this is not ok
         */
         ok = FALSE;

         for( j = 0; j < n_children; ++j )
         {
            child = children[j];
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               var = get_edge(psd,psd->ParentSets[i][k][l],child);
               if( var == NULL || SCIPisLT(scip, SCIPgetSolVal(scip, sol, var), 0.5) )
               {
                  /* reversing i->child would create child->i<-(psd->ParentSets[i][k][l]) immorality */
                  ok = TRUE;
                  break;
               }
            }
            
            if( ok )
               /* not a sink, but cannot be made one without creating an immorality
                  no need to inspect further children */
               break;
            
            for( l = 0; l < psd->nParents[child][pa_idx[child]]; ++l )
            {
               childpa = psd->ParentSets[child][pa_idx[child]][l];
               if( childpa == i )
                  continue;
               var = get_edge(psd,childpa,i);
               if( var == NULL || SCIPisLT(scip, SCIPgetSolVal(scip, sol, var), 0.5) )
               {
                  /* reversing i->child would remove i->child<-childpa immorality */
                  ok = TRUE;
                  break;
               }
            }

            if( ok )
               /* not a sink, but cannot be made one without removing an immorality
                  no need to inspect further children */
               break;

            for( l = j+1; l < n_children; ++l )
            {
               var = get_edge(psd,child,children[l]);
               if( var == NULL || SCIPisLT(scip, SCIPgetSolVal(scip, sol, var), 0.5) )
               {
                  /* reversing i->child and i->children[l] would create an immorality */
                  ok = TRUE;
                  break;
               }
            }

            if( ok )
               /* not a sink, but cannot be made one without creating an immorality
                  no need to inspect further children */

               break;
         }
         
         if( !ok )
         {
            /* This DAG is not a distinguished rep
               since i is not a sink, but could be made into one without alterning the MEC
            */

            if( constructing )
            {
               /* in the dr_sol make i a sink in the subgraph of remaining vertices */
               /* new parent set is all undeleted original parents plus all undeleted children */
               n_newparents = 0;
               for( j = 0; j < n_children; ++j )
                  newparents[n_newparents++] = children[j];
               for( l = 0; l < psd->nParents[i][k]; ++l )
                  if( todo[psd->ParentSets[i][k][l]] )
                     newparents[n_newparents++] = psd->ParentSets[i][k][l];
                     
               SCIPsortInt(newparents,n_newparents);

               /* find correct new parent set (if any) */
               for( kk = 0; kk < psd->nParentSets[i]; ++kk )
               {
                  if( psd->nParents[i][kk] != n_newparents )
                     continue;

                  all_found = TRUE; 
                  for( l = 0; l < n_newparents; ++l )
                     if( newparents[l] != psd->ParentSets[i][kk][l] )
                     {
                        all_found = FALSE;
                        break;
                     }

                  if( all_found )
                  {
                     SCIP_CALL( SCIPsetSolVal(scip, dr_sol, psd->PaVars[i][kk], 1.0) );
                     /* printf("Setting: %s\n", SCIPvarGetName(psd->PaVars[i][kk]));  */
                     for( l = 0; l < psd->nParents[i][kk]; l++ )
                     {
                        /* set both arrow and edge indicator variables, appropriately */
                        for(h = 0; h < 2; h++)
                        {
                           var = get_arrowedge(psd, i, psd->ParentSets[i][kk][l], (h==0));
                           assert( var != NULL );
                           if( SCIPvarIsActive(var) &&
                              ( SCIPvarGetStatus(SCIPvarGetTransVar(var)) != SCIP_VARSTATUS_MULTAGGR
                                 || SCIPvarGetMultaggrNVars(SCIPvarGetTransVar(var)) == 1 ) )
                              SCIP_CALL( SCIPsetSolVal(scip, dr_sol, var, 1.0) );
                        }
                     }
                     break;
                  }
               }
               if( kk == psd->nParentSets[i] )
               {
                  /* did not find the required parent set 
                   so give up
                  */
                  /* printf("%d: ",i); */
                  /* for( l = 0; l < n_newparents; ++l ) */
                  /*    printf("%d ",newparents[l]); */
                  /* printf("\n"); */
                  result = FALSE;
               }
               todo[i] = FALSE;
               n_todo--;
            }
            else
            {
               /* if not constructing just record that this is not a distinguished representative */
               result = FALSE;
            }
            /* have now found highest indexed vertex so exit for loop */
            break;
         }
      }
      if( old_n_todo == n_todo )
         /* could not find a sink, must be cyclic */
         result = FALSE;
   }

ABORT:
   for( i = 0; i < n; ++i )
      SCIPfreeMemoryArray(scip, &(ch[i]));

   SCIPfreeMemoryArray(scip, &todo);
   SCIPfreeMemoryArray(scip, &ch);
   SCIPfreeMemoryArray(scip, &n_ch);
   SCIPfreeMemoryArray(scip, &children);
   SCIPfreeMemoryArray(scip, &pa_idx);
   if( constructing )
      SCIPfreeMemoryArray(scip, &newparents);

   if( !constructing )
   {
      if( result )
         SCIPdebugMessage("Solution feasible\n");
      else
         SCIPdebugMessage("Solution infeasible\n");
   }
   return result;
}
      
