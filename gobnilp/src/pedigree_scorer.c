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

#include "pedigree_scorer.h"
#include "pedigree_data.h"
#include "metadata.h"
#include "utils.h"
#include "versiongit.h"
#include <string.h>

/** @file
 *  Contains all the functions related to create scores for pedigrees in the program.
 */

static
SCIP_Bool strequ(
   const char* str1,
   const char* str2
   )
{
   return ( strcmp(str1, str2) == 0 );
}

static
SCIP_RETCODE getColumnIndices(
   SCIP* scip,
   int* name_column,
   int* sex_column,
   int* age_column,
   int* genome_start_column,
   int* genome_end_column,
   int* adult_column
   )
{
   SCIPgetIntParam(scip, "gobnilp/pedigree/inputfile/namecolumn", name_column);
   SCIPgetIntParam(scip, "gobnilp/pedigree/inputfile/sexcolumn", sex_column);
   SCIPgetIntParam(scip, "gobnilp/pedigree/inputfile/agecolumn", age_column);
   SCIPgetIntParam(scip, "gobnilp/pedigree/inputfile/allelestartcolumn", genome_start_column);
   SCIPgetIntParam(scip, "gobnilp/pedigree/inputfile/alleleendcolumn", genome_end_column);
   SCIPgetIntParam(scip, "gobnilp/pedigree/inputfile/adultcolumn", adult_column);
   return SCIP_OKAY;
}

static
int findNumColumnsNeeded(
   int name_column,
   int sex_column,
   int age_column,
   int genome_start_column,
   int genome_end_column,
   int adult_column
   )
{
   return MAX(name_column, MAX(sex_column, MAX(age_column, MAX(genome_start_column, genome_end_column))));
}

static
SCIP_RETCODE openFile(
   const char* filename,
   FILE** file
   )
{
   if( strequ(filename, "-") )
      (*file) = stdin;
   else
      (*file) = fopen(filename, "r");

   if( (*file) == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", filename);
      return SCIP_NOFILE;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE checkLineLengths(
   int num_lines,
   int* line_lengths,
   int min_columns_needed
   )
{
   int i;
   for( i = 1; i < num_lines; i++ )
   {
      if( line_lengths[0] != line_lengths[i] )
      {
         SCIPerrorMessage("Wrong number of data items on line %d.  Found %d when %d were expected.\n", i + 1, line_lengths[i], line_lengths[0]);
         return SCIP_READERROR;
      }
   }

   if( line_lengths[0] <= min_columns_needed )
   {
      SCIPerrorMessage("Settings specify that data is needed from column %d, but only %d columns were found.\n", min_columns_needed, line_lengths[0]);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE readFromFile(
   SCIP* scip,
   const char* filename,
   char* delims,
   int num_delims,
   SCIP_Bool merge_delims,
   char**** lines,
   int* num_lines,
   int** line_lengths,
   int min_columns_needed
   )
{
   FILE* file;
   SCIP_CALL( openFile(filename, &file) );
   SCIP_CALL( UT_readFileAndSplit(scip, file, delims, num_delims, merge_delims, lines, num_lines, line_lengths) );
   fclose(file);
   SCIP_CALL( checkLineLengths((*num_lines), (*line_lengths), min_columns_needed) );
   return SCIP_OKAY;
}

static
SCIP_RETCODE extractNameColumn(
   SCIP* scip,
   char*** names,
   char*** lines,
   int num_lines,
   int first_row,
   int column
   )
{
   int i;
   SCIP_CALL( SCIPallocMemoryArray(scip, names, num_lines - first_row) );
   for( i = first_row; i < num_lines; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*names)[i - first_row]), SCIP_MAXSTRLEN) );
      if( column > -1 )
         sprintf((*names)[i - first_row], "%s", lines[i][column]);
      else
         sprintf((*names)[i - first_row], "%d", i);
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE extractSexColumn(
   SCIP* scip,
   char** sexes,
   char*** lines,
   int num_lines,
   int first_row,
   int column
   )
{
   int i;
   SCIP_CALL( SCIPallocMemoryArray(scip, sexes, num_lines - first_row) );
   for( i = first_row; i < num_lines; i++ )
   {
      if( lines[i][column][0] == 'm' || lines[i][column][0] == 'M' )
         (*sexes)[i - first_row] = 'M';
      else if( lines[i][column][0] == 'f' || lines[i][column][0] == 'F' )
         (*sexes)[i - first_row] = 'F';
      else
         (*sexes)[i - first_row] = 'U';
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE extractAgeColumn(
   SCIP* scip,
   int** ages,
   char*** lines,
   int num_lines,
   int first_row,
   int column
   )
{
   int i;
   SCIP_CALL( SCIPallocMemoryArray(scip, ages, num_lines - first_row) );
   for( i = first_row; i < num_lines; i++ )
      (*ages)[i - first_row] = atoi(lines[i][column]);
   return SCIP_OKAY;
}

static
SCIP_RETCODE extractGenomeColumns(
   SCIP* scip,
   char**** genomes,
   char*** lines,
   int num_lines,
   int first_row,
   int start_column,
   int end_column
   )
{
   int i;
   int j;
   int num_markers = end_column - start_column + 1;
   SCIP_CALL( SCIPallocMemoryArray(scip, genomes, num_lines - first_row) );
   for( i = first_row; i < num_lines; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*genomes)[i - first_row], num_markers) );
      for( j = start_column; j < end_column + 1; j++ )
         (*genomes)[i - first_row][j - start_column] = strdup(lines[i][j]);
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE extractAdultColumn(
   SCIP* scip,
   SCIP_Bool** adult,
   char*** lines,
   int num_lines,
   int first_row,
   int column
   )
{
   int i;
   SCIP_CALL( SCIPallocMemoryArray(scip, adult, num_lines - first_row) );
   for( i = first_row; i < num_lines; i++ )
   {
      if( lines[i][column][0] == 't' || lines[i][column][0] == 'T' )
         (*adult)[i - first_row] = TRUE;
      else
         (*adult)[i - first_row] = FALSE;
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE extractColumns(
   SCIP* scip,
   int num_lines,
   int first_row,
   char*** lines,
   int name_column,
   int sex_column,
   int age_column,
   int genome_start_column,
   int genome_end_column,
   int adult_column,
   char*** names,
   char** sexes,
   int** ages,
   char**** genomes,
   SCIP_Bool** adult
   )
{
   SCIP_CALL( extractNameColumn(scip, names, lines, num_lines, first_row, name_column) );

   SCIP_CALL( extractGenomeColumns(scip, genomes, lines, num_lines, first_row, genome_start_column, genome_end_column) );

   if( sex_column > -1 )
      SCIP_CALL( extractSexColumn(scip, sexes, lines, num_lines, first_row, sex_column) );

   if( age_column > -1 )
      SCIP_CALL( extractAgeColumn(scip, ages, lines, num_lines, first_row, age_column) );

   if( adult_column > -1 )
      SCIP_CALL( extractAdultColumn(scip, adult, lines, num_lines, first_row, adult_column) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE freeLines(
   char**** lines,
   int* num_lines,
   int** line_lengths
   )
{
   int i;
   int j;
   for( i = 0; i < (*num_lines); i++ )
   {
      for( j = 0; j < (*line_lengths)[i]; j++ )
         SCIPfreeMemoryArray(scip, &((*lines)[i][j]));
      SCIPfreeMemoryArray(scip, &((*lines)[i]));
   }
   SCIPfreeMemoryArray(scip, lines);
   SCIPfreeMemoryArray(scip, line_lengths);
   return SCIP_OKAY;
}

static
SCIP_RETCODE readDataFromFile(
   SCIP* scip,
   const char* filename,
   char* delims,
   int num_delims,
   SCIP_Bool merge_delims,
   int* first_row,
   int* num_lines,
   int* num_markers,
   char*** names,
   char** sexes,
   int** ages,
   char**** genomes,
   SCIP_Bool** adult,
   SCIP_Bool* have_ages,
   SCIP_Bool* have_sexes,
   SCIP_Bool* have_adult
   )
{
   char*** lines;
   int* line_lengths;

   SCIP_Bool has_column_headers;

   int name_column;
   int sex_column;
   int age_column;
   int genome_start_column;
   int genome_end_column;
   int adult_column;

   int min_columns_needed;

   SCIP_CALL( getColumnIndices(scip, &name_column, &sex_column, &age_column, &genome_start_column, &genome_end_column, &adult_column) );

   min_columns_needed = findNumColumnsNeeded(name_column, sex_column, age_column, genome_start_column, genome_end_column, adult_column);

   SCIPgetBoolParam(scip, "gobnilp/pedigree/inputfile/columnheaders", &has_column_headers);
   if( has_column_headers )
      (*first_row) = 1;

   (*num_markers) = genome_end_column - genome_start_column + 1;
   if( (*num_markers) / 2 * 2 != (*num_markers) )
   {
      SCIPerrorMessage("Must give an even number of columns for the genome.\n");
      return SCIP_READERROR;
   }

   SCIP_CALL( readFromFile(scip, filename, delims, num_delims, merge_delims, &lines, num_lines, &line_lengths, min_columns_needed) );

   SCIP_CALL( extractColumns(scip, (*num_lines), (*first_row), lines, name_column, sex_column, age_column, genome_start_column, genome_end_column, adult_column, names, sexes, ages, genomes, adult) );

   SCIP_CALL( freeLines(&lines, num_lines, &line_lengths) );

   (*have_ages) = (age_column != -1);
   (*have_sexes) = (sex_column != -1);
   (*have_adult) = (adult_column != -1);

   return SCIP_OKAY;
}

/** Looks up the index of an allele in a list of alleles.
 *
 *  @param allele The allele to look up.
 *  @param alleles The list of alleles to search.
 *  @param num_alleles The number of alleles in the list.
 *  @return The index of the allele in the list or -1 if it doesn't appear.
 */
static
int lookupAlleleIndex(
   char* allele,
   char** alleles,
   int num_alleles
   )
{
   int i;
   for( i = 0; i < num_alleles; i++ )
      if( strequ(allele, alleles[i]) )
         return i;
   return -1;
}
/** Finds the number of each allele of each marker in the dataset.
 *
 *  @param num_lines The number of genetic profiles in the dataset.
 *  @param num_markers The number of markers in each profile.
 *  @param genomes The genetic profiles.
 *  @param num_alleles A pointer to return the number of alleles for each marker.
 *  @param counts A pointer to return the count of each alleles for each marker.
 *  @param alleles A pointer to return the name of each alleles for each marker.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error code otherwise.
 */
static
SCIP_RETCODE countFrequencies(
   SCIP* scip,
   int num_lines,
   int num_markers,
   char*** genomes,
   char* unknown,
   int** num_alleles,
   SCIP_Real*** counts,
   char**** alleles
   )
{
   int i;
   int j;
   int* num_nonmissing_alleles;


   SCIP_CALL( SCIPallocMemoryArray(scip, &num_nonmissing_alleles, num_markers / 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, num_alleles, num_markers / 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, counts, num_markers / 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, alleles, num_markers / 2) );
   for( i = 0; i < num_markers / 2; i++ )
   {
      num_nonmissing_alleles[i] = 0;
      (*num_alleles)[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*counts)[i]), num_lines * 2) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*alleles)[i]), num_lines * 2) );
   }

   for( i = 0; i < num_lines; i++ )
   {
      for( j = 0; j < num_markers; j++ )
      {
         if (strcmp(unknown, genomes[i][j]) != 0)
         { /* i.e. skip over missing data */
            int pos = lookupAlleleIndex(genomes[i][j], (*alleles)[j / 2], (*num_alleles)[j / 2]);
            if( pos == -1 )
            {
               SCIP_CALL( SCIPallocMemoryArray(scip, &((*alleles)[j / 2][(*num_alleles)[j / 2]]), strlen(genomes[i][j]) + 1) );
               strcpy((*alleles)[j / 2][(*num_alleles)[j / 2]], genomes[i][j]);
               (*counts)[j / 2][(*num_alleles)[j / 2]] = 1;
               (*num_alleles)[j / 2] += 1;
            }
            else
            {
               (*counts)[j / 2][pos] += 1;
            }
            num_nonmissing_alleles[j / 2] += 1;
         }
      }
   }

   for( i = 0; i < num_markers / 2; i++ )
      for( j = 0; j < (*num_alleles)[i / 2]; j++ )
         (*counts)[i][j] = (*counts)[i][j] / num_nonmissing_alleles[j / 2];

   for( i = 0; i < num_markers / 2; i++ )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((*counts)[i]), (*num_alleles)[i]) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((*alleles)[i]), (*num_alleles)[i]) );
   }

   SCIPfreeMemoryArray(scip, &num_nonmissing_alleles);

   return SCIP_OKAY;
}
/** Reads in the frequencies of alleles for each marker from a file
 *
 *  @param scip The SCIP instance being used.
 *  @param frequency_file The file to read from.
 *  @param num_delims The number of field delimiters that are to be used when reading the file.
 *  @param delims The delimiters to use when reading the file.
 *  @param merge_delims Whether consecutive delimiters should be merged when reading the file.
 *  @param num_markers The number of markers in each profile.
 *  @param num_alleles A pointer to return the number of alleles for each marker.
 *  @param counts A pointer to return the count of each alleles for each marker.
 *  @param alleles A pointer to return the name of each alleles for each marker.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error code otherwise.
 *
*/
static
SCIP_RETCODE readFrequencyFile(
   SCIP* scip,
   const char* frequency_file,
   int num_delims,
   char* delims,
   SCIP_Bool merge_delims,
   int num_markers,
   int** num_alleles,
   SCIP_Real*** counts,
   char**** alleles
   )
{
   int i;
   int j;
   FILE* file;

   char*** lines;
   int num_lines;
   int* line_lengths;

   int marker_entry = 0;
   SCIP_Bool names_next = TRUE;

   SCIP_CALL( SCIPallocMemoryArray(scip, num_alleles, num_markers / 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, counts, num_markers / 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, alleles, num_markers / 2) );

   file = fopen(frequency_file, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", frequency_file);
      return SCIP_NOFILE;
   }
   SCIP_CALL( UT_readFileAndSplit(scip, file, delims, num_delims, merge_delims, &lines, &num_lines, &line_lengths) );
   fclose(file);

   for( i = 0; i < num_lines; i++ )
   {
      if( line_lengths[i] == 1 && strequ(lines[i][0], "") )
      {
         /* Blank line */
      }
      else if( names_next )
      {
         if( marker_entry == num_markers / 2 )
         {
            SCIPerrorMessage("Data for more markers in file %s than expected.\n", frequency_file);
            return SCIP_READERROR;
         }
         (*num_alleles)[marker_entry] = line_lengths[i];
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*counts)[marker_entry]), line_lengths[i]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*alleles)[marker_entry]), line_lengths[i]) );
         for( j = 0; j < line_lengths[i]; j++ )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*alleles)[marker_entry][j]), strlen(lines[i][j]) + 1) );
            strcpy((*alleles)[marker_entry][j], lines[i][j]);
         }
         names_next = FALSE;
      }
      else
      {
         if( line_lengths[i] != line_lengths[i - 1] )
         {
            SCIPerrorMessage("Error on line %d of %s.  Expected %d items but found %d.\n", i + 1, frequency_file, line_lengths[i - 1], line_lengths[i]);
            return SCIP_READERROR;
         }
         for( j = 0; j < line_lengths[i]; j++ )
            (*counts)[marker_entry][j] = atof(lines[i][j]);
         marker_entry += 1;
         names_next = TRUE;
      }
   }

   if( !names_next )
   {
      SCIPerrorMessage("Different number of names and frequency lines in %s.\n", frequency_file);
      return SCIP_READERROR;
   }

   for( i = 0; i < num_lines; i++ )
   {
      for( j = 0; j < line_lengths[i]; j++ )
         SCIPfreeMemoryArray(scip, &(lines[i][j]));
      SCIPfreeMemoryArray(scip, &(lines[i]));
   }
   SCIPfreeMemoryArray(scip, &lines);
   SCIPfreeMemoryArray(scip, &line_lengths);

   return SCIP_OKAY;
}

/** Determines whether an individual could be the parent of another individual.
 *
 *  @param child The index of the potential child.
 *  @param parent The index of the potential parent.
 *  @param num_markers The number of markers in the genetic profile.
 *  @param genomes The genetic profiles.
 *  @param genomes The ages of each of the individuals.
 *  @return True if @c parent is genetically capable of being the parent of @c child, or False otherwise.
 */
static
SCIP_Bool isPossibleParentOf(
   int child,
   int parent,
   int num_markers,
   char*** genomes,
   int* ages,
   SCIP_Bool* adult,
   int max_mismatches,
   char* unknown,
   int maxparentchildagegap
   )
{
   int i;
   int mismatches_sofar = 0;
   
   if( adult != NULL && adult[parent] == FALSE )
      /* The parent is not an adult */
      return FALSE;
   if( ages != NULL )
   {
      int agegap;
      agegap = ages[parent] - ages[child];
      /* The parent is not older than the child */
      /* The parent is not more than maxparentchildagegarp older than the child */
      if( agegap <= 0 || (maxparentchildagegap != -1 && agegap > maxparentchildagegap) )
         return FALSE;
   }
   for( i = 0; i < num_markers; i += 2 )
   {
      char* child1 = genomes[child][i];
      char* child2 = genomes[child][i + 1];
      char* parent1 = genomes[parent][i];
      char* parent2 = genomes[parent][i + 1];
      if ( strequ(child1, unknown) || strequ(child2, unknown) )
      {
         /* Child has an unknown allele - could be parent */
      }
      else if( strequ(child1, parent1) || strequ(child1, parent2) || strequ(child2, parent1) || strequ(child2, parent2))
      {
         /* An allele matched - could be parent */
      }
      else
      {
         /* Neither allele matched - have to consider mutation or genotyping error */
         mismatches_sofar += 1;
         if (mismatches_sofar > max_mismatches)
            return FALSE;
      }
   }
   return TRUE;
}

static
int get_idx(
   char** names,
   int num_lines,
   char* child
   )
{
   int i;

   for( i = 0; i < num_lines; i++ )
   {
      if( strcmp(names[i], child) == 0 )
         return i;
   }
   return -1;
}

static 
SCIP_RETCODE findAllowedParents(
   SCIP* scip,
   int num_lines,
   char* pedigreeconstraintsfile,
   char** names,
   SCIP_Bool*** allowed_parent_ptr,
   int** num_allowed_parent_ptr
   )
{
   
   int status;
   char line[SCIP_MAXSTRLEN];
   FILE* pedigreeconstraints;

   int i;

   char child[SCIP_MAXSTRLEN];
   char parent1[SCIP_MAXSTRLEN];
   char parent2[SCIP_MAXSTRLEN];
   
   int child_index;
   int parent1_index;
   int parent2_index;


   pedigreeconstraints = fopen(pedigreeconstraintsfile, "r");
   if( pedigreeconstraints == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", pedigreeconstraintsfile);
      return SCIP_NOFILE;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, num_allowed_parent_ptr, num_lines) );
   for( i = 0; i < num_lines; ++i )
      (*num_allowed_parent_ptr)[i] = -1;      


   status = fscanf(pedigreeconstraints, "%[^\n]%*c", line);
   while( status == 1 )
   {
      if( sscanf(line, "%[^<]<=", child) == 1 )
      {
         child_index = get_idx(names, num_lines, child);
         
         assert(child_index > -1);

         if( (*allowed_parent_ptr) == NULL )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, allowed_parent_ptr, num_lines) );
            for( i = 0; i < num_lines; ++i )
               (*allowed_parent_ptr)[i] = NULL;
         }

         if( (*allowed_parent_ptr)[child_index] == NULL )
         {
            (*num_allowed_parent_ptr)[child_index] = 0;
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*allowed_parent_ptr)[child_index]), num_lines) );
            for( i = 0; i < num_lines; ++i )
               (*allowed_parent_ptr)[child_index][i] = FALSE;
         }
         if( sscanf(line, "%[^<]<=%[^,],%s", child, parent1, parent2) == 3 )
         {

            (*num_allowed_parent_ptr)[child_index] = 2;
            parent1_index = get_idx(names, num_lines, parent1);
            parent2_index = get_idx(names, num_lines, parent2);

            assert(parent1_index > -1);
            assert(parent2_index > -1);

            (*allowed_parent_ptr)[child_index][parent1_index] = TRUE;
            (*allowed_parent_ptr)[child_index][parent2_index] = TRUE;
         }
         else if( sscanf(line, "%[^<]<=%s", child, parent1) == 2 )
         {

            (*num_allowed_parent_ptr)[child_index] = 1;
            parent1_index = get_idx(names, num_lines, parent1);
            assert(parent1_index > -1);
            (*allowed_parent_ptr)[child_index][parent1_index] = TRUE;
         }
      }
      status = fscanf(pedigreeconstraints, "%[^\n]%*c", line);
   }

   fclose(pedigreeconstraints);
   return SCIP_OKAY;
}
    

/** Finds all the individuals capable of being the parent of each individual in the dataset.
 *
 *  @param num_lines The number of individuals in the dataset.
 *  @param num_markers The number of markers in each genetic profile.
 *  @param genomes The genetic profiles.
 *  @param genomes The ages of each of the individuals.
 *  @param num_possible_parents A pointer to return the number of possible parents for each individual in the dataset.
 *  @param possible_parents A pointer to return the indices of each potential parent of each individual.
 *
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error code otherwise.
 */
static SCIP_RETCODE findPossibleParents(
   SCIP* scip,
   int num_lines, 
   int num_markers, 
   char*** genomes, 
   int* ages, 
   SCIP_Bool* adult, 
   int max_mismatches, 
   char* unknown, 
   int** num_possible_parents, 
   int*** possible_parents,
   SCIP_Bool** allowed_parent
   )
{
   int i;
   int j;

   int maxparentchildagegap;

   SCIPgetIntParam(scip, "gobnilp/pedigree/maxparentchildagegap", &maxparentchildagegap);
   
   SCIP_CALL( SCIPallocMemoryArray(scip, num_possible_parents, num_lines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, possible_parents, num_lines) );
   for( i = 0; i < num_lines; i++ )
   {
      (*num_possible_parents)[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*possible_parents)[i]), num_lines) );
   }

   for( i = 0; i < num_lines; i++ )
      for( j = 0; j < num_lines; j++ )
         if( (i != j) && (allowed_parent == NULL || allowed_parent[i] == NULL || allowed_parent[i][j]) && isPossibleParentOf(i, j, num_markers, genomes, ages, adult, max_mismatches, unknown, maxparentchildagegap) )
         {
            (*possible_parents)[i][(*num_possible_parents)[i]] = j;
            (*num_possible_parents)[i] += 1;
         }

   for( i = 0; i < num_lines; i++ )
   {
      /* printf("%d %d\n", i, (*num_possible_parents)[i]); */
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((*possible_parents)[i]), (*num_possible_parents)[i]) );
   }

   return SCIP_OKAY;
}
/** Finds all possible parent sets for each individual.
 *
 *  @param num_lines The number of individuals in the dataset.
 *  @param num_possible_parents The number of possible parents for each individual in the dataset.
 *  @param possible_parents The indices of each potential parent of each individual.
 *  @param allow_single_parents Whether parent sets of size 1 are permitted.
 *  @param has_sex_information Whether there is sex information to use.
 *  @param sexes The sex information of each individual, if known.
 *  @param num_parent_sets A pointer to return the number of parent sets for each individual
 *  @param num_parents_in_set A pointer to return the number of parents in each of an individual's parent set
 *  @param parent_sets The parents sets for each individual.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error code otherwise.
 */
static
SCIP_RETCODE findPossibleParentSets(
   SCIP* scip,
   int num_lines,
   int* num_possible_parents,
   int** possible_parents,
   SCIP_Bool allow_single_parents,
   SCIP_Bool has_sex_information,
   char* sexes,
   int** num_parent_sets,
   int*** num_parents_in_set,
   int**** parent_sets,
   int* num_allowed_parent
   )
{
   int i;
   int j;
   int k;

   SCIP_CALL( SCIPallocMemoryArray(scip, num_parent_sets, num_lines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, num_parents_in_set, num_lines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, parent_sets, num_lines) );
   for( i = 0; i < num_lines; i++ )
   {
      (*num_parent_sets)[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*num_parents_in_set)[i]), 1 + num_possible_parents[i] + num_possible_parents[i]*num_possible_parents[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*parent_sets)[i]), 1 + num_possible_parents[i] + num_possible_parents[i]*num_possible_parents[i]) );
   }

   for( i = 0; i < num_lines; i++ )
   {
      if( num_allowed_parent == NULL || num_allowed_parent[i] < 1 )
      {
         /* Could be a founder */
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*parent_sets)[i][(*num_parent_sets)[i]]), 0) );
         (*num_parents_in_set)[i][(*num_parent_sets)[i]] = 0;
         (*num_parent_sets)[i] += 1;
      }
      for( j = 0; j < num_possible_parents[i]; j++ )
      {
         if( allow_single_parents && (num_allowed_parent == NULL || num_allowed_parent[i] == -1 || num_allowed_parent[i] == 1))
         {
            /* Could have just j as a parent */
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*parent_sets)[i][(*num_parent_sets)[i]]), 1) );
            (*parent_sets)[i][(*num_parent_sets)[i]][0] = possible_parents[i][j];
            (*num_parents_in_set)[i][(*num_parent_sets)[i]] = 1;
            (*num_parent_sets)[i] += 1;
         }
         for( k = 0; k < num_possible_parents[i]; k++ )
         {
            if( (j < k) && (num_allowed_parent == NULL || num_allowed_parent[i] == -1 || num_allowed_parent[i] == 2) && (!has_sex_information || sexes[possible_parents[i][j]] == 'U' || sexes[possible_parents[i][j]] != sexes[possible_parents[i][k]]) )
            {
               /* Could have just j and k as parents */
               SCIP_CALL( SCIPallocMemoryArray(scip, &((*parent_sets)[i][(*num_parent_sets)[i]]), 2) );
               (*parent_sets)[i][(*num_parent_sets)[i]][0] = possible_parents[i][j];
               (*parent_sets)[i][(*num_parent_sets)[i]][1] = possible_parents[i][k];
               (*num_parents_in_set)[i][(*num_parent_sets)[i]] = 2;
               (*num_parent_sets)[i] += 1;
            }
         }
      }
   }

   for( i = 0; i < num_lines; i++ )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((*num_parents_in_set)[i]), (*num_parent_sets)[i]) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((*parent_sets)[i]), (*num_parent_sets)[i]) );
   }

   return SCIP_OKAY;
}

static
SCIP_Real getCountWithErrorMessage(
   int line,
   char* allele,
   int* num_alleles,
   SCIP_Real** allele_counts,
   char*** allele_names,
   int index,
   int position,
   char* unknown
   )
{
   if (strequ(allele, unknown))
   {
      return 1;
   }
   else
   {
      int lookup = lookupAlleleIndex(allele, allele_names[index], num_alleles[index]);
      if (lookup == -1)
      {
         SCIPerrorMessage("Error: Unknown allele frequency at line %d for %s at profile position %d\n", line, allele, position);
         exit(-1);
      }
      return allele_counts[index][lookup];
   }
}

static
SCIP_Bool isHomozygotic(
   char* allele1,
   char* allele2
   )
{
   return strequ(allele1, allele2);
}

static
void swapString(
   char** str1,
   char** str2
   )
{
   char* tmp = *str1;
   *str1 = *str2;
   *str2 = tmp;
}

static
SCIP_Bool sortAlleles(
   char** allele1,
   char** allele2,
   char* unknown
   )
{
   if ( strequ(*allele1, unknown) )
   {
      if ( strequ(*allele2, unknown) )
      {
         /* Both unknown - no need to do anything */
         return FALSE;
      }
      else
      {
         /* allele1 unknown, allele2 known - need to swap */
         swapString(allele1, allele2);
         return TRUE;
      }
   }
   else
   {
      if ( strequ(*allele2, unknown) )
      {
         /* allele1 known, allele2 unknown - no need to do anything */
         return FALSE;
      }
      else
      {
         /* both known - need to swap if out of order*/
         if ( strcmp(*allele1, *allele2) > 0 )
         {
            swapString(allele1, allele2);
            return TRUE;
         }
         else
         {
            return FALSE;
         }
      }
   }
}

static
void sortAllelesAndCounts(
   char** allele1,
   char** allele2,
   SCIP_Real* count1,
   SCIP_Real* count2,
   char* unknown
   )
{
   SCIP_Bool swapped = sortAlleles(allele1, allele2, unknown);
   if ( swapped )
   {
      SCIP_Real tmp = *count1;
      *count1 = *count2;
      *count2 = tmp;
   }
}

static
SCIP_Real NONEtoAA(
   SCIP_Real probA
   )
{
   /* = p(pick A then pick A) */
   return probA*probA;
}

static
SCIP_Real NONEtoAB(
   SCIP_Real probA,
   SCIP_Real probB
   )
{
   /* = p(pick A then pick B, or pick B then pick A) */
   return 2*probA*probB;
}

static
SCIP_Real AAtoAA(
   SCIP_Real probA,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A without mutation and pick A) */

   return no_mutation * probA ;
}

static
SCIP_Real AAtoAB(
   SCIP_Real probA,
   SCIP_Real probB,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A without mutation and pick B, or
          inherit A with mutation to B and pick A) */

   return no_mutation * probB +
          mutation * probA ;
}

static
SCIP_Real AAtoBB(
   SCIP_Real probA,
   SCIP_Real probB,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A with mutation to B and pick B) */

   return mutation * probB ;
}

static
SCIP_Real AAtoBC(
   SCIP_Real probA,
   SCIP_Real probB,
   SCIP_Real probC,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A with mutation to B and pick C, or
          inherit A with mutation to C and pick B) */

   /*return mutation * probC +
          mutation * probB ;*/
   return mutation * (probC + probB) ;
}

static
SCIP_Real ABtoAA(
   SCIP_Real probA,
   SCIP_Real probB,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A with no mutation and pick A, or
          inherit B with mutation to A and pick A) */

   /*return 0.5 * no_mutation * probA +
          0.5 * mutation * probA ; */
   return 0.5 * (no_mutation + mutation) * probA ;
}

static
SCIP_Real ABtoCC(
   SCIP_Real probA,
   SCIP_Real probB,
   SCIP_Real probC,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A with mutation to C and pick C, or
          inherit B with mutation to C and pick C) */

   /*return 0.5 * mutation * probC +
          0.5 * mutation * probC ;*/
   return mutation * probC ;
}

static
SCIP_Real ABtoAB(
   SCIP_Real probA,
   SCIP_Real probB,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A with no mutation and pick B, or
          inherit A with mutation to B and pick A, or
          inherit B with no mutation and pick A, or
          inherit B with mutation to A and pick B) */

   /*return 0.5 * no_mutation * probB +
          0.5 * mutation * probA +
          0.5 * no_mutation * probA +
          0.5 * mutation * probB ;*/
   return 0.5 * (mutation + no_mutation) * (probA + probB) ;
}

static
SCIP_Real ABtoAC(
   SCIP_Real probA,
   SCIP_Real probB,
   SCIP_Real probC,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A with no mutation and pick C, or
          inherit A with mutation to C and pick A, or
          inherit B with mutation to A and pick C, or
          inherit B with mutation to C and pick A) */

   /*return 0.5 * no_mutation * probC +
          0.5 * mutation * probA +
          0.5 * mutation * probC +
          0.5 * mutation * probA ;*/
   return 0.5 * (no_mutation + mutation) * probC + mutation * probA ;
}

static
SCIP_Real ABtoCD(
   SCIP_Real probA,
   SCIP_Real probB,
   SCIP_Real probC,
   SCIP_Real probD,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* = p(inherit A with mutation to C and pick D, or
          inherit A with mutation to D and pick C, or
          inherit B with mutation to C and pick D, or
          inherit B with mutation to D and pick C) */

   /*return 0.5 * mutation * probD +
          0.5 * mutation * probC +
          0.5 * mutation * probD +
          0.5 * mutation * probC ;*/
   return mutation * (probC + probD) ;
}


static
SCIP_Real AAbecomesAA(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 stays as A and 2 stays as A */
   return no_mutation * no_mutation;
}

static
SCIP_Real AAbecomesBB(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 becomes B and 2 becomes B */
   return mutation * mutation;
}

static
SCIP_Real AAbecomesAB(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 stays as A and 2 becomes B, or
      1 becomes B and 2 stays as A */
   return 2 * no_mutation * mutation;
}

static
SCIP_Real AAbecomesBC(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 becomes B and 2 becomes C, or
      1 becomes C and 2 becomes B */
   return 2 * mutation * mutation ;
}

static
SCIP_Real ABbecomesAA(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 stays as A and 2 becomes A */
   return no_mutation * mutation ;
}

static
SCIP_Real ABbecomesCC(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 becomes C and 2 becomes C */
   return mutation * mutation ;
}

static
SCIP_Real ABbecomesAB(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 stays as A and 2 stays as B, or
      1 becomes B and 2 becomes A */
   return (no_mutation * no_mutation) + (mutation * mutation) ;
}

static
SCIP_Real ABbecomesAC(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 stays as A and 2 becomes C, or
      1 becomes C and 2 becomes A */
   return (no_mutation * mutation) + (mutation * mutation) ;
}

static
SCIP_Real ABbecomesCD(
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   /* 1 becomes C and 2 becomes D, or
      1 becomes D and 2 becomes C */
   return 2 * mutation * mutation ;
}


static
SCIP_Real mutationProbability(
   char* inherited1,
   char* inherited2,
   char* observed1,
   char* observed2,
   SCIP_Real mutation,
   SCIP_Real no_mutation
   )
{
   if ( isHomozygotic(inherited1, inherited2) )
      if ( isHomozygotic(observed1, observed2) )
         /* homo -> homo */
         if ( strequ(inherited1, observed1) )
            return AAbecomesAA(mutation, no_mutation);
         else
            return AAbecomesBB(mutation, no_mutation);
      else
         /* homo -> hetro */
         if ( strequ(inherited1, observed1) || strequ(inherited1, observed2) )
            return AAbecomesAB(mutation, no_mutation);
         else
            return AAbecomesBC(mutation, no_mutation);
   else
      if ( isHomozygotic(observed1, observed2) )
         /* hetro -> homo */
         if ( strequ(inherited1, observed1) || strequ(inherited2, observed1) )
            return ABbecomesAA(mutation, no_mutation);
         else
            return ABbecomesCC(mutation, no_mutation);
      else
         /* hetro -> hetro */
         if ( strequ(inherited1, observed1) || strequ(inherited1, observed2) )
            if ( strequ(inherited2, observed1) || strequ(inherited2, observed2) )
               return ABbecomesAB(mutation, no_mutation);
            else
               return ABbecomesAC(mutation, no_mutation);
         else
            if ( strequ(inherited2, observed1) || strequ(inherited2, observed2) )
               return ABbecomesAC(mutation, no_mutation);
            else
               return ABbecomesCD(mutation, no_mutation);
}


static
SCIP_Real oneMarkerLikelihoodNoParents(
   char* allele1,
   char* allele2,
   SCIP_Real count1,
   SCIP_Real count2,
   char* unknown
   )
{
   if( isHomozygotic(allele1, allele2) )
      return log( NONEtoAA(count1) );
   else
      return log( NONEtoAB(count1,count2) );
}

static
SCIP_Real oneMarkerLikelihoodOneParent(
   char* child1,
   char* child2,
   SCIP_Real count1,
   SCIP_Real count2,
   char* parent1,
   char* parent2,
   SCIP_Real pcount1,
   SCIP_Real pcount2,
   char* unknown,
   SCIP_Real mutation,
   int num_alleles
   )
{
   SCIP_Real total;
   SCIP_Real no_mutation = 1 - mutation ;
   SCIP_Real specific_mutation = mutation / (num_alleles-1) ;

   if ( strequ(parent1, unknown) )
   {
      /* Parent has unknown genotype, so count as if child was founder */
      total = oneMarkerLikelihoodNoParents(child1, child2, count1, count2, unknown);
   }
   else if (strequ(child1, unknown) )
   {
      /* Child has unknown genotype */
      total = log(1);
   }
   else if ( isHomozygotic(parent1, parent2) )
   {
      if ( isHomozygotic(child1, child2) )
      {
         if ( strequ(parent1, child1) )
         {
            /* AA -> AA */
            total = log( AAtoAA(count1,specific_mutation,no_mutation) );
         }
         else
         {
            /* AA -> BB */
            total = log( AAtoBB(pcount1,count1,specific_mutation,no_mutation) );
         }
      }
      else
      {
         if ( strequ(parent1, child1) )
         {
            /* AA -> AB */
            total = log( AAtoAB(count1,count2,specific_mutation,no_mutation) );
         }
         else if ( strequ(parent1, child2) )
         {
            /* AA -> BA */
            total = log( AAtoAB(count2,count1,specific_mutation,no_mutation) );
         }
         else
         {
            /* AA -> BC */
            total = log( AAtoBC(pcount1,count1,count2,specific_mutation,no_mutation) );
         }
      }
   }
   else
   {
      if ( isHomozygotic(child1, child2) )
      {
         if ( strequ(parent1, child1) )
         {
            /* AB -> AA */
            total = log( ABtoAA(count1,pcount2,specific_mutation,no_mutation) );
         }
         else if ( strequ(parent2, child1) )
         {
            /* AB -> BB */
            total = log( ABtoAA(count1,pcount1,specific_mutation,no_mutation) );
         }
         else
         {
            /* AB -> CC */
            total = log( ABtoCC(pcount1,pcount2,count1,specific_mutation,no_mutation) );
         }
      }
      else
      {
         if ( strequ(parent1, child1) )
         {
            if ( strequ(parent2, child2) )
            {
               /* AB -> AB */
               total = log( ABtoAB(count1,count2,specific_mutation,no_mutation) );
            }
            else
            {
               /* AB -> AC */
               total = log( ABtoAC(count1,pcount2,count2,specific_mutation,no_mutation) );
            }
         }
         else if ( strequ(parent1, child2) )
         {
            /* AB -> CA */
            total = log( ABtoAC(count2,pcount2,count1,specific_mutation,no_mutation) );
         }
         else if ( strequ(parent2, child1) )
         {
            /* AB -> BC */
            total = log( ABtoAC(count1,pcount1,count2,specific_mutation,no_mutation) );
         }
         else if ( strequ(parent2, child2) )
         {
            /* AB -> CB */
            total = log( ABtoAC(count2,pcount1,count1,specific_mutation,no_mutation) );
         }
         else
         {
            /* AB -> CD */
            total = log( ABtoCD(pcount1,pcount2,count1,count2,specific_mutation,no_mutation) );
         }
      }
   }
   return total;
}

static
SCIP_Real oneMarkerLikelihoodTwoParents(
   char* childa,
   char* childb,
   SCIP_Real count1,
   SCIP_Real count2,
   char* parent1a,
   char* parent1b,
   SCIP_Real pcount1a,
   SCIP_Real pcount1b,
   char* parent2a,
   char* parent2b,
   SCIP_Real pcount2a,
   SCIP_Real pcount2b,
   char* unknown,
   SCIP_Real mutation,
   int num_alleles
   )
{
   SCIP_Real total;

   if ( strequ(parent1a, unknown) )
   {
      if ( strequ(parent2a, unknown) )
      {
         /* Both parents have unknown genotype, so count as if child was founder */
         total = oneMarkerLikelihoodNoParents(childa, childb, count1, count2, unknown);
      }
      else
      {
         /* Parent 1 has unknown genotype, so count as if child only had parent 2 */
         total = oneMarkerLikelihoodOneParent(childa, childb, count1, count2, parent2a, parent2b, pcount2a, pcount2b, unknown, mutation, num_alleles);
      }
   }
   else if ( strequ(parent2a, unknown) )
   {
      /* Parent 2 has unknown genotype, so count as if child only had parent 1 */
      total = oneMarkerLikelihoodOneParent(childa, childb, count1, count2, parent1a, parent1b, pcount1a, pcount1b, unknown, mutation, num_alleles);
   }
   else if (strequ(childa, unknown) )
   {
      /* Child has unknown genotype */
      total = log(1);
   }
   else
   {
      /* Work out the Punnett square */
      SCIP_Real no_mutation = 1 - mutation ;
      SCIP_Real specific_mutation = mutation / (num_alleles-1) ;
      SCIP_Real p_child1 = 0.25 * mutationProbability(parent1a, parent2a, childa, childb, specific_mutation, no_mutation);
      SCIP_Real p_child2 = 0.25 * mutationProbability(parent1a, parent2b, childa, childb, specific_mutation, no_mutation);
      SCIP_Real p_child3 = 0.25 * mutationProbability(parent1b, parent2a, childa, childb, specific_mutation, no_mutation);
      SCIP_Real p_child4 = 0.25 * mutationProbability(parent1b, parent2b, childa, childb, specific_mutation, no_mutation);
      total = log(p_child1 + p_child2 + p_child3 + p_child4);
   }

   return total;
}

/** Computes the log likelihood of an individual having no parents.
 *
 *  @param num_markers The number of markers in the genetic profile.
 *  @param genome The genetic profile.
 *  @param num_alleles The number of alleles for each marker.
 *  @param allele_counts The count of each allele of each marker.
 *  @param allele_names The name of each allele of each marker.
 *  @return The log likelihood of the individual having no parents.
 */
static
SCIP_Real computeLikelihoodNoParents(
   int line,
   int num_markers,
   char** genome,
   int* num_alleles,
   SCIP_Real** allele_counts,
   char*** allele_names,
   char* unknown
   )
{
   int i;
   SCIP_Real total = 0;
   for( i = 0; i < num_markers / 2; i++ )
   {
      char* allele1 = genome[2 * i];
      char* allele2 = genome[2 * i + 1];
      SCIP_Real count1 = getCountWithErrorMessage(line, allele1, num_alleles, allele_counts, allele_names, i, 2*i, unknown);
      SCIP_Real count2 = getCountWithErrorMessage(line, allele2, num_alleles, allele_counts, allele_names, i, 2*i+1, unknown);
      sortAllelesAndCounts(&allele1, &allele2, &count1, &count2, unknown);

      if ( !strequ(allele1, unknown) && strequ(allele2, unknown) )
      {
         SCIPerrorMessage("Error: Cannot currently handle genomes where only one of the allele pair is missing.\n");
         exit(-1);
      }

      total += oneMarkerLikelihoodNoParents(allele1, allele2, count1, count2, unknown);
      if (!isfinite(total))
         return total;
   }
   return total;
}

/** Computes the log likelihood of an individual having one parent.
 *
 *  @param num_markers The number of markers in the genetic profile.
 *  @param child_genome The genetic profile of the child.
 *  @param parent_genome The genetic profile of the parent.
 *  @param num_alleles The number of alleles for each marker.
 *  @param allele_counts The count of each allele of each marker.
 *  @param allele_names The name of each allele of each marker.
 *  @return The log likelihood of the individual having the given parent.
 */
static
SCIP_Real computeLikelihoodOneParent(
   int line,
   int num_markers,
   char** child_genome,
   char** parent_genome,
   int* num_alleles,
   SCIP_Real** allele_counts,
   char*** allele_names,
   int max_mismatches,
   char* unknown,
   int mutation
   )
{
   int i;
   SCIP_Real total = 0;
   for( i = 0; i < num_markers / 2; i++ )
   {
      char* child1 = child_genome[2 * i];
      char* child2 = child_genome[2 * i + 1];
      char* parent1 = parent_genome[2 * i];
      char* parent2 = parent_genome[2 * i + 1];
      SCIP_Real count1 = getCountWithErrorMessage(line, child1, num_alleles, allele_counts, allele_names, i, 2*i, unknown);
      SCIP_Real count2 = getCountWithErrorMessage(line, child2, num_alleles, allele_counts, allele_names, i, 2*i+1, unknown);
      SCIP_Real pcount1 = getCountWithErrorMessage(line, parent1, num_alleles, allele_counts, allele_names, i, 2*i, unknown);
      SCIP_Real pcount2 = getCountWithErrorMessage(line, parent2, num_alleles, allele_counts, allele_names, i, 2*i+1, unknown);
      sortAllelesAndCounts(&child1, &child2, &count1, &count2, unknown);
      sortAllelesAndCounts(&parent1, &parent2, &pcount1, &pcount2, unknown);

      if ( (!strequ(parent1, unknown) && strequ(parent2, unknown)) || (!strequ(child1, unknown) && strequ(child2, unknown)) )
      {
         SCIPerrorMessage("Error: Cannot currently handle genomes where only one of the allele pair is missing.\n");
         exit(-1);
      }

      total += oneMarkerLikelihoodOneParent(child1, child2, count1, count2, parent1, parent2, pcount1, pcount2, unknown, mutation, num_alleles[i]);
      if (!isfinite(total))
         return total;

   }
   return total;
}

/** Computes the log likelihood of an individual having two parents.
 *
 *  @param num_markers The number of markers in the genetic profile.
 *  @param child_genome The genetic profile of the child.
 *  @param parent1_genome The genetic profile of the first parent.
 *  @param parent2_genome The genetic profile of the second parent.
 *  @return The log likelihood of the individual having the given parents.
 */
static
SCIP_Real computeLikelihoodTwoParents(
   int line,
   int num_markers,
   char** child_genome,
   char** parent1_genome,
   char** parent2_genome,
   int* num_alleles,
   SCIP_Real** allele_counts,
   char*** allele_names,
   int max_mismatches,
   char* unknown,
   int mutation
   )
{
   int i;
   SCIP_Real total = 0;
   for( i = 0; i < num_markers / 2; i++ )
   {
      char* childa = child_genome[2 * i];
      char* childb = child_genome[2 * i + 1];
      char* parent1a = parent1_genome[2 * i];
      char* parent1b = parent1_genome[2 * i + 1];
      char* parent2a = parent2_genome[2 * i];
      char* parent2b = parent2_genome[2 * i + 1];
      SCIP_Real counta = getCountWithErrorMessage(line, childa, num_alleles, allele_counts, allele_names, i, 2*i, unknown);
      SCIP_Real countb = getCountWithErrorMessage(line, childa, num_alleles, allele_counts, allele_names, i, 2*i+1, unknown);
      SCIP_Real pcount1a = getCountWithErrorMessage(line, parent1a, num_alleles, allele_counts, allele_names, i, 2*i, unknown);
      SCIP_Real pcount1b = getCountWithErrorMessage(line, parent1b, num_alleles, allele_counts, allele_names, i, 2*i+1, unknown);
      SCIP_Real pcount2a = getCountWithErrorMessage(line, parent2a, num_alleles, allele_counts, allele_names, i, 2*i, unknown);
      SCIP_Real pcount2b = getCountWithErrorMessage(line, parent2b, num_alleles, allele_counts, allele_names, i, 2*i+1, unknown);
      sortAllelesAndCounts(&childa, &childa, &counta, &countb, unknown);
      sortAllelesAndCounts(&parent1a, &parent1b, &pcount1a, &pcount1b, unknown);
      sortAllelesAndCounts(&parent2a, &parent2b, &pcount2a, &pcount2b, unknown);

      sortAlleles(&childa, &childb, unknown);
      sortAlleles(&parent1a, &parent1b, unknown);
      sortAlleles(&parent2a, &parent2b, unknown);

      if ( (!strequ(parent1a, unknown) && strequ(parent1b, unknown)) || (!strequ(parent2a, unknown) && strequ(parent2b, unknown)) || (!strequ(childa, unknown) && strequ(childb, unknown)) )
      {
         SCIPerrorMessage("Error: Cannot currently handle genomes where only one of the allele pair is missing.\n");
         exit(-1);
      }

      total += oneMarkerLikelihoodTwoParents(childa, childb, counta, countb, parent1a, parent1b, pcount1a, pcount1b, parent2a, parent2b, pcount2a, pcount2b, unknown, mutation, num_alleles[i]);
      if (!isfinite(total))
         return total;
   }
   return total;
}

/** Computes the log likelihoods of each possible parent set for each individual.
 *
 *  @param num_lines The number of individuals in the dataset.
 *  @param num_markers The number of markers in the genetic profile.
 *  @param genomes The genetic profiles.
 *  @param num_parent_sets The number of parent sets for each individual.
 *  @param num_parents_in_set The number of parents in each of an individual's parent set.
 *  @param parent_sets The parents sets for each individual.
 *  @param num_alleles The number of alleles for each marker.
 *  @param allele_counts The count of each allele of each marker.
 *  @param allele_names The name of each allele of each marker.
 *  @param likelihoods A pointer to return the likelihood of each possible parent set.
 *  @param new_num_parent_sets A pointer to return the updated number of parent sets for each individual after those which are not possible have been removed.
 *  @param new_num_parents_in_set A pointer to return the updated number of parents in each of an individual's parent set after those which are not possible have been removed.
 *  @param new_parent_sets A pointer to return the updated parents sets for each individual after those which are not possible have been removed.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error code otherwise.
 */
static
SCIP_RETCODE computeLikelihoods(
   int num_lines,
   int num_markers,
   char*** genomes,
   int max_mismatches,
   char* unknown,
   int* num_parent_sets,
   int** num_parents_in_set,
   int*** parent_sets,
   int* num_alleles,
   SCIP_Real** allele_counts,
   char*** allele_names,
   SCIP_Real*** likelihoods,
   int** new_num_parent_sets,
   int*** new_num_parents_in_set,
   int**** new_parent_sets,
   int mutation
   )
{
   int i;
   int j;
   int k;

   SCIP_CALL( SCIPallocMemoryArray(scip, new_num_parent_sets, num_lines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, new_num_parents_in_set, num_lines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, new_parent_sets, num_lines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, likelihoods, num_lines) );
   for( i = 0; i < num_lines; i++ )
   {
      (*new_num_parent_sets)[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*likelihoods)[i]), num_parent_sets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*new_num_parents_in_set)[i]), num_parent_sets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*new_parent_sets)[i]), num_parent_sets[i]) );
   }

   for( i = 0; i < num_lines; i++ )
   {
      for( j = 0; j < num_parent_sets[i]; j++ )
      {
         SCIP_Real likelihood = 0;
         if( num_parents_in_set[i][j] == 0 )
            likelihood = computeLikelihoodNoParents(i,num_markers, genomes[i], num_alleles, allele_counts, allele_names, unknown);
         else if( num_parents_in_set[i][j] == 1 )
            likelihood = computeLikelihoodOneParent(i,num_markers, genomes[i], genomes[parent_sets[i][j][0]], num_alleles, allele_counts, allele_names, max_mismatches, unknown, mutation);
         else
            likelihood = computeLikelihoodTwoParents(i,num_markers, genomes[i], genomes[parent_sets[i][j][0]], genomes[parent_sets[i][j][1]], num_alleles, allele_counts, allele_names, max_mismatches, unknown, mutation);
         if( isfinite(likelihood) )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*new_parent_sets)[i][(*new_num_parent_sets)[i]]), num_parents_in_set[i][j]) );
            if( num_parents_in_set[i][j] > 0 )
               for( k = 0; k < num_parents_in_set[i][j]; k++ )
                  (*new_parent_sets)[i][(*new_num_parent_sets)[i]][k] = parent_sets[i][j][k];
            (*new_num_parents_in_set)[i][(*new_num_parent_sets)[i]] = num_parents_in_set[i][j];
            (*likelihoods)[i][(*new_num_parent_sets)[i]] = likelihood;
            (*new_num_parent_sets)[i] += 1;
         }
      }
   }

   for( i = 0; i < num_lines; i++ )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((*likelihoods)[i]), (*new_num_parent_sets)[i]) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((*new_num_parents_in_set)[i]), (*new_num_parent_sets)[i]) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((*new_parent_sets)[i]), (*new_num_parent_sets)[i]) );
   }

   return SCIP_OKAY;
}


static
SCIP_RETCODE setParentSetData(
   ParentSetData* psd,
   int num_items,
   int* num_parent_sets,
   int** num_parents_in_set,
   int*** parent_sets,
   char** names
   )
{
   psd->n = num_items;
   psd->nParentSets = num_parent_sets;
   psd->nParents = num_parents_in_set;
   psd->ParentSets = parent_sets;
   psd->nodeNames = names;
   return SCIP_OKAY;
}

static
SCIP_RETCODE setPedigreeData(
   SCIP* scip,
   int num_items,
   int* ages,
   char* sexes
   )
{
   PedigreeData peddata;
   peddata.n = num_items;
   peddata.ages = ages;
   peddata.sexes = sexes;
   peddata.SexVars = NULL;
   SCIP_CALL( MD_setPedigreeData(scip, &peddata) );
   return SCIP_OKAY;
}

static
SCIP_RETCODE setPropertyData(
   SCIP* scip,
   PropertyData* prop,
   int num_items,
   const char* filename,
   const char* frequency_file,
   int num_markers,
   SCIP_Bool allow_single_parents,
   SCIP_Bool have_ages,
   SCIP_Bool have_sexes,
   SCIP_Bool have_adults,
   int* ages,
   char* sexes,
   SCIP_Bool* adult,
   char*** genomes,
   int* num_alleles,
   SCIP_Real** allele_counts,
   char*** allele_names
   )
{
   int i;
   int j;

   prop->n = num_items;

   SCIP_CALL( PR_setGlobalProperty(scip, prop, "scorer_name", "GOBNILP") );
   SCIP_CALL( PR_setGlobalProperty(scip, prop, "scorer_version", GOBNILP_VERSION) );
   SCIP_CALL( PR_setGlobalProperty(scip, prop, "scorer_url", "www.cs.york.ac.uk/aig/sw/gobnilp/") );
   SCIP_CALL( PR_setGlobalProperty(scip, prop, "input_file", filename) );

   SCIP_CALL( PR_setGlobalProperty(scip, prop, "score_type", "Mendelian") );
   SCIP_CALL( PR_setGlobalPropertyFromInt(scip, prop, "num_markers", num_markers / 2) );
   SCIP_CALL( PR_setGlobalPropertyFromBool(scip, prop, "single_parents", allow_single_parents) );
   SCIP_CALL( PR_setGlobalPropertyFromBool(scip, prop, "filtered_by_age", have_ages) );
   SCIP_CALL( PR_setGlobalPropertyFromBool(scip, prop, "filtered_by_sex", have_sexes) );
   SCIP_CALL( PR_setGlobalPropertyFromBool(scip, prop, "filtered_by_sexual_maturity", have_adults) );

   if( strcmp(frequency_file, "") != 0 )
      for( i = 0; i < num_markers / 2; i++ )
      {
         char name_string[SCIP_MAXSTRLEN];
         char frequency_string[SCIP_MAXSTRLEN];
         sprintf(name_string, "marker_%d_alleles", i + 1);
         sprintf(frequency_string, "marker_%d_frequencies", i + 1);
         SCIP_CALL( PR_setGlobalPropertyFromArray(scip, prop, name_string, (const char**)allele_names[i], num_alleles[i]) );
         SCIP_CALL( PR_setGlobalPropertyFromRealArray(scip, prop, frequency_string, allele_counts[i], num_alleles[i]) );
      }

   for( i = 0; i < prop->n; i++ )
   {
      char genome_string[10000];
      if( have_ages )
         SCIP_CALL( PR_setPropertyFromInt(scip, prop, i, "age", ages[i]) );
      if( have_adults )
         SCIP_CALL( PR_setPropertyFromInt(scip, prop, i, "is_adult", adult[i]) );
      if( have_sexes )
      {
         if( sexes[i] == 'M' )
            SCIP_CALL( PR_setProperty(scip, prop, i, "sex", "Male") );
         if( sexes[i] == 'F' )
            SCIP_CALL( PR_setProperty(scip, prop, i, "sex", "Female") );
         if( sexes[i] == 'U' )
            SCIP_CALL( PR_setProperty(scip, prop, i, "sex", "Unknown") );
      }
      sprintf(genome_string, "%s,%s", genomes[i][0], genomes[i][1]);
      for( j = 1; j < num_markers / 2; j++ )
      {
         strcat(genome_string, " ");
         strcat(genome_string, genomes[i][2 * j]);
         strcat(genome_string, ",");
         strcat(genome_string, genomes[i][2 * j + 1]);
      }
      SCIP_CALL( PR_setProperty(scip, prop, i, "genome", genome_string) );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE storeResults(
   SCIP* scip,
   ParentSetData* psd,
   PropertyData* prop,
   SCIP_Real*** scores,
   int num_items,
   int* num_parent_sets,
   int** num_parents_in_set,
   int*** parent_sets,
   char** names,
   int* ages,
   char* sexes,
   SCIP_Bool* adult,
   char*** genomes,
   const char* filename,
   const char* frequency_file,
   int num_markers,
   SCIP_Bool allow_single_parents,
   SCIP_Bool have_ages,
   SCIP_Bool have_sexes,
   SCIP_Bool have_adults,
   int* num_alleles,
   SCIP_Real** allele_counts,
   char*** allele_names,
   SCIP_Real** likelihoods
   )
{
   SCIP_CALL( setParentSetData(psd, num_items, num_parent_sets, num_parents_in_set, parent_sets, names) );
   SCIP_CALL( setPedigreeData(scip, num_items, ages, sexes) );
   SCIP_CALL( setPropertyData(scip, prop, num_items, filename, frequency_file, num_markers, allow_single_parents, have_ages, have_sexes, have_adults, ages, sexes, adult, genomes, num_alleles, allele_counts, allele_names) );
   (*scores) = likelihoods;
   return SCIP_OKAY;
}

static
SCIP_RETCODE freeMemory(
   SCIP* scip,
   char**** genomes,
   char*** names,
   int** ages,
   char** sexes,
   SCIP_Bool** adult,
   int** num_alleles,
   SCIP_Real*** allele_counts,
   char**** allele_names,
   int** num_possible_parents,
   int*** possible_parents,
   int** num_parent_sets,
   int*** num_parents_in_set,
   int**** parent_sets,
   int** new_num_parent_sets,
   int*** new_num_parents_in_set,
   int**** new_parent_sets,
   SCIP_Real*** likelihoods,
   SCIP_Bool have_ages,
   SCIP_Bool have_sexes,
   SCIP_Bool have_adults,
   int first_row,
   int num_markers,
   int num_lines,
   char** unknown,
   SCIP_Bool*** allowed_parent_ptr,
   int** num_allowed_parent_ptr
   )
{
   int i;
   int j;

   for( i = 0; i < num_markers / 2; i++ )
   {
      for( j = 0; j < (*num_alleles)[i]; j++ )
         SCIPfreeMemoryArray(scip, &((*allele_names)[i][j]));
      SCIPfreeMemoryArray(scip, &((*allele_names)[i]));
      SCIPfreeMemoryArray(scip, &((*allele_counts)[i]));
   }
   SCIPfreeMemoryArray(scip, allele_counts);
   SCIPfreeMemoryArray(scip, allele_names);
   SCIPfreeMemoryArray(scip, num_alleles);

   for( i = 0; i < num_lines - first_row; i++ )
      SCIPfreeMemoryArray(scip, &((*possible_parents)[i]));
   SCIPfreeMemoryArray(scip, possible_parents);
   SCIPfreeMemoryArray(scip, num_possible_parents);

   for( i = 0; i < num_lines - first_row; i++ )
   {
      for( j = 0; j < (*num_parent_sets)[i]; j++ )
         SCIPfreeMemoryArray(scip, &((*parent_sets)[i][j]));
      SCIPfreeMemoryArray(scip, &((*parent_sets)[i]));
      SCIPfreeMemoryArray(scip, &((*num_parents_in_set)[i]));
   }
   SCIPfreeMemoryArray(scip, num_parent_sets);
   SCIPfreeMemoryArray(scip, num_parents_in_set);
   SCIPfreeMemoryArray(scip, parent_sets);


   for( i = 0; i < num_lines - first_row; i++ ) {
      for( j = 0; j < num_markers; j++ )
         free( (*genomes)[i][j] );
      SCIPfreeMemoryArray(scip, &((*genomes)[i]));
   }
   SCIPfreeMemoryArray(scip, genomes);

   if( have_sexes )
      SCIPfreeMemoryArray(scip, sexes);
   if( have_ages )
      SCIPfreeMemoryArray(scip, ages);
   if( have_adults )
      SCIPfreeMemoryArray(scip, adult);

   if( (*allowed_parent_ptr) != NULL )
   {
      for( i = 0; i < num_lines - first_row; ++i )
         if( (*allowed_parent_ptr)[i] != NULL )
            SCIPfreeMemoryArray(scip, &((*allowed_parent_ptr)[i])); 
      SCIPfreeMemoryArray(scip, allowed_parent_ptr);
   }
   
   if( (*num_allowed_parent_ptr) != NULL )
      SCIPfreeMemoryArray(scip, num_allowed_parent_ptr);

   return SCIP_OKAY;
}

/** Reads data in from a genome file and computes scores.
 *
 *  @param scip The SCIP instance being used.
 *  @param filename The file to read from.
 *  @param num_delims The number of field delimiters that are to be used when reading the file.
 *  @param delims The delimiters to use when reading the file.
 *  @param merge_delims Whether consecutive delimiters should be merged when reading the file.
 *  @param psd The parentage data that is read from the file.
 *  @param scores The score of each of the parent set combinations.
 *  @param frequency_file The file to read frequencies from or the empty string to calculate for the dataset.
 *  @return SCIP_OKAY if reading was successful, or an appropriate error code otherwise.
 */
SCIP_RETCODE PS_readProblemInGenomeFormat(
   SCIP* scip,
   const char* filename,
   int num_delims,
   char* delims,
   SCIP_Bool merge_delims,
   const char* frequency_file,
   ParentSetData* psd,
   SCIP_Real*** scores,
   PropertyData* prop
   )
{
   int first_row = 0;
   int num_markers;
   int num_lines;

   char***     genomes = NULL;
   char**      names   = NULL;
   int*        ages    = NULL;
   char*       sexes   = NULL;
   SCIP_Bool*  adult   = NULL;

   SCIP_Bool have_ages   = FALSE;
   SCIP_Bool have_sexes  = FALSE;
   SCIP_Bool have_adults = FALSE;

   int*          num_alleles   = NULL;
   SCIP_Real**   allele_counts = NULL;
   char***       allele_names  = NULL;

   int*    num_possible_parents = NULL;
   int**   possible_parents     = NULL;

   int*    num_parent_sets      = NULL;
   int**   num_parents_in_set   = NULL;
   int***  parent_sets          = NULL;

   int*    new_num_parent_sets      = NULL;
   int**   new_num_parents_in_set   = NULL;
   int***  new_parent_sets          = NULL;

   SCIP_Real** likelihoods = NULL;


   SCIP_Bool allow_single_parents;
   char* unknown;
   int maxmismatches;
   SCIP_Real mutation_probability;

   SCIP_Bool** allowed_parent = NULL;
   int* num_allowed_parent = NULL;
   char* pedigreeconstraintsfile;

   /* int i,j; */

   SCIPgetBoolParam(scip, "gobnilp/pedigree/singleparents", &allow_single_parents);
   SCIPgetStringParam(scip, "gobnilp/pedigree/realdata/unknown", &unknown);
   SCIPgetIntParam(scip, "gobnilp/pedigree/realdata/maxmismatches", &maxmismatches);
   SCIPgetRealParam(scip, "gobnilp/pedigree/realdata/mutationprob", &mutation_probability);
   SCIPgetStringParam(scip, "gobnilp/pedigree/constraintsfile", &pedigreeconstraintsfile);

   /* Get the raw data from the file */
   SCIP_CALL( readDataFromFile(scip, filename, delims, num_delims, merge_delims, &first_row, &num_lines, &num_markers, &names, &sexes, &ages, &genomes, &adult, &have_ages, &have_sexes, &have_adults) );

   /* Get allele frequency counts. */
   if( strequ(frequency_file, "") )
      SCIP_CALL( countFrequencies(scip, num_lines - first_row, num_markers, genomes, unknown, &num_alleles, &allele_counts, &allele_names) );
   else
      SCIP_CALL( readFrequencyFile(scip, frequency_file, num_delims, delims, merge_delims, num_markers, &num_alleles, &allele_counts, &allele_names) );

  /* Additional pedigree constraints */
   if( strcmp(pedigreeconstraintsfile, "") != 0 )
      SCIP_CALL( findAllowedParents(scip, num_lines - first_row, pedigreeconstraintsfile, names, &allowed_parent, &num_allowed_parent) );

   /* if( num_allowed_parent != NULL ) */
   /* { */
   /*    for( i = 0; i < num_lines - first_row; i++ ) */
   /*       printf("%d %d\n", i, num_allowed_parent[i]);  */
   /* } */
   /* else */
   /*    printf("is NULL\n"); */

   /* if( allowed_parent != NULL ) */
   /*    for( i = 0; i < num_lines - first_row; i++ ) */
   /*       if( allowed_parent[i] != NULL ) */
   /*          for( j = 0; j < num_lines - first_row; j++ ) */
   /*             printf("%d %d %u\n", i, j, allowed_parent[i][j]); */

   /* Find possible family triples. */
   SCIP_CALL( findPossibleParents(scip, num_lines - first_row, num_markers, genomes, ages, adult, maxmismatches, unknown, &num_possible_parents, &possible_parents, allowed_parent) );

   SCIP_CALL( findPossibleParentSets(scip, num_lines - first_row, num_possible_parents, possible_parents, allow_single_parents, have_sexes, sexes, &num_parent_sets, &num_parents_in_set, &parent_sets, num_allowed_parent) );

   /* Compute likelihoods. */
   SCIP_CALL( computeLikelihoods(num_lines - first_row, num_markers, genomes, maxmismatches, unknown, num_parent_sets, num_parents_in_set, parent_sets, num_alleles, allele_counts, allele_names, &likelihoods, &new_num_parent_sets, &new_num_parents_in_set, &new_parent_sets, mutation_probability) );

   /* Store the results. */
   SCIP_CALL( storeResults(scip, psd, prop, scores, num_lines-first_row, new_num_parent_sets, new_num_parents_in_set, new_parent_sets, names, ages, sexes, adult, genomes, filename, frequency_file, num_markers, allow_single_parents, have_ages, have_sexes, have_adults, num_alleles, allele_counts, allele_names, likelihoods) );

   /* Free memory. */
   SCIP_CALL( freeMemory(scip, &genomes, &names, &ages, &sexes, &adult, &num_alleles, &allele_counts, &allele_names, &num_possible_parents, &possible_parents, &num_parent_sets, &num_parents_in_set, &parent_sets, &new_num_parent_sets, &new_num_parents_in_set, &new_parent_sets, &likelihoods, have_ages, have_sexes, have_adults, first_row, num_markers, num_lines, &unknown, &allowed_parent, &num_allowed_parent) );

   return SCIP_OKAY;
}

