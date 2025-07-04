/*
*
*  msi2lmp.exe  V3.1
*
*   v3.1 JEC - changed IO interface to standard in/out, forcefield file
*              location can be indicated by environmental variable; added 
*              printing options, consistency checks and forcefield
*              parameter versions sensitivity (highest one used)
*
*   v3.0 JEC - program substantially rewritten to reduce execution time
*              and be 98 % dynamic in memory use (still fixed limits on
*              number of parameter types for different internal coordinate
*              sets)
*
*   v2.0 MDP - got internal coordinate information from mdf file and
*              forcefield parameters from frc file thus eliminating 
*              need for Discover
*
*   V1.0 SL  - original version. Used .car file and internal coordinate
*              information from Discover to produce LAMMPS data file.
*
*  This program uses the .car and .mdf files from MSI/Biosyms's INSIGHT
*  program to produce a LAMMPS data file.
*
*  The program is started by supplying information at the command prompt
* according to the usage described below.  
*
*  USAGE: msi2lmp3 ROOTNAME {-print #} {-class #} {-frc FRC_FILE}
*
*  -- msi2lmp3 is the name of the executable
*  -- ROOTNAME is the base name of the .car and .mdf files
*
*  -- -print
*        # is the print level:  0  - silent except for errors
*                               1  - minimal (default)
*                               2  - more verbose 
*  -- -class 
*        # is the class of forcefield to use (I  = Class I e.g., CVFF)
*		  			     (II = Class II e.g., CFFx )
*     default is -class I
*
*  -- -frc   - specifies name of the forcefield file (e.g., cff91)
* 
*     If the name includes a hard wired directory (i.e., if the name
*     starts with . or /), then the name is used alone. Otherwise,
*     the program looks for the forcefield file in $BIOSYM_LIBRARY.
*     If $BIOSYM_LIBRARY is not set, then the current directory is 
*     used.
*
*     If the file name does not include a dot after the first
*     character, then .frc is appended to the name.
*
*     For example,  -frc cvff (assumes cvff.frc is in $BIOSYM_LIBRARY
*                              or .)
*
*                   -frc cff/cff91 (assumes cff91.frc is in 
*                                   $BIOSYM_LIBRARY/cff or ./cff)
*
*                   -frc /usr/local/biosym/forcefields/cff95 (absolute
*                                                             location)
*
*     By default, the program uses $BIOSYM_LIBRARY/cvff.frc
*
*  -- output is written to a file called data.ROOTNAME
*
*
****************************************************************
*
* Msi2lmp3
*
* This is the third version of a program that generates a LAMMPS
* data file based on the information in a MSI car file (atom
* coordinates) and mdf file (molecular topology). A key part of
* the program looks up forcefield parameters from an MSI frc file.
*
* The first version was written by Steve Lustig at Dupont, but
* required using Discover to derive internal coordinates and
* forcefield parameters
*
* The second version was written by Michael Peachey while an
* in intern in the Cray Chemistry Applications Group managed
* by John Carpenter. This version derived internal coordinates
* from the mdf file and looked up parameters in the frc file
* thus eliminating the need for Discover.
*
* The third version was written by John Carpenter to optimize
* the performance of the program for large molecular systems
* (the original  code for deriving atom numbers was quadratic in time)
* and to make the program fully dynamic. The second version used
* fixed dimension arrays for the internal coordinates.
*
* John Carpenter can be contacted by sending email to
* jec374@earthlink.net
*
* November 2000
*/

#define MAIN

#include "Msi2LMP2.h"

int main (int argc, char *argv[])
{
   int n,i,found_dot;		/* Counter */
   char *string;
   char *frc_dir_name;
   char *frc_file_name;
   FILE *DatF;

 /* Functions called from within main */

/* All code is located in .c file with function name */
   extern void FrcMenu();
   extern void ReadCarFile();
   extern void ReadMdfFile();
   extern void ReadFrcFile();
   extern void MakeLists();
   extern void GetParameters(int);
   extern void CheckLists();
   extern void WriteDataFile(FILE *,char *,int);



   pflag = 1;
   forcefield = 1;		/* Variable that identifies forcefield to use */

   frc_file_name = (char *) calloc(160,sizeof(char));
   frc_dir_name = (char *) calloc(160,sizeof(char));

   frc_dir_name = getenv("BIOSYM_LIBRARY");

   if (frc_dir_name == NULL) {
     frc_file_name = strcpy(frc_file_name,"./cvff.frc");
   }
   else {
     for (i=0; i < strlen(frc_dir_name); i++)
       frc_file_name[i] = frc_dir_name[i];
     frc_file_name = strcat(frc_file_name,"/cvff.frc");
   }

   if (argc < 2) { /* If no rootname was supplied, prompt for it */
     fprintf(stderr,"The rootname of the .car and .mdf files must be entered\n");
   }
   else /* rootname was supplied as first argument, copy to rootname */
     sprintf(rootname,"%s",argv[1]);

   n = 2;
   while (n < argc) {
     if (strcmp(argv[n],"-class") == 0) {
       if (strcmp(argv[n+1],"I") == 0) {
	 forcefield = 1;
	 n++;
       }
       else if (strcmp(argv[n+1],"II") == 0) {
	 forcefield = 2;
	 n++;
       }
       else {
	 fprintf(stderr,"Unrecognized Forcefield class: %s\n",
		argv[n+1]);
	 n++;
       }
     }
     else if (strcmp(argv[n],"-frc") == 0) {
       if ((frc_dir_name != NULL) && (argv[n+1][0] != '.')) {
	 for (i=0; i < strlen(frc_dir_name); i++) {
	   frc_file_name[i] = frc_dir_name[i];
	 }
	 frc_file_name[strlen(frc_dir_name)] = '\0';
	 frc_file_name = strcat(frc_file_name,"/");
	 frc_file_name = strcat(frc_file_name,argv[n+1]);
       }
       else {
	 frc_file_name = strcpy(frc_file_name,argv[n+1]);
       }
       found_dot = 0;
       for (i=1; i < strlen(frc_file_name); i++) {
	 if (frc_file_name[i] == '.') found_dot = 1;
       }
       if (found_dot == 0) 
	 frc_file_name = strcat(frc_file_name,".frc");
       n++;
     }
     else if (strstr(argv[n],"-p") != NULL) {
       pflag = atoi(argv[n+1]);
       n++;
     }
     else {
       fprintf(stderr,"Unrecognized option: %s\n",argv[n]);
     }
     n++;
   }
   for (i=0; i < strlen(frc_file_name); i++) 
     FrcFileName[i] = frc_file_name[i];
   free(frc_file_name);

   if (pflag > 0) {
     fprintf(stderr,"\nRunning Msi2lmp.....\n\n");
     fprintf(stderr," Forcefield file name: %s\n",FrcFileName);
     fprintf(stderr," Forcefield class: %d\n\n",forcefield);
   }

   if (((forcefield == 1) && (strstr(FrcFileName,"cff") != NULL) ||
	(forcefield == 2) && (strstr(FrcFileName,"cvff") != NULL))) {
     fprintf(stderr," WARNING - forcefield name and class appear to\n");
     fprintf(stderr,"           be inconsistent - Errors may result\n\n");
   }

 /* Read in .car file */

   ReadCarFile();

 /*Read in .mdf file */

   ReadMdfFile();

 /* Define bonds, angles, etc...*/

   if (pflag > 0) fprintf(stderr,"\n Building internal coordinate lists \n");
   MakeLists();

 /* Read .frc file into memory */

   if (pflag > 0) fprintf(stderr,"\n Reading forcefield file \n");
   ReadFrcFile();

 /* Get forcefield parameters */

   if (pflag > 0) fprintf(stderr,"\n Get parameters for this molecular system\n");
   GetParameters(forcefield);

  /* Do internal check of internal coordinate lists */

   if (pflag > 0) fprintf(stderr,"\n Check parameters for internal consistency\n");
   CheckLists();

   if (pflag > 0) fprintf(stderr,"\n Writing LAMMPS data file\n");
   WriteDataFile(stdout,rootname,forcefield);

   if (pflag > 0) fprintf(stderr,"\nNormal program termination\n");
}
