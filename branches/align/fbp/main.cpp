#include "defs.h"
#include "time_stamp.h"
#include "handlers.h"

// Used during error handling. See timestamp.cpp.
//extern const char *myname;


/* ----------------------------------------------------------------------
 * main()
 * ----------------------------------------------------------------------
 * Application main(). Pass in an XML filename as first arg.
 *
 * Returns 0 on success, -1 on error;
 * ----------------------------------------------------------------------
 */


int main(int argc, char *argv[]){
	char input_file_xml[MAX_FOLDER];
	xHMset *xset;
	allData *aD;
	bool tres;
#ifdef __linux__
#ifdef GNU_TRAP_FPE
    struct sigaction new_action;
    enable_fpe_traps ();
    //signal (SIGFPE, float_error);
    new_action.sa_sigaction=float_action;
    new_action.sa_flags=SA_SIGINFO;
    sigaction(SIGFPE,&new_action,NULL);
#endif
#endif /*__linux__*/

	// Get the app name for error reporting
	myname = strdup(basename(argv[0]));

	// Store the (global) starting time-of-day 
	init_timestamp();

	timestamp("starting at ",1);

	// Get command-line arg xml filename and read it
	if(argc < 2){
		printusage(argc, argv);
		return -1;
	}
	sprintf(input_file_xml,"%s",argv[argc-1]);
	printf("XML file:\n\"%s\"\n\n",input_file_xml);		
	xset = new xHMset();
	if(!xset->readXml(input_file_xml)){
		printf("Error:%s", xset->errorString);
		return -1;
	}

	// Main class to hold application data
	aD = new allData(xset);
	
	// Process the XML tree
	for(int i = 0; i < xset->nchild; i++){
		printf("======= Child structure #%i (of %i).\n\n",i+1,xset->nchild);
		switch(xset->child_name[i]){
			case NAME_FBP:
				printf(">>>> Filtered backprojection <<<<\n\n");

				timestamp("Calling xset->fbp",4);
				xset->fbp = (xFBP *)(xset->pointers[i]);
				timestamp("Returned from xset->fbp",4);

				timestamp("Calling main_fbp",4);
				tres = main_fbp(aD);
				timestamp("Returned from main_fbp",4);

				break;
			case NAME_XRM2TIF:
				printf(">>>> Conversion xRadia to TIFF + XML files <<<<\n\n");
				xset->x2t = (xXrm2tif *)(xset->pointers[i]);
				tres = main_xrm2tif(aD);
				break;
			case NAME_ALLFDK:
				printf(">>>> Cone beam, sub-volume reconstruction <<<<\n\n");
				xset->af = (xAllFDK *)(xset->pointers[i]);
				tres = main_allfdk(aD);
				break;
			default:
				printf(">>>> Unknown structure <<<<\n\n");
				break;
		}
		if(aD->fp != NULL){
			closeAllTags(aD,tres);
			fclose(aD->fp);
			aD->fp = NULL;
		}
		if(tres){
			printf("The application has finished its work successfully.\n\n");
		}else{
			printf("Error in the application.\n\n");
		}
		aD->data->freeMemory();
		printf("======= end of structure ===========\n\n");
	}

	delete xset;
	delete aD;
	return 0;
}

/* ----------------------------------------------------------------------
 * printusage()
 * ----------------------------------------------------------------------
 * Add command line usage and other help text here.
 * ----------------------------------------------------------------------
 */
void printusage(int argc, char **argv){
   printf("%s: Tomographic reconstruction program. Usage:\n", argv[0]);
   printf("%s  xmlfile\n", argv[0]);
}

