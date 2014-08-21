#include "defs.h"

bool printTag(allData *aD, char *tag, int i, char *comm){
	fprintf(aD->fp,"%s<%s info=\"%s\">%i</%s>\n", aD->cspace, tag, comm, i, tag);
	return true;
}

bool printTag(allData *aD, char *tag, int i){
	fprintf(aD->fp,"%s<%s>%i</%s>\n", aD->cspace, tag, i, tag);
	return true;
}

bool printTag(allData *aD, char *tag, unsigned int i, char *comm){
	fprintf(aD->fp,"%s<%s info=\"%s\">%u</%s>\n", aD->cspace, tag, comm, i, tag);
	return true;
}

bool printTag(allData *aD, char *tag, unsigned int i){
	fprintf(aD->fp,"%s<%s>%u</%s>\n", aD->cspace, tag, i, tag);
	return true;
}

bool printTag(allData *aD, char *tag, long i, char *comm){
	fprintf(aD->fp,"%s<%s info=\"%s\">%li</%s>\n", aD->cspace, tag, comm, i, tag);
	return true;
}

bool printTag(allData *aD, char *tag, long i){
	fprintf(aD->fp,"%s<%s>%li</%s>\n", aD->cspace, tag, i, tag);
	return true;
}


bool printTag(allData *aD, char *tag, float i, char *comm){
	fprintf(aD->fp,"%s<%s info=\"%s\">%f</%s>\n", aD->cspace, tag, comm, i, tag);
	return true;
}

bool printTag(allData *aD, char *tag, float i){
	fprintf(aD->fp,"%s<%s>%f</%s>\n", aD->cspace, tag, i, tag);
	return true;
}

bool printTag(allData *aD, char *tag, char *text, char *comm){
	fprintf(aD->fp,"%s<%s info=\"%s\">%s</%s>\n", aD->cspace, tag, comm, text, tag);
	return true;
}

bool printTag(allData *aD, char *tag, char *text){
	fprintf(aD->fp,"%s<%s>%s</%s>\n", aD->cspace, tag, text, tag);
	return true;
}

bool closeAllTags(allData *aD, bool tres){
	int nt = aD->tagNumber;
	for(int i=0; i<nt-1; i++){
		printTagEnd(aD);
	}
	print_global_time_stamp(aD);
	if(tres){
		printTag(aD,"Result","Success");
	}else{
		printTag(aD,"Result","Fail");
	}
	printTagEnd(aD);
	return true;
}


bool printTagStart(allData *aD, char *tag){
	fprintf(aD->fp,"%s<%s>\n",aD->cspace,tag);
	sprintf(aD->tagName + aD->tagNumber * MAX_TAG_LENGTH,"%s",tag);
	aD->tagNumber++;
	for(int i=0;i<aD->tagNumber;i++){
		aD->cspace[i] = '\t';
	}
	aD->cspace[aD->tagNumber]='\0';
	aD->time->local_start[aD->tagNumber-1] = clock();
	return true;
}

bool printTagStart2(allData *aD, char *tag){
	fprintf(aD->fp,"%s<%s>\n",aD->cspace,tag);
	sprintf(aD->tagName + aD->tagNumber * MAX_TAG_LENGTH,"%s",tag);
	aD->tagNumber++;
	for(int i=0;i<aD->tagNumber;i++){
		aD->cspace[i] = '\t';
	}
	aD->cspace[aD->tagNumber]='\0';
	return true;
}

bool printTagStart(allData *aD, char *tag, char *comm){
	fprintf(aD->fp,"%s<%s info=\"%s\">\n",aD->cspace, tag, comm);
	sprintf(aD->tagName + aD->tagNumber * MAX_TAG_LENGTH,"%s",tag);
	aD->tagNumber++;
	for(int i=0;i<aD->tagNumber;i++){
		aD->cspace[i] = '\t';
	}
	aD->cspace[aD->tagNumber]='\0';
	aD->time->local_start[aD->tagNumber-1] = clock();
	return true;
}


bool printTagStart2(allData *aD, char *tag, char *comm){
	fprintf(aD->fp,"%s<%s info=\"%s\">\n", aD->cspace, tag, comm);
	sprintf(aD->tagName + aD->tagNumber * MAX_TAG_LENGTH,"%s",tag);
	aD->tagNumber++;
	for(int i=0;i<aD->tagNumber;i++){
		aD->cspace[i] = '\t';
	}
	aD->cspace[aD->tagNumber]='\0';
	return true;
}


bool printTagEnd(allData *aD){
	print_local_time_stamp(aD);
	aD->tagNumber--;
	for(int i=0;i<aD->tagNumber;i++){
		aD->cspace[i] = '\t';
	}
	aD->cspace[aD->tagNumber]='\0';
	fprintf(aD->fp,"%s</%s>\n", aD->cspace, aD->tagName + aD->tagNumber * MAX_TAG_LENGTH);
	return true;
}


bool printTagEnd2(allData *aD){
	aD->tagNumber--;
	for(int i=0;i<aD->tagNumber;i++){
		aD->cspace[i] = '\t';
	}
	aD->cspace[aD->tagNumber]='\0';
	fprintf(aD->fp,"%s</%s>\n", aD->cspace, aD->tagName + aD->tagNumber * MAX_TAG_LENGTH);
	return true;
}


bool printInfo(allData *aD){
	fprintf(aD->fp,"%s<Info>%s</Info>\n",aD->cspace, aD->message);
	return true;
}

bool printError(allData *aD){
	printf("%s<Error>%s</Error>\n",aD->cspace, aD->message);
	fprintf(aD->fp,"%s<Error>%s</Error>\n",aD->cspace, aD->message);
	return true;
}

bool printWarning(allData *aD){
	fprintf(aD->fp,"%s<Warning>%s</Warning>\n",aD->cspace, aD->message);
	return true;
}

bool printInfo(allData *aD, char *text){
	fprintf(aD->fp,"%s<Info>%s</Info>\n",aD->cspace, text);
	return true;
}

bool printError(allData *aD, char *text){
	printf("%s<Error>%s</Error>\n",aD->cspace, text);
	fprintf(aD->fp,"%s<Error>%s</Error>\n",aD->cspace, text);
	return true;
}

bool printWarning(allData *aD, char *text){
	fprintf(aD->fp,"%s<Warning>%s</Warning>\n",aD->cspace, text);
	return true;
}


bool print_global_time_stamp(allData *aD){
	char sf[100], strDest[200];

#if defined(_WIN32) || defined(_WIN64)
	__time64_t long_time;
	_time64( &long_time );
	tm *newtime = _localtime64(&long_time); 
#else
	time_t long_time;
	time( &long_time );
	struct tm *newtime = localtime( &long_time );
#endif
	printTagStart2(aD,"TimeStamp");
	
	sprintf(sf,"%%d %%B %%Y");
	strftime(strDest,200,sf,newtime);
	printTag(aD,"Date",strDest);

	sprintf(sf,"%%A");
	strftime(strDest,200,sf,newtime);
	printTag(aD,"Weekday",strDest);

	sprintf(sf,"%%H:%%M:%%S");
	strftime(strDest,200,sf,newtime);
	printTag(aD,"Time",strDest);

	sprintf(sf,"%%Z");
	strftime(strDest,200,sf,newtime);
	printTag(aD,"TimeZone",strDest);

	printTagEnd2(aD);	   
	return true;
}



bool print_local_time_stamp(allData *aD){
	aD->time->local_end = clock();
	double duration;
	duration = (double)(aD->time->local_end - aD->time->local_start[aD->tagNumber-1])/CLOCKS_PER_SEC;
	printTag(aD,"Duration",(float)duration);
	return true;
}

