/* ----------------------------------------------------------------------
 * xmlStructure.cpp()
 * ----------------------------------------------------------------------
 * Read the XML file given reconstruction parameters and store the
 * various params in the general storage class. The XML structure is
 * hard-coded in here.
 *
 * TODO: Add interface to external XML library
 * ----------------------------------------------------------------------
 */
#include "defs.h"


bool getValue(Ipp8u *vf, int ch2, int ch3, char *value){
	int len = ch3-ch2-1;
	if(len>0)ippsCopy_8u(vf+ch2+1,(Ipp8u *)value,len);
	value[len]='\0';
//	printf("Value: \"%s\".\n",value);
	return true;
}


/* ----------------------------------------------------------------------
 * removeContiguousSpace()
 * ----------------------------------------------------------------------
 * Removes in-place a contiguous block of white-space starting from the
 * given position in the string. The incoming length of the string is
 * modified to hold the new length of the string. If the given start
 * position does not represent a white-space character nothing is done
 * to the string.
 *
 * Args:
 * In/Out
 *	Ipp8u	*str:	String to modify
 *	int	*len:	Length of current and modified string
 * In
 *	int	start_pos: index in to str from where we start
 *
 * Returns number of white space chars removed
 * ----------------------------------------------------------------------
 */

int removeContiguousSpace(Ipp8u *str, int *len, int start_pos){
  int space_count = 0;

  if ( start_pos >= *len )
    return 0;

  // Step over white-space until we hit a non white-space char
  for(int i=start_pos; i<*len; i++){
    if( !isspace(str[i]) ) 
      break;
    ++space_count;
  }

  // Remove that white-space from the string
  ippsRemove_8u_I(str, len, start_pos, space_count);

  // Return number of chars removed
  return space_count;
}

// Reverse direction of above
int removeContiguousSpaceRev(Ipp8u *str, int *len, int start_pos){
  int space_count = 0;

  if ( start_pos >= *len )
    return 0;

  // Step backwards over white-space until we hit a non white-space char
  for(int i=start_pos; i>=0; i--){
    if( !isspace(str[i]) ) 
      break;
    ++space_count;
  }

  // Remove that white-space from the string
  ippsRemove_8u_I(str, len, start_pos-space_count+1, space_count);

  // Return number of chars removed
  return space_count;
}


/* ----------------------------------------------------------------------
 * findNextTag()
 * ----------------------------------------------------------------------
 * Search for a string between < and > and return a copy of that string.
 * Note that we are NOT searching for a specific tag name. We are simply
 * getting the string found between the next < and > characters.
 *
 * The last parameter is optional and will default to a NULL pointer.
 * If a valid pointer is passed in a test will be performed to see if the
 * tag is of the form: <tag/>. This is an alternative way of specifying an
 * empty tag (rather than doing <tag></tag>). The caller can use this to
 * determine whether to carry on looking for the end tag, for example. 
 *
 * Args:
 * In
 *	Ipp8u	*vf:		Entire XML file flattened to an array
 *	int	start_pos:	Current starting position in array
 *	int	end_pos:	Position of last point up to which we will search
 * Out
 *	int	*tag_left_pos:	Position of the '<' in the vf array
 *	int	*tag_right_pos:	Position of the '>' in the vf array
 *	char	*tagName:	Will get a copy of the string between the < >
 *	bool	*isEmpty:	If a valid ptr, set to true if tag is an empty
 *				tag of the form <tag/>, otherwise set false;
 *			  
 *
 * Returns true if *any* tag name was found, false on error.
 * ----------------------------------------------------------------------
 */
bool findNextTag(Ipp8u *vf, int start_pos, int end_pos, int *tag_left_pos, int *tag_right_pos, char *tagName, bool *isEmpty)
{
	Ipp8u str_tag_left[2]  = "<";
	Ipp8u str_tag_right[2] = ">";
	int i1;			// Position of found '<' (within entire vf string)
	int i2;			// Position of found '>' (within entire vf string)
	int offset;		// Offset of found '<' or '>' from starting point
	int max_search_len;	// Number of chars of vf that we can search through
	int tag_name_len;	// Length of string between '<' and '>'

	// Only search between the given start and end positions
	max_search_len = end_pos-start_pos+1;

	// Search for '<' and get its position
	ippsFind_8u(vf+start_pos, max_search_len, str_tag_left, 1, &offset);
	if( offset<0 )
	  return false;
	i1 = start_pos+offset;

	// Search for '>' (starting from '<') and get its position
	ippsFind_8u(vf+i1, max_search_len, str_tag_right, 1, &offset);
	if( offset<0 )
	  return false;
	i2 = i1+offset;

	tag_name_len = offset-1;

	// Copy the tag name and overwrite the ending '>'
	ippsCopy_8u(vf+i1+1,(Ipp8u *)tagName, tag_name_len);
	tagName[tag_name_len]='\0';

	// Check if last char of tagName is '/' meaning the tag was
	// an empty tag of the form <tagName/> and so doesn't have
	// an end tag of the form <tagName></tagName>. Only do this
	// test if a non-null ptr is passed in.
	if ( isEmpty != 0 ) {
	  if ( tagName[tag_name_len-1] == '/' ) {

	    // Overwrite the '/'
	    tagName[tag_name_len-1] = '\0';
	    *isEmpty=1;
	    // printf( "Empty tag: %s\n", tagName );
	  }
	  else
	    *isEmpty=0;
	}

	// Return the positions of the < and > characters
	*tag_left_pos  = i1;
	*tag_right_pos = i2;

	return true;
}

/* ----------------------------------------------------------------------
 * findEndTag()
 * ----------------------------------------------------------------------
 * Search for a tag of the form "</tagName>" where tagName is known.
 * Unlike findNextTag(), we now know the tag name to search for so can
 * search for that specific name.
 *
 * We do not handle empty tags of the form <tagName/>. Those are handled
 * by the findNextTag() function when looking for opening tags.
 *
 * Args:
 * In
 *	Ipp8u	*vf:		Entire XML file flattened to an array
 *	int	start_pos:	Current starting position in array
 *	int	end_pos:	Position of last point up to which we will search
 * Out
 *	int	*tag_left_pos:	Position of the '<' in the vf array
 *	int	*tag_right_pos:	Position of the '>' in the vf array
 *	char	*tagName:	Will get a copy of the string between the < >
 *	bool	*isEmpty:	If a valid ptr, set to true if tag is an empty
 *				tag of the form <tag/>, otherwise set false;
 *			  
 *
 * Returns true if the closing tag name was found, false on error.
 * ----------------------------------------------------------------------
 */

bool findEndTag(Ipp8u *vf, int start_pos, int end_pos, int *tag_left_pos, int *tag_right_pos, char *tagName)
{
	Ipp8u etag[100];	// Complete tag (with </ > added) to search for
	int etag_len;		// Length of tag (with </ > added).
	int offset;		// Offset of found '<' or '>' from starting point
	int max_search_len;	// Number of chars of vf that we can search through

	// Construct end tag with given name
	sprintf((char*)etag, "</%s>", tagName);
	etag_len = (int)strlen((char *)etag);

	// Only search between the given start and end positions
	max_search_len = end_pos-start_pos+1;

	// Search for tag from given starting point up to given end point.
	ippsFind_8u(vf+start_pos, max_search_len, etag, etag_len, &offset);
	if( offset<0 )
	  return false;

	// Return the positions of the < and > characters
	*tag_left_pos  = start_pos+offset;
	*tag_right_pos = *tag_left_pos+etag_len-1;
	
	return true;
}

xBackprojection::xBackprojection(){
	filter = new xFilter();
	tilt = new xTilt();
	coordinateSystem = new xCoordinateSystem();
	circles = new xCircles();
	roi = new xRoi();
	pc_interpolation = PC_INTERPOLATION_LINEAR;
}

xBackprojection::~xBackprojection(){
	delete filter; filter = NULL;
	delete tilt; tilt = NULL;
	delete coordinateSystem; coordinateSystem = NULL;
	delete circles; circles = NULL;
	delete roi; roi = NULL;
}
bool xBackprojection::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"Backprojection");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
				
		
		if(strcmp(childTag,"Filter") == 0){
			if(!this->filter->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"Tilt") == 0){
			//if(!this->flatField->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"CoordinateSystem") == 0){
			//if(!this->flatField->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"Circles") == 0){
			if(!this->circles->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		
		if(strcmp(childTag,"ROI") == 0){
			if(!this->roi->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		

		if(strcmp(childTag,"ImageCentre") == 0){
		  getValue(vf,ch2,ch3,value); 
		  this->imageCentre = float(atof(value)); 
		  izt = true;
		}

		if(strcmp(childTag,"ClockwiseRotation") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				pc_interpolation = CLOCKWISE_ROTATION_NO;
			}else if(strcmp(value,"yes") == 0){
				pc_interpolation = CLOCKWISE_ROTATION_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\" or \"Yes\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"PolarCartesianInterpolation") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"nearestneighbour") == 0){
				pc_interpolation = PC_INTERPOLATION_NN;
			}else if(strcmp(value,"linear") == 0){
				pc_interpolation = PC_INTERPOLATION_LINEAR;
			}else if(strcmp(value,"cubic") == 0){
				pc_interpolation = PC_INTERPOLATION_CUBIC;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"NearestNeighbour\", \"Linear\" or \"Cubic\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	
	return true;
}





xCircles::xCircles(){
	radiusMin = new xRadius();
	radiusMax = new xRadius();
	radiusStep = new xRadius();
}

xCircles::~xCircles(){
	delete radiusMin; radiusMin = NULL;
	delete radiusMax; radiusMax = NULL;
	delete radiusStep; radiusStep = NULL;
}


bool xCircles::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100];
	sprintf(mTag,"Circles");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}

	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
		
		if(strcmp(childTag,"ValueMin") == 0){
			if(!this->radiusMin->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"ValueMax") == 0){
			if(!this->radiusMax->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"ValueStep") == 0){
			if(!this->radiusStep->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}

	return true;
}



xFBP::xFBP(){
	GPUDeviceNumber = 0;
	beamlineUser = new xBeamlineUser();
	inputData = new xInputData();
	flatDarkFields = new xFlatDarkFields();
	preprocessing = new xPreprocessing();
	transform = new xTransform();
	backprojection = new xBackprojection();
	outputData = new xOutputData();
}

xFBP::~xFBP(){
	delete beamlineUser; beamlineUser = NULL;
	delete inputData; inputData = NULL;
	delete flatDarkFields; flatDarkFields = NULL;
	delete preprocessing; preprocessing = NULL;
	delete transform; transform = NULL;
	delete backprojection; backprojection = NULL;
	delete outputData; outputData = NULL;

}

bool xFBP::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"FBP");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}

	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;

		if(strcmp(childTag,"GPUDeviceNumber") == 0){getValue(vf,ch2,ch3,value); this->GPUDeviceNumber = atoi(value); izt = true;}

		if(strcmp(childTag,"DefaultXml") == 0){
			getValue(vf,ch2,ch3,value);
			sprintf(this->defaultXml,value);
			izt = true;
		}
		if(strcmp(childTag,"LogFile") == 0){
			getValue(vf,ch2,ch3,value);
			replaceSlash(value);
			sprintf(this->logFile,value);
			izt = true;
		}
		if(strcmp(childTag,"BeamlineUser") == 0){
			izt = true;
		}
		if(strcmp(childTag,"InputData") == 0){
			if(!inputData->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		if(strcmp(childTag,"FlatDarkFields") == 0){
			if(!this->flatDarkFields->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		if(strcmp(childTag,"Preprocessing") == 0){
			if(!this->preprocessing->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		if(strcmp(childTag,"Transform") == 0){
			if(!this->transform->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		if(strcmp(childTag,"Backprojection") == 0){
			if(!this->backprojection->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		if(strcmp(childTag,"OutputData") == 0){
			if(!this->outputData->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}

	
	return true;
}



xFlatDarkFields::xFlatDarkFields(){
	this->flatField = new xFDField();
	this->darkField = new xFDField();
}
xFlatDarkFields::~xFlatDarkFields(){
	delete flatField; flatField = NULL;
	delete darkField; darkField = NULL;
}

bool xFlatDarkFields::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100];
	sprintf(mTag,"FlatDarkFields");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}

	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
		
		if(strcmp(childTag,"FlatField") == 0){
			if(!this->flatField->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"DarkField") == 0){
			if(!this->darkField->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}

	return true;
}






xHMset::xHMset(){
	pointers = (void **)malloc(20*sizeof(void *));
	outState = true;
	errorString[0] = '\0';
	fbp = NULL;
	#ifdef _HMXIF
	x2t = NULL;
	af = NULL;
	#endif // _HMXIF
	//this->fbp = new xFBP();
	//this->x2t = new xXrm2tif();
	lenf = 0;
	vf = NULL;
	nchild = 0;
}

xHMset::~xHMset(){
	for(int i=0;i<nchild;i++){
		switch(child_name[i]){
			case NAME_FBP:
				delete (xFBP *)(pointers[i]);
				break;
			#ifdef _HMXIF
			case NAME_XRM2TIF:
				delete (xXrm2tif *)(pointers[i]);
				break;
			case NAME_ALLFDK:
				delete (xAllFDK *)(pointers[i]);
				break;
			#endif // _HMXIF
			default:
				break;
		}
	}
	free(pointers); pointers = NULL;
	fbp = NULL;
	#ifdef _HMXIF
	x2t = NULL;
	af = NULL;
	#endif // _HMXIF
	if(vf != NULL){
		ippsFree(vf); vf = NULL;
	}
}


bool xHMset::readXml(char *inputFile){
	char tagName[100], childTag[100];
	int pos1, pos2, pos;
	int s1, s2, s3, s4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;	
	Ipp8u str_comm_left[5], str_comm_right[5];
	Ipp8u str_xml_left[3], str_xml_right[3];
	Ipp8u str_tag_left[2], str_tag_right[2];
	Ipp8u str_space[2], str_any_space[5];
	FILE *ff;

	sprintf((char *)str_comm_left,"<!--");
	sprintf((char *)str_comm_right,"-->");
	sprintf((char *)str_tag_left,"<");
	sprintf((char *)str_tag_right,">");
	sprintf((char *)str_xml_left,"<?");
	sprintf((char *)str_xml_right,"?>");
	sprintf((char *)str_space," ");
	sprintf((char *)str_any_space, " \t\n\r");

	// Read XML file in to array
	ff = fopen(inputFile,"r");
	if(ff == NULL){
		sprintf(errorString,"Can not open file \"%s\".\n",inputFile);
		return false;
	}
	fseek(ff,0,SEEK_END);
	lenf = ftell(ff);
	fseek(ff,0,SEEK_SET);
	vf = ippsMalloc_8u(lenf);
	if(vf == NULL){
		printf("Can not allocate memory to read xml file.");
		return false;
	}
	fread(vf,sizeof(Ipp8u),lenf,ff);
	fclose(ff); ff = NULL;

	// GWL: Why do the following?
	ippsAddC_8u_ISfs(128, vf, lenf, 0);
	ippsSubC_8u_ISfs(128, vf, lenf, 0);
	ippsSubC_8u_ISfs(32, vf, lenf, 0);
	ippsAddC_8u_ISfs(32, vf, lenf, 0);

	// ippsReplaceC_8u() crashes on icpc v11.0.074
#if 0	
	vt = ippsMalloc_8u(lenf);
	ippsReplaceC_8u(vf, vt, lenf, 127, 32);
	ippsCopy_8u(vt,vf,lenf);
	ippsFree(vt); vt = NULL;
#endif

	// Strip various unwanted items from the XML
	// -----------------------------------------
	
	// Remove the header: <? blah ?>
	for(;;){
		// Search for a '<?'
		ippsFind_8u(vf,lenf,str_xml_left, 2, &pos1);
		if(pos1<0)
		  break;

		// Search for a '?>'
		ippsFind_8u(vf,lenf,str_xml_right, 2, &pos2);
		if(pos2<pos1){
			printf("Error (xml file, header)\n");
			return false;
		}

		// Remove the entire <? ... ?> string
		ippsRemove_8u_I(vf, &lenf, pos1, pos2-pos1+2);
	}

	// Remove any comments: <-- blah -->
	for(;;){
		// Search for a '<--'
		ippsFind_8u(vf,lenf, str_comm_left, 3, &pos1);
		if(pos1<0)
		  break;

		// Search for a '-->'
		ippsFind_8u(vf,lenf, str_comm_right, 3, &pos2);
		if(pos2<pos1){
			printf("Error (xml file, comments)\n");
			return false;
		}

		// Remove entire <-- ... --> string
		ippsRemove_8u_I(vf, &lenf, pos1, pos2-pos1+3);
	}

	// I think this is here to remove any attributes within a tag: <tag attr="prop">.
	// It removes *all* chars starting at the first space after a tag's '<' to the tag's end '>'.
	// If no attributes, it will remove any unwanted space after the tag name: <tag   > 
	// but fails if space exists before tag name <   tag>. 
	// Now also supports tags of the form <tag attr="prop"/>.
	pos = 0;
	for(;;){
		// Search for a '<'
		ippsFind_8u(vf+pos, lenf-pos, str_tag_left, 1, &pos2);
		if(pos2 <0)
		  break;
		pos+=pos2;	// Jump to the '<'. All searches now start from here.

		// Search for a '>' but step back if '/>' closes the tag
		ippsFind_8u(vf+pos, lenf-pos, str_tag_right, 1, &pos1);
		if ( pos1 > 0 && vf[pos+pos1-1] == '/' ) {
		  --pos1;		// Use the '/' as the end marker
		}
		
		// Search for a ' ' which marks start of attributes after tag name
		ippsFind_8u(vf+pos, lenf-pos, str_space, 1, &pos2);

		// Would like to search for any of the white-space in str_any_space
		// but this crashes (in icpc v11.0.074)
		// ippsFindCAny_8u(vf+pos, lenf-pos, str_any_space, 4, &pos2);

		// Error or no space in <tag> so gone beyond the tag's ending '/' or '>'
		if(pos2 < 0 || pos2 > pos1){
			// Jump to '/' or '>' and repeat
			pos += pos1;
			continue;
		}

		// Remove char from first white-space to '/' or '>'
		ippsRemove_8u_I(vf, &lenf, pos+pos2, pos1-pos2);

		// Step over the initial '>' for next search
		++pos;

	}


	// Remove space from very start of XML to '<' of first tag
	if (isspace(vf[0])){
		ippsFind_8u(vf,lenf,str_tag_left, 1, &pos1);
		ippsRemove_8u_I(vf, &lenf, 0, pos1);
	}

	// Remove space from '>' of last tag to end of string
	if (isspace(vf[lenf-1])) {
		ippsFindRev_8u(vf,lenf,str_tag_right, 1, &pos1);
		ippsRemove_8u_I(vf, &lenf, pos1+1, lenf-pos1-1);
	}

	// Remove spaces after a '>': ...blah>   blah
	pos = 0;
	do{
		// Search for a '>' 
		ippsFind_8u(vf+pos,lenf,str_tag_right, 1, &pos1);
		if( pos1<0 )
		  break;
		pos += pos1+1;	// Jump to char after the '>'

		// If there's a block of white space, remove it.
		removeContiguousSpace(vf, &lenf, pos );

	}while(pos<lenf);

	// Remove spaces before a '<': blah    <blah...
	pos = lenf-1;		// Start at end of string
	do{
		ippsFindRev_8u(vf,pos,str_tag_left, 1, &pos1);
		if(pos1<0)break;
		pos=pos1-1;	// Jump to the char before the  '<'

		// If there's a block of white space before us, remove it
		pos -= removeContiguousSpaceRev( vf, &lenf, pos );
	}while(pos>=0);

	// We've done a lot of inplace removal so if you want to print the
	// final string then hide the junk left over.
	vf[lenf] = '\0';

	// Start parsing tag names
	// -----------------------

	// printf( "\nWorking with the following:\n%sEND OF INPUT STRING\n", vf );
	
	t = findNextTag(vf, 0, lenf-1, &s1, &s2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the root start-tag.\n");
		return false;
	}
	
	t = findEndTag(vf, s2, lenf-1, &s3, &s4, tagName);
	if(!t){
		printf("Can not find the root end-tag.\n");
		return false;
	}

	if(strcmp(tagName,"HMxml") != 0){
		printf("Wrong root tag: \"%s\" (should be \"HMxml\").\n",tagName);
		return false;
	}

	for(;;){
	  t = findNextTag(vf, s2, s3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  s2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, s3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
		if(strcmp(childTag,"FBP") == 0){
			pointers[nchild] = (void *) (new xFBP());
			this->child_name[nchild] = NAME_FBP;
			if(!((xFBP *)pointers[nchild])->readXml(vf,ch1,ch4))return false;
			
			izt = true;
		}
		#ifdef _HMXIF
		if(strcmp(childTag,"Xrm2tif") == 0){
			pointers[nchild] = (void *) (new xXrm2tif());
			this->child_name[nchild] = NAME_XRM2TIF;
			if(!((xXrm2tif *)pointers[nchild])->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		if(strcmp(childTag,"AllFDK") == 0){
			pointers[nchild] = (void *) (new xAllFDK());
			this->child_name[nchild] = NAME_ALLFDK;
			if(!((xAllFDK *)pointers[nchild])->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		#endif // _HMXIF
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"HMxml\".\n",childTag);
			return false;
		}
		nchild++;

		s2 = ch4+1;

	}
	return true;
}


xInputData::xInputData(){
	this->raw = new xRaw();
	folder[0] = '\0';
	prefix[0] = '\0';
	suffix[0] = '\0';
	extension[0] = '\0';
}

xInputData::~xInputData(){
	delete raw; raw = NULL;
}

bool xInputData::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"InputData");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag=0;

	// Search between s1 and s2 for any tag name.
	// t1 will be position of the '<', t2 position of '>'
	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}

	
	for(;;){
	  t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;

			
		if(strcmp(childTag,"Folder") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->folder,value); izt = true;}
		if(strcmp(childTag,"Prefix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->prefix,value); izt = true;}
		if(strcmp(childTag,"Suffix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->suffix,value); izt = true;}
		if(strcmp(childTag,"Extension") == 0){getValue(vf,ch2,ch3,value); sprintf(this->extension,value); izt = true;}
		if(strcmp(childTag,"NOD") == 0){getValue(vf,ch2,ch3,value); this->NOD = atoi(value); izt = true;}
		if(strcmp(childTag,"FileFirst") == 0){getValue(vf,ch2,ch3,value); this->fileFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"FileLast") == 0){getValue(vf,ch2,ch3,value); this->fileLast = atoi(value); izt = true;}
		if(strcmp(childTag,"FileStep") == 0){getValue(vf,ch2,ch3,value); this->fileStep = atoi(value); izt = true;}
		if(strcmp(childTag,"ImageFirst") == 0){getValue(vf,ch2,ch3,value); this->imageFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"ImageLast") == 0){getValue(vf,ch2,ch3,value); this->imageLast = atoi(value); izt = true;}
		if(strcmp(childTag,"ImageStep") == 0){getValue(vf,ch2,ch3,value); this->imageStep = atoi(value); izt = true;}
		if(strcmp(childTag,"FirstImageIndex") == 0){getValue(vf,ch2,ch3,value); this->firstImageIndex = atoi(value); izt = true;}
		if(strcmp(childTag,"ImagesPerFile") == 0){getValue(vf,ch2,ch3,value); this->imagesPerFile = atoi(value); izt = true;}

		if(strcmp(childTag,"MemorySizeMin") == 0){getValue(vf,ch2,ch3,value); this->memorySizeMin = atoi(value); izt = true;}
		if(strcmp(childTag,"MemorySizeMax") == 0){getValue(vf,ch2,ch3,value); this->memorySizeFMax = float(atof(value)); izt = true;}

		if(strcmp(childTag,"ValueMin") == 0){getValue(vf,ch2,ch3,value); this->valueMin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ValueMax") == 0){getValue(vf,ch2,ch3,value); this->valueMax = float(atof(value)); izt = true;}
		if(strcmp(childTag,"PixelParam") == 0){getValue(vf,ch2,ch3,value); this->pixelParam = float(atof(value)); izt = true;}
				
		if(strcmp(childTag,"Raw") == 0){
			if(!this->raw->readXml(vf,ch1,ch4))
			  return false;
			izt = true;
		}

		if(strcmp(childTag,"Orientation") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"projection") == 0){
				orientation = INPUT_ORIENTATION_PROJECTION;
			}else if(strcmp(value,"slice") == 0){
				orientation = INPUT_ORIENTATION_SLICE;
			}else{
				printf("Unknown value \"%s\" for the tag \"Orientation\" (should be either \"Projection\" or \"Slice\").\n",value);
				return false;
			}
			izt = true;

		}

		if(strcmp(childTag,"Shape") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"point") == 0){
				shape = SHAPE_POINT;
			}else if(strcmp(value,"pixel") == 0){
				shape = SHAPE_PIXEL;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Point\" or \"Pixel\").\n",value,childTag);
				return false;
			}
			izt = true;

		}

		if(strcmp(childTag,"Restrictions") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"no") == 0){
				restrictions = INPUT_RESTRICTIONS_NO;
			}else if(strcmp(value,"yes") == 0){
				restrictions = INPUT_RESTRICTIONS_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"Restrictions\" (should be either \"Yes\" or \"No\").\n",value);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"intensity") == 0){
				type = INPUT_TYPE_INTENSITY;
			}else if(strcmp(value,"attenuation") == 0){
				type = INPUT_TYPE_ATTENUATION;
			}else{
				printf("Unknown value \"%s\" for the tag \"Type\" (should be either \"Intensity\" or \"Attenuation\").\n",value);
				return false;
			}
			izt = true;
		}
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}

	return true;
}




xOutputData::xOutputData(){
	folder[0] = '\0';
	prefix[0] = '\0';
	suffix[0] = '\0';
	extension[0] = '\0';
	NOD = 0;
	fileFirst = 0;
	fileStep = 1;
	bits = 16;
	bitsType = OUTPUT_BITS_GIVEN;
	valueMin = 0.0;
	valueMax = 0.0;
	restrictions = OUTPUT_RESTRICTIONS_NO;
}

bool xOutputData::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"OutputData");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;			
		
		if(strcmp(childTag,"Folder") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->folder,value); izt = true;}
		if(strcmp(childTag,"Prefix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->prefix,value); izt = true;}
		if(strcmp(childTag,"Suffix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->suffix,value); izt = true;}
		if(strcmp(childTag,"Extension") == 0){getValue(vf,ch2,ch3,value); sprintf(this->extension,value); izt = true;}
		if(strcmp(childTag,"NOD") == 0){getValue(vf,ch2,ch3,value); this->NOD = atoi(value); izt = true;}
		if(strcmp(childTag,"FileFirst") == 0){getValue(vf,ch2,ch3,value); this->fileFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"FileStep") == 0){getValue(vf,ch2,ch3,value); this->fileStep = atoi(value); izt = true;}

		if(strcmp(childTag,"Bits") == 0){getValue(vf,ch2,ch3,value); this->bits = atoi(value); izt = true;}
		if(strcmp(childTag,"ValueMin") == 0){getValue(vf,ch2,ch3,value); this->valueMin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ValueMax") == 0){getValue(vf,ch2,ch3,value); this->valueMax = float(atof(value)); izt = true;}

		if(strcmp(childTag,"BitsType") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"input") == 0){
				bitsType = OUTPUT_BITS_INPUT;
			}else if(strcmp(value,"given") == 0){
				bitsType = OUTPUT_BITS_GIVEN;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Input\" or \"Given\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"Restrictions") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"no") == 0){
				restrictions = OUTPUT_RESTRICTIONS_NO;
			}else if(strcmp(value,"yes") == 0){
				restrictions = OUTPUT_RESTRICTIONS_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"Restrictions\" (should be either \"Yes\" or \"No\").\n",value);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);

			if(strcmp(value,"solution") == 0){
				type = OUTPUT_TYPE_SOLUTION;
			}else if(strcmp(value,"ffcorrection") == 0){
				type = OUTPUT_TYPE_FFCORRECTION;
			}else if(strcmp(value,"preprocessed") == 0){
				type = OUTPUT_TYPE_PREPROCESSED;
			}else if(strcmp(value,"transformed") == 0){
				type = OUTPUT_TYPE_TRANSFORMED;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Solution\", \"FFCorrection\", \"Preprocessed\" or \"Transformed\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"State") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);

			if(strcmp(value,"intensity") == 0){
				state = OUTPUT_STATE_INTENSITY;
			}else if(strcmp(value,"attenuation") == 0){
				state = OUTPUT_STATE_ATTENUATION;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Intensity\" or \"Attenuation\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"Shape") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"point") == 0){
				shape = SHAPE_POINT;
			}else if(strcmp(value,"pixel") == 0){
				shape = SHAPE_PIXEL;
			}else if(strcmp(value,"input") == 0){
				shape = SHAPE_INPUT;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Point\", \"Pixel\" or \"Input\").\n",value,childTag);
				return false;
			}
			izt = true;

		}

		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	return true;
}




xPreprocessing::xPreprocessing(){
	highPeaksBefore = new xHighPeaks();
	highPeaksCols = new xHighPeaks();
	highPeaksRows = new xHighPeaks();
	ringArtefacts = new xRingArtefacts();
	intensity = new xIntensity();
}

xPreprocessing::~xPreprocessing(){
	delete highPeaksBefore; highPeaksBefore = NULL;
	delete highPeaksCols; highPeaksCols = NULL;
	delete highPeaksRows; highPeaksRows = NULL;
	delete ringArtefacts; ringArtefacts = NULL;
	delete intensity; intensity = NULL;
}



bool xPreprocessing::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100];
	sprintf(mTag,"Preprocessing");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}

	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
		
		if(strcmp(childTag,"HighPeaksBefore") == 0){
			if(!this->highPeaksBefore->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"RingArtefacts") == 0){
			if(!this->ringArtefacts->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"Intensity") == 0){
			if(!this->intensity->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"HighPeaksAfterColumns") == 0){
			if(!this->highPeaksCols->readXml(vf,ch1,ch4))return false;
			izt = true;
		}

		if(strcmp(childTag,"HighPeaksAfterRows") == 0){
			if(!this->highPeaksRows->readXml(vf,ch1,ch4))return false;
			izt = true;
		}
		
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}

	return true;
}





bool xRaw::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"Raw");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}

	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
			
		
		if(strcmp(childTag,"Bits") == 0){getValue(vf,ch2,ch3,value); this->bits = atoi(value); izt = true;}
		if(strcmp(childTag,"Offset") == 0){getValue(vf,ch2,ch3,value); this->offset = atoi(value); izt = true;}
		if(strcmp(childTag,"Gap") == 0){getValue(vf,ch2,ch3,value); this->gap = atoi(value); izt = true;}
		if(strcmp(childTag,"Xlen") == 0){getValue(vf,ch2,ch3,value); this->xlen = atoi(value); izt = true;}
		if(strcmp(childTag,"Ylen") == 0){getValue(vf,ch2,ch3,value); this->ylen = atoi(value); izt = true;}
		if(strcmp(childTag,"Zlen") == 0){getValue(vf,ch2,ch3,value); this->zlen = atoi(value); izt = true;}
		
		if(strcmp(childTag,"ByteOrder") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"be") == 0){
				byteOrder = BYTE_ORDER_BIG;
			}else if(strcmp(value,"le") == 0){
				byteOrder = BYTE_ORDER_LITTLE;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"BE\" or \"LE\").\n",value,childTag);
				return false;
			}
			izt = true;
		}
		
		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"yes") == 0){
				type = RAW_TYPE_YES;
			}else if(strcmp(value,"no") == 0){
				type = RAW_TYPE_NO;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Yes\" or \"No\").\n",value,childTag);
				return false;
			}
			izt = true;
		}
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}

	return true;
}


bool xXrm2tif::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"Xrm2tif");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;

		if(strcmp(childTag,"InputFile") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->inputFile,value); izt = true;}			
		if(strcmp(childTag,"Folder") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->folder,value); izt = true;}
		if(strcmp(childTag,"LogFile") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->logFile,value); izt = true;}
		if(strcmp(childTag,"TagFile") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->tagFile,value); izt = true;}

		if(strcmp(childTag,"Prefix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->prefix,value); izt = true;}
		if(strcmp(childTag,"Suffix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->suffix,value); izt = true;}
		if(strcmp(childTag,"Extension") == 0){getValue(vf,ch2,ch3,value); sprintf(this->extension,value); izt = true;}

		if(strcmp(childTag,"FlatName") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->flatName,value); izt = true;}
		if(strcmp(childTag,"DarkName") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->darkName,value); izt = true;}

		if(strcmp(childTag,"NOD") == 0){getValue(vf,ch2,ch3,value); this->NOD = atoi(value); izt = true;}
		if(strcmp(childTag,"FileFirst") == 0){getValue(vf,ch2,ch3,value); this->fileFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"FileStep") == 0){getValue(vf,ch2,ch3,value); this->fileStep = atoi(value); izt = true;}

		if(strcmp(childTag,"ImageFirst") == 0){getValue(vf,ch2,ch3,value); this->imageFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"ImageLast") == 0){getValue(vf,ch2,ch3,value); this->imageLast = atoi(value); izt = true;}
		if(strcmp(childTag,"ImageStep") == 0){getValue(vf,ch2,ch3,value); this->imageStep = atoi(value); izt = true;}

		if(strcmp(childTag,"RowFirst") == 0){getValue(vf,ch2,ch3,value); this->rowFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"RowLast") == 0){getValue(vf,ch2,ch3,value); this->rowLast = atoi(value); izt = true;}
		
		if(strcmp(childTag,"ImagesPerFile") == 0){getValue(vf,ch2,ch3,value); this->imagesPerFile = atoi(value); izt = true;}

		if(strcmp(childTag,"MemorySizeMin") == 0){getValue(vf,ch2,ch3,value); this->memorySizeMin = atoi(value); izt = true;}
		if(strcmp(childTag,"MemorySizeMax") == 0){getValue(vf,ch2,ch3,value); this->memorySizeFMax = float(atof(value)); izt = true;}

		if(strcmp(childTag,"OutputType") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"all") == 0){
				this->output_type = XRM_OUTPUT_ALL;
			}else if(strcmp(value,"xml") == 0){
				this->output_type = XRM_OUTPUT_XML;
			}else{
				printf("Unknown value \"%s\" for the tag \"Type\" (should be either \"All\" or \"Xml\").\n",value);
				return false;
			}
			izt = true;
		}
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}

	return true;
}


bool xBeamlineUser::readXml(Ipp8u *vf, int s1, int s2){
	return true;
}




bool xFDField::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"FDField");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
		
		if(strcmp(childTag,"FileProfile") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->fileProfile,value); izt = true;}
		if(strcmp(childTag,"FileBefore") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->fileBefore,value); izt = true;}
		if(strcmp(childTag,"FileAfter") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->fileAfter,value); izt = true;}

		
		if(strcmp(childTag,"ValueBefore") == 0){getValue(vf,ch2,ch3,value); this->valueBefore = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ValueAfter") == 0){getValue(vf,ch2,ch3,value); this->valueAfter = float(atof(value)); izt = true;}

		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"user") == 0){
				type = FIELD_TYPE_USER;

			}else if(strcmp(value,"row") == 0){
				type = FIELD_TYPE_ROW;
			}else{
				printf("Unknown value \"%s\" for the tag \"Type\" (should be either \"User\" or \"Row\").\n",value);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"ProfileType") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"constant") == 0){
				typeProfile = FIELD_PROFILE_CONSTANT;
			}else if(strcmp(value,"linear") == 0){
				typeProfile = FIELD_PROFILE_LINEAR;
			}else if(strcmp(value,"profile") == 0){
				typeProfile = FIELD_PROFILE_PROFILE;
			}
			else{
				printf("Unknown value \"%s\" for the tag \"Type\" (should be either \"Constant\", \"Linear\" or \"Profile\").\n",value);
				return false;
			}
			izt = true;
		}

		if(!izt){

			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	
	return true;

}



bool xRingArtefacts::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"RingArtefacts");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
				
		if(strcmp(childTag,"ParameterN") == 0){getValue(vf,ch2,ch3,value); this->parameterN = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ParameterR") == 0){getValue(vf,ch2,ch3,value); this->parameterR = float(atof(value)); izt = true;}
		if(strcmp(childTag,"NumSeries") == 0){getValue(vf,ch2,ch3,value); this->num_series = atoi(value); izt = true;}

		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				type = RING_ARTEFACTS_NO;
			}else if(strcmp(value,"column") == 0){
				type = RING_ARTEFACTS_COLUMN;
			}else if(strcmp(value,"aml") == 0){
				type = RING_ARTEFACTS_AML;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\", \"Column\" or \"Aml\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	
	return true;

}


bool xTransform::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"Transform");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;

		if(strcmp(childTag,"RotationAngle") == 0){getValue(vf,ch2,ch3,value); this->rotationAngle = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ReCentreAngle") == 0){getValue(vf,ch2,ch3,value); this->reCentreAngle = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ReCentreRadius") == 0){getValue(vf,ch2,ch3,value); this->reCentreRadius = float(atof(value)); izt = true;}

		if(strcmp(childTag,"CropTop") == 0){getValue(vf,ch2,ch3,value); this->cropTop = atoi(value); izt = true;}
		if(strcmp(childTag,"CropBottom") == 0){getValue(vf,ch2,ch3,value); this->cropBottom = atoi(value); izt = true;}
		if(strcmp(childTag,"CropLeft") == 0){getValue(vf,ch2,ch3,value); this->cropLeft = atoi(value); izt = true;}
		if(strcmp(childTag,"CropRight") == 0){getValue(vf,ch2,ch3,value); this->cropRight = atoi(value); izt = true;}

		if(strcmp(childTag,"ScaleWidth") == 0){getValue(vf,ch2,ch3,value); this->scaleWidth = atoi(value); izt = true;}
		if(strcmp(childTag,"ScaleHeight") == 0){getValue(vf,ch2,ch3,value); this->scaleHeight = atoi(value); izt = true;}
		if(strcmp(childTag,"ExtrapolationWidth") == 0){getValue(vf,ch2,ch3,value); this->extrapolationWidth = atoi(value); izt = true;}
		if(strcmp(childTag,"ExtrapolationPixels") == 0){getValue(vf,ch2,ch3,value); this->extrapolationPixels = atoi(value); izt = true;}

		if(strcmp(childTag,"MissedProjections") == 0){getValue(vf,ch2,ch3,value); sprintf(this->missedProjections,value); izt = true;}

		if(strcmp(childTag,"RotationAngleType") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"180") == 0){
				this->rotationAngleType = ROTATION_ANGLE_180;
			}else if(strcmp(value,"other") == 0){
				this->rotationAngleType = ROTATION_ANGLE_OTHER;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"180\" or \"Other\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"Interpolation") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"nearestneighbour") == 0){
				this->interpolation = IPPI_INTER_NN;
			}else if(strcmp(value,"linear") == 0){
				this->interpolation = IPPI_INTER_LINEAR;
			}else if(strcmp(value,"cubic") == 0){
				this->interpolation = IPPI_INTER_CUBIC;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"NearestNeighbour\", \"Linear\" or \"Cubic\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"MissedProjectionsType") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				missedProjectionsType = MISSED_PROJECTIONS_NO;
			}else if(strcmp(value,"linear") == 0){
				missedProjectionsType = MISSED_PROJECTIONS_LINEAR;
			}else if(strcmp(value,"zero") == 0){
				missedProjectionsType = MISSED_PROJECTIONS_ZERO;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\", \"Linear\" or \"Zero\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"RotationAngleEndPoints") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				this->rotationAngleEndPoints = ROTATION_ANGLE_END_NO;
			}else if(strcmp(value,"yes") == 0){
				this->rotationAngleEndPoints = ROTATION_ANGLE_END_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\" or \"Yes\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"ScaleType") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				this->scaleType = SCALE_NO;
			}else if(strcmp(value,"yes") == 0){
				this->scaleType = SCALE_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\" or \"Yes\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"ExtrapolationType") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				this->extrapolationType = EXTRAPOLATION_NO;
			}else if(strcmp(value,"linear") == 0){
				this->extrapolationType = EXTRAPOLATION_ZERO;
			}else if(strcmp(value,"constant") == 0){
				this->extrapolationType = EXTRAPOLATION_CONSTANT;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\", \"Zero\" or \"Constant\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	
	return true;

}


bool xHighPeaks::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"HighPeaks");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
				
		if(strcmp(childTag,"Jump") == 0){getValue(vf,ch2,ch3,value); this->jump = float(atof(value)); izt = true;}

		if(strcmp(childTag,"NumberPixels") == 0){getValue(vf,ch2,ch3,value); this->numberPixels = atoi(value); izt = true;}

		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				type = HIGH_PEAKS_NO;
			}else if(strcmp(value,"yes") == 0){
				type = HIGH_PEAKS_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\" or \"Yes\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	
	return true;

}


bool xRadius::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"Radius");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
				
		if(strcmp(childTag,"Percent") == 0){getValue(vf,ch2,ch3,value); this->percent = float(atof(value)); izt = true;}

		if(strcmp(childTag,"Pixel") == 0){getValue(vf,ch2,ch3,value); this->pixel = float(atof(value)); izt = true;}


		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"pixel") == 0){
				type = RADIUS_PIXEL;
			}else if(strcmp(value,"percent") == 0){
				type = RADIUS_PERCENT;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Pixel\" or \"Percent\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	
	return true;

}

bool xFilter::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"Filter");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
			
		

		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"yes") == 0){
				type = FILTER_YES;
			}else if(strcmp(value,"no") == 0){
				type = FILTER_NO;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Yes\" or \"No\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"Bandwidth") == 0){getValue(vf,ch2,ch3,value); this->bandwidth = float(atof(value)); izt = true;}
		
		if(strcmp(childTag,"PixelSize") == 0){getValue(vf,ch2,ch3,value); this->pixelSize = float(atof(value)); izt = true;}

		if(strcmp(childTag,"Normalisation") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"cc") == 0){
				this->norm = FILTER_NORM_CC;
			}else if(strcmp(value,"ca") == 0){
				this->norm = FILTER_NORM_CA;
			}else if(strcmp(value,"pc") == 0){
				this->norm = FILTER_NORM_PC;
			}else if(strcmp(value,"pa") == 0){
				this->norm = FILTER_NORM_PA;
			}else if(strcmp(value,"sc") == 0){
				this->norm = FILTER_NORM_SC;
			}else if(strcmp(value,"sa") == 0){
				this->norm = FILTER_NORM_SA;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"CC\", \"CA\", \"PC\", \"PA\", \"SC\" or \"SA\").\n",value,childTag);
				return false;
			}
			izt = true;
		}


		if(strcmp(childTag,"Name") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"rl") == 0){
				name = FILTER_NAME_RL;
			}else if(strcmp(value,"sl") == 0){
				name = FILTER_NAME_SL;
			}else if(strcmp(value,"cos") == 0){
				name = FILTER_NAME_COS;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"RL\", \"SL\" or \"Cos\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"WindowName") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				windowName = FILTER_WINDOW_NO;
			}else if(strcmp(value,"hann") == 0){
				windowName = FILTER_WINDOW_HANN;
			}else if(strcmp(value,"hamming") == 0){
				windowName = FILTER_WINDOW_HAMMING;
			}else if(strcmp(value,"blackman") == 0){
				windowName = FILTER_WINDOW_BLACKMAN;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\", \"Hann\", \"Namming\" or \"Blackman\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	
	return true;

}


bool xIntensity::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"Intensity");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;
				
		

		if(strcmp(childTag,"ColumnLeft") == 0){getValue(vf,ch2,ch3,value); this->columnLeft = atoi(value); izt = true;}
		if(strcmp(childTag,"ColumnRight") == 0){getValue(vf,ch2,ch3,value); this->columnRight = atoi(value); izt = true;}
		if(strcmp(childTag,"ZeroLeft") == 0){getValue(vf,ch2,ch3,value); this->zeroLeft = atoi(value); izt = true;}
		if(strcmp(childTag,"ZeroRight") == 0){getValue(vf,ch2,ch3,value); this->zeroRight = atoi(value); izt = true;}

		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"no") == 0){
				type = INTENSITY_NO;
			}else if(strcmp(value,"row") == 0){
				type = INTENSITY_ROW;
			}else if(strcmp(value,"column") == 0){
				type = INTENSITY_COLUMN;
			}else if(strcmp(value,"zero") == 0){
				type = INTENSITY_ZERO;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\", \"Column\", \"Row\" or \"Zero\").\n",value,childTag);
				return false;
			}
			izt = true;
		}


		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	
	return true;
}




bool xRoi::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"ROI");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;			
		

		if(strcmp(childTag,"Xmin") == 0){getValue(vf,ch2,ch3,value); this->xmin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"Xmax") == 0){getValue(vf,ch2,ch3,value); this->xmax = float(atof(value)); izt = true;}
		if(strcmp(childTag,"Ymin") == 0){getValue(vf,ch2,ch3,value); this->ymin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"Ymax") == 0){getValue(vf,ch2,ch3,value); this->ymax = float(atof(value)); izt = true;}
		if(strcmp(childTag,"Angle") == 0){getValue(vf,ch2,ch3,value); this->angle = float(atof(value)); izt = true;}

		if(strcmp(childTag,"OutputWidth") == 0){getValue(vf,ch2,ch3,value); this->outputWidth = atoi(value); izt = true;}

		if(strcmp(childTag,"Type") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"standard") == 0){
				type = ROI_STANDARD;
			}else if(strcmp(value,"rectangle") == 0){
				type = ROI_RECTANGLE;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Standard\" or \"Rectangle\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		

		if(strcmp(childTag,"OutputWidthType") == 0){
			getValue(vf,ch2,ch3,value);			
			toLower(value);
			if(strcmp(value,"standard") == 0){
				outputWidthType = OUTPUT_WIDTH_STANDARD;
			}else if(strcmp(value,"given") == 0){
				outputWidthType = OUTPUT_WIDTH_GIVEN;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Standard\" or \"Given\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}
	return true;
}



bool xAllFDK::readXml(Ipp8u *vf, int s1, int s2){
	char tagName[100], childTag[100], mTag[100], value[MAX_FOLDER];
	sprintf(mTag,"AllFDK");
	int t1, t2, t3, t4;
	int ch1, ch2, ch3, ch4;
	bool t, izt, isEmptyTag;

	t = findNextTag(vf, s1, s2, &t1, &t2, tagName, &isEmptyTag);
	if(!t){
		printf("Can not find the start-tag \"%s\".\n",mTag);
		return false;
	}

	t = findEndTag(vf, t1, s2, &t3, &t4, tagName);
	if(!t){
		printf("Can not find the end-tag \"%s\".\n",mTag);
		return false;
	}

	if(strcmp(tagName,mTag) != 0){
		printf("Wrong tag: \"%s\" (should be \"%s\").\n",tagName,mTag);
		return false;
	}
	
	for(;;){
		t = findNextTag(vf, t2, t3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;

		// Skip an empty <tag/>
		if( isEmptyTag ) {
		  t2 = ch2+1;
		  continue;
		}

		t = findEndTag(vf, ch2, t3, &ch3, &ch4, childTag);
		if(!t){
			printf("Can not find the child end-tag \"%s\".\n",childTag);
			return false;
		}
		izt  = false;

		if(strcmp(childTag,"InputFolder") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->inputFolder,value); izt = true;}
		if(strcmp(childTag,"OutputFolder") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->outputFolder,value); izt = true;}

		if(strcmp(childTag,"LogFile") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->logFile,value); izt = true;}
		
		if(strcmp(childTag,"InputPrefix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->inputPrefix,value); izt = true;}
		if(strcmp(childTag,"InputSuffix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->inputSuffix,value); izt = true;}
		if(strcmp(childTag,"InputExtension") == 0){getValue(vf,ch2,ch3,value); sprintf(this->inputExtension,value); izt = true;}

		if(strcmp(childTag,"OutputPrefix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->outputPrefix,value); izt = true;}
		if(strcmp(childTag,"OutputSuffix") == 0){getValue(vf,ch2,ch3,value); sprintf(this->outputSuffix,value); izt = true;}
		if(strcmp(childTag,"OutputExtension") == 0){getValue(vf,ch2,ch3,value); sprintf(this->outputExtension,value); izt = true;}

		if(strcmp(childTag,"FlatName") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->flatName,value); izt = true;}
		if(strcmp(childTag,"DarkName") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->darkName,value); izt = true;}
		if(strcmp(childTag,"FilterName") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->filterName,value); izt = true;}

		if(strcmp(childTag,"XShiftsFile") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->xShiftsFile,value); izt = true;}
		if(strcmp(childTag,"YShiftsFile") == 0){getValue(vf,ch2,ch3,value); replaceSlash(value); sprintf(this->yShiftsFile,value); izt = true;}

		if(strcmp(childTag,"NumberOfCircles") == 0){getValue(vf,ch2,ch3,value); this->numberOfCircles = atoi(value); izt = true;}
		if(strcmp(childTag,"InputNOD") == 0){getValue(vf,ch2,ch3,value); this->inputNOD = atoi(value); izt = true;}
		if(strcmp(childTag,"OutputNOD") == 0){getValue(vf,ch2,ch3,value); this->outputNOD = atoi(value); izt = true;}
	
		if(strcmp(childTag,"InputFileFirst") == 0){getValue(vf,ch2,ch3,value); this->inputFileFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"InputFileStep") == 0){getValue(vf,ch2,ch3,value); this->inputFileStep = atoi(value); izt = true;}

		if(strcmp(childTag,"InputImageFirst") == 0){getValue(vf,ch2,ch3,value); this->inputImageFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"InputImageLast") == 0){getValue(vf,ch2,ch3,value); this->inputImageLast = atoi(value); izt = true;}
		if(strcmp(childTag,"InputImageStep") == 0){getValue(vf,ch2,ch3,value); this->inputImageStep = atoi(value); izt = true;}
		if(strcmp(childTag,"InputImagesPerFile") == 0){getValue(vf,ch2,ch3,value); this->inputImagesPerFile = atoi(value); izt = true;}

		if(strcmp(childTag,"RowFirst") == 0){getValue(vf,ch2,ch3,value); this->rowFirst = atoi(value); izt = true;}

		if(strcmp(childTag,"OutputIndexFirst") == 0){getValue(vf,ch2,ch3,value); this->outputIndexFirst = atoi(value); izt = true;}
		if(strcmp(childTag,"OutputIndexStep") == 0){getValue(vf,ch2,ch3,value); this->outputIndexStep = atoi(value); izt = true;}
		if(strcmp(childTag,"OutputWidth") == 0){getValue(vf,ch2,ch3,value); this->outputWidth = atoi(value); izt = true;}

		if(strcmp(childTag,"FlatValue") == 0){getValue(vf,ch2,ch3,value); this->flatValue = float(atof(value)); izt = true;}
		if(strcmp(childTag,"DarkValue") == 0){getValue(vf,ch2,ch3,value); this->darkValue = float(atof(value)); izt = true;}

		if(strcmp(childTag,"FlatMin") == 0){getValue(vf,ch2,ch3,value); this->flatMin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"FlatMax") == 0){getValue(vf,ch2,ch3,value); this->flatMax = float(atof(value)); izt = true;}
		if(strcmp(childTag,"DarkMin") == 0){getValue(vf,ch2,ch3,value); this->darkMin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"DarkMax") == 0){getValue(vf,ch2,ch3,value); this->darkMax = float(atof(value)); izt = true;}
		if(strcmp(childTag,"InputMin") == 0){getValue(vf,ch2,ch3,value); this->inputMin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"InputMax") == 0){getValue(vf,ch2,ch3,value); this->inputMax = float(atof(value)); izt = true;}
		if(strcmp(childTag,"OutputMin") == 0){getValue(vf,ch2,ch3,value); this->outputMin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"OutputMax") == 0){getValue(vf,ch2,ch3,value); this->outputMax = float(atof(value)); izt = true;}

		if(strcmp(childTag,"SliceFirst") == 0){getValue(vf,ch2,ch3,value); this->sliceFirst = float(atof(value)); izt = true;}
		if(strcmp(childTag,"SliceLast") == 0){getValue(vf,ch2,ch3,value); this->sliceLast = float(atof(value)); izt = true;}
		if(strcmp(childTag,"SliceStep") == 0){getValue(vf,ch2,ch3,value); this->sliceStep = float(atof(value)); izt = true;}

		if(strcmp(childTag,"ROIXmin") == 0){getValue(vf,ch2,ch3,value); this->ROIXmin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ROIXmax") == 0){getValue(vf,ch2,ch3,value); this->ROIXmax = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ROIYmin") == 0){getValue(vf,ch2,ch3,value); this->ROIYmin = float(atof(value)); izt = true;}
		if(strcmp(childTag,"ROIYmax") == 0){getValue(vf,ch2,ch3,value); this->ROIYmax = float(atof(value)); izt = true;}

		if(strcmp(childTag,"SourceXShift") == 0){getValue(vf,ch2,ch3,value); this->sourceXShift = float(atof(value)); izt = true;}
		if(strcmp(childTag,"SourceYShift") == 0){getValue(vf,ch2,ch3,value); this->sourceYShift = float(atof(value)); izt = true;}
		
		if(strcmp(childTag,"AxisXShift") == 0){getValue(vf,ch2,ch3,value); this->axisXShift = float(atof(value)); izt = true;}
		if(strcmp(childTag,"PixelSize") == 0){getValue(vf,ch2,ch3,value); this->pixelSize = float(atof(value)); izt = true;}
		if(strcmp(childTag,"AngleStep") == 0){getValue(vf,ch2,ch3,value); this->angleStep = float(atof(value)); izt = true;}

		if(strcmp(childTag,"RotationAngle") == 0){getValue(vf,ch2,ch3,value); this->rotationAngle = float(atof(value)); izt = true;}
		if(strcmp(childTag,"SourceToObject") == 0){getValue(vf,ch2,ch3,value); this->sourceToObject = float(atof(value)); izt = true;}
		if(strcmp(childTag,"DetectorToObject") == 0){getValue(vf,ch2,ch3,value); this->detectorToObject = float(atof(value)); izt = true;}
		if(strcmp(childTag,"RingArtefactsParamR") == 0){getValue(vf,ch2,ch3,value); this->ringArtefactsParamR = float(atof(value)); izt = true;}
		if(strcmp(childTag,"RingArtefactsParamN") == 0){getValue(vf,ch2,ch3,value); this->ringArtefactsParamN = float(atof(value)); izt = true;}

		if(strcmp(childTag,"OriginalImageHeight") == 0){getValue(vf,ch2,ch3,value); this->originalImageHeight = float(atof(value)); izt = true;}
		if(strcmp(childTag,"FilterParam") == 0){getValue(vf,ch2,ch3,value); this->filterParam = float(atof(value)); izt = true;}
	
		
		if(strcmp(childTag,"OutputBits") == 0){
			getValue(vf,ch2,ch3,value);
			this->outputBits = atoi(value);
			if(this->outputBits != 8 && outputBits != 16 && outputBits != 32){
				printf("Wrong value \"%s\" for the tag \"OutputBits\" (should be either \"8\", \"16\" or \"32\").\n",value);
				return false;
			}
			izt = true;
		}
		

		if(strcmp(childTag,"XShiftsAbs") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"absolute") == 0){
				this->xShiftsAbs = SHIFTS_ABSOLUTE;
			}else if(strcmp(value,"relative") == 0){
				this->xShiftsType = SHIFTS_RELATIVE;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Absolute\"  or \"Relative\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"YShiftsAbs") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"absolute") == 0){
				this->yShiftsAbs = SHIFTS_ABSOLUTE;
			}else if(strcmp(value,"relative") == 0){
				this->yShiftsType = SHIFTS_RELATIVE;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Absolute\"  or \"Relative\").\n",value,childTag);
				return false;
			}
			izt = true;
		}


		if(strcmp(childTag,"XShiftsType") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"no") == 0){
				this->xShiftsType = SHIFTS_NO;
			}else if(strcmp(value,"positive") == 0){
				this->xShiftsType = SHIFTS_POSITIVE;
			}else if(strcmp(value,"negative") == 0){
				this->xShiftsType = SHIFTS_NEGATIVE;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\", \"Positive\" or \"Negative\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"YShiftsType") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"no") == 0){
				this->yShiftsType = SHIFTS_NO;
			}else if(strcmp(value,"positive") == 0){
				this->yShiftsType = SHIFTS_POSITIVE;
			}else if(strcmp(value,"negative") == 0){
				this->yShiftsType = SHIFTS_NEGATIVE;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\", \"Positive\" or \"Negative\").\n",value,childTag);
				return false;
			}
			izt = true;
		}


		if(strcmp(childTag,"OutputWidthType") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"standard") == 0){
				this->outputWidthType = CLOCKWISE_ROTATION_NO;
			}else if(strcmp(value,"given") == 0){
				this->outputWidthType = CLOCKWISE_ROTATION_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Standard\" or \"Given\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"ClockwiseRotation") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"no") == 0){
				this->clockwiseRotation = CLOCKWISE_ROTATION_NO;
			}else if(strcmp(value,"yes") == 0){
				this->clockwiseRotation = CLOCKWISE_ROTATION_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\" or \"Yes\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"RingArtefactsType") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"no") == 0){
				this->ringArtefactsType = RING_ARTEFACTS_NO;
			}else if(strcmp(value,"column") == 0){
				this->ringArtefactsType = RING_ARTEFACTS_COLUMN;
			}else if(strcmp(value,"aml") == 0){
				this->ringArtefactsType = RING_ARTEFACTS_AML;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\", \"Column\" or \"AML\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"FlatType") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"value") == 0){
				this->flatType = FIELD_TYPE_VALUE;
			}else if(strcmp(value,"file") == 0){
				this->flatType = FIELD_TYPE_FILE;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Value\" or \"File\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"DarkType") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"value") == 0){
				this->darkType = FIELD_TYPE_VALUE;
			}else if(strcmp(value,"file") == 0){
				this->darkType = FIELD_TYPE_FILE;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"Value\" or \"File\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		if(strcmp(childTag,"InputRestrictions") == 0){
			getValue(vf,ch2,ch3,value);
			toLower(value);
			if(strcmp(value,"no") == 0){
				this->inputRestrictions = INPUT_RESTRICTIONS_NO;
			}else if(strcmp(value,"yes") == 0){
				this->inputRestrictions = INPUT_RESTRICTIONS_YES;
			}else{
				printf("Unknown value \"%s\" for the tag \"%s\" (should be either \"No\" or \"Yes\").\n",value,childTag);
				return false;
			}
			izt = true;
		}

		
		if(!izt){
			printf("Unknown child tag \"%s\" of the tag \"%s\".\n",childTag,mTag);
			return false;
		}
		t2 = ch4+1;
	}

	return true;
}


