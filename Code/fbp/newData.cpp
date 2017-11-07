/* ----------------------------------------------------------------------
 * newData.cpp
 * ----------------------------------------------------------------------
 * Implement the methods in the xData class.
 *
 * ----------------------------------------------------------------------
 */

#include "defs.h"


xData::xData(){
	matd1 = NULL;
	matd2 = NULL;
	matf1 = NULL;
	matf2 = NULL;
	profd = NULL;
	proff = NULL;
	veci = NULL;
	vecio = NULL;
	vecd = NULL;
	vecf = NULL;
	vect1 = NULL;
	vect2 = NULL;
	veco = NULL;
	vecp1 = NULL;
	vecp2 = NULL;
	vec32 = NULL;
	vec64 = NULL;
	vec_res = NULL;
	mata = NULL;
	vmiss = NULL;
	vtb = NULL;
	vta = NULL;
	mapx = NULL;
	mapy = NULL;
	shift_x = NULL;
	shift_y = NULL;
	vecSinCos = NULL;
	vta_warp = NULL;
}

xData::~xData(){
	freeMemory();
}

void xData::freeMemory(){
	if(matd1 != NULL){ippsFree(matd1); matd1 = NULL;}
	if(matd2 != NULL){ippsFree(matd2); matd2 = NULL;}
	if(matf1 != NULL){ippsFree(matf1); matf1 = NULL;}
	if(matf2 != NULL){ippsFree(matf2); matf2 = NULL;}
	if(profd != NULL){ippsFree(profd); profd = NULL;}
	if(proff != NULL){ippsFree(proff); proff = NULL;}
	if(veci != NULL){ippsFree(veci); veci = NULL;}
	if(vecio != NULL){ippsFree(vecio); vecio = NULL;}
	if(vecd != NULL){ippsFree(vecd); vecd = NULL;}
	if(vecf != NULL){ippsFree(vecf); vecf = NULL;}
	if(vect1 != NULL){ippsFree(vect1); vect1 = NULL;}
	if(vect2 != NULL){ippsFree(vect2); vect2 = NULL;}
	veco = NULL;
	//if(veco != NULL){ippsFree(veco); veco = NULL;}
	if(vecp1 != NULL){ippsFree(vecp1); vecp1 = NULL;}
	if(vecp2 != NULL){ippsFree(vecp2); vecp2 = NULL;}
	if(vec32 != NULL){ippsFree(vec32); vec32 = NULL;}
	if(vec64 != NULL){ippsFree(vec64); vec64 = NULL;}
	if(vec_res != NULL){ippsFree(vec_res); vec_res = NULL;}
	if(mata != NULL){ippsFree(mata); mata = NULL;}
	if(vmiss != NULL){ippsFree(vmiss); vmiss = NULL;}
	if(vtb != NULL){ippsFree(vtb); vtb = NULL;}
	if(vta != NULL){ippsFree(vta); vta = NULL;}
	if(mapx != NULL){ippsFree(mapx); mapx = NULL;}
	if(mapy != NULL){ippsFree(mapy); mapy = NULL;}
	if(shift_x != NULL){ippsFree(shift_x); shift_x = NULL;}
	if(shift_y != NULL){ippsFree(shift_y); shift_y = NULL;}
	if(vecSinCos != NULL){ippsFree(vecSinCos); vecSinCos = NULL;}
	if(vta_warp != NULL){ippsFree(vta_warp); vta_warp = NULL;}
}
