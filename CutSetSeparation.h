#ifndef CUT_SET_SEPARATION_H
#define CUT_SET_SEPARATION_H

#include "Util.h"
#include "Separator.h"
#include "../include/Graph.h"
#include "../include/Variable.h"
#include "../include/CallbackData.h"

#include <cplex.h>

int
CPXPUBLIC separateIntegerCutSet(
	CPXCENVptr _env,
	void       *_cbdata,
	int        _wherefrom,
	void       *_cbhandle,
	int        *_useraction_p);

int
CPXPUBLIC separateFractionalCutSet(
	CPXCENVptr _env,
	void       *_cbdata,
	int        _wherefrom,
	void       *_cbhandle,
	int        *_useraction_p);

int __stdcall integerCutSet(Graph* _g, TypeVariableHashPtr _vHashPtr, double * _x, std::set<int> & _s);

int 
__stdcall fractionalCutSet(
	Graph*							_g,
	TypeVariableHashPtr				_vHashPtr,
	double*							_x,
	std::vector<std::set<int>* > &	_cutSets,
	bool							_findAllCutsSets);

int __stdcall fractionalCutSet(Graph* _g, TypeVariableHashPtr _vHashPtr, double * _x, std::set<int> & _s);

#endif