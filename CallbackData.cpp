#include "../include/CallbackData.h"

int CBData::totalNbCalls = 0;
int CBData::totalNbCuts = 0;
int CBData::totalNbIntCuts = 0;
int CBData::totalNbFracCuts = 0;
int CBData::maxSep = 0;

CBData::CBData()
{
	seprounds = 0;
	numCols = 0;
	numCuts = 0;
	currnode = 0;
	last_node = 0;
	currobj = 0.;
	lastobj = 0.;
	graph = NULL;
	vHashPtr = NULL;
	hashTable = NULL;

	return;
}

CBData::CBData(Graph * _g, TypeVariableHashPtr _vHashPtr, int _nCols)
	: CBData()
{
	graph = _g;
	vHashPtr = _vHashPtr;
	numCols = _nCols;
	nodeSolX.clear();
	nodeSolX.resize(_nCols, 0.);
	hashTable = new std::set<long>;
}

CBData::~CBData()
{
	if (hashTable != NULL)
	{
		delete hashTable;
		hashTable = NULL;
	}

	return;
}

void CBData::reset()
{
	seprounds = 0;
	numCols = 0;
	numCuts = 0;
	currnode = 0;
	last_node = 0;
	currobj = 0.;
	lastobj = 0.;
	totalNbCalls = 0;
	totalNbCuts = 0;
	totalNbIntCuts = 0;
	totalNbFracCuts = 0;
}