#ifndef CALLBACK_DATA_H
#define CALLBACK_DATA_H

#include <set>
#include <vector>
#include <string>

#include "../include/Graph.h"
#include "../include/Variable.h"

class CBData
{
public:
	CBData();
	CBData(Graph* _g, TypeVariableHashPtr _vHashPtr, int _nCols);
	virtual ~CBData();

	void reset();

	static int totalNbCalls;
	static int totalNbIntCuts;	// Total amount of integer cuts separated so far.
	static int totalNbFracCuts;	// Total amount of fractional cuts separated so far.
	static int totalNbCuts;		// Total amount of cuts separated so far.
	static int maxSep;			// Maximum number of total separation loop per bound of BB tree

	int seprounds;				// Counter of total separation loop per bound of BB tree

	int numCols;				// Amount of columns in problem matrix
	int numCuts;				// Amount of generated cuts by separation routine

	int currnode;				// Stores the number of current explored node
	int last_node;				// Stores the number of the last immediate node

	double currobj;
	double lastobj;

	std::vector<double> nodeSolX;	// Node solution vector

	Graph* graph;				// Pointer to graph structure

	TypeVariableHashPtr vHashPtr;

	std::set<long>* hashTable;
};

#endif