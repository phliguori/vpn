#ifndef FACILITY_LOCATION_H
#define FACILITY_LOCATION_H

#include <fstream>

#include "SolverMIP.h"
#include "Separator.h"
#include "../include/Graph.h"
#include "../include/CallbackData.h"
#include "../include/Variable.h"
#include "../include/Constraint.h"

typedef int(__stdcall *SingleCutSetSeparation)(Graph* _g, TypeVariableHashPtr _vHashPtr, double* _x, std::set<int> & _s);
typedef int(__stdcall *MulltipleCutSetSeparation)(Graph* _g, TypeVariableHashPtr _vHashPtr, double * _x, std::vector<std::set<int>* > & _cutSets, bool _findAllCutsSets);

class FacilityLocation : public SolverMIP
{
public:
	FacilityLocation(Graph* _g);

	virtual ~FacilityLocation();

	void setSeparationRoutine(void* _cbdata = NULL, Separator* _separator = NULL);

	SOLVERSTAT solveLRExtentedFormulations(MulltipleCutSetSeparation _sepFunc, double _tol = 1e-5);

	void builSolutionGraph();

	std::string printXSol();

	void solutionToDot(std::ostream & _output = std::cout) const;

	int createVariables();

	int createConstraints();

	int getActiveRouters() const { return activeRouters; }

	int getActiveArcs() const { return activeArcs; }

	// Variable hash
	TypeVariableHash vHash;

	// Constraint hash
	ConstraintHash cHash;

protected:
	int createVarX();

	int createVarY();

	int createVarZ();

	int createConsAssignment();

	int createConsFacilityActivation();

	int createConsFacilityActivation(int);

	int createConsTree();

	//TODO
	int createConsRadiusDistance();

	//TODO
	int createConsEdgeActivation();

	int activeRouters;

	int activeArcs;

	Graph* g;	// Base graph

	Graph* s;	// Solution graph
};

#endif