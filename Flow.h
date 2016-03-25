#ifndef FLOW_H
#define FLOW_H

#include "SolverMIP.h"
#include "Separator.h"

#include "../include/Graph.h"
#include "../include/Variable.h"
#include "../include/Constraint.h"

class Flow : public SolverMIP
{
public:
	Flow(Graph * _g);

	int createVariables();
	int createConstraints();

	std::string  printXSol();
	void builSolutionGraph();
	void solutionToDot(std::ostream & _output) const;

	// Variable hash
	TypeVariableHash vHash;

	// Constraint hash
	ConstraintHash cHash;

protected:
	int createVarX();

	int createVarY();

	int createVarW();

	int createConsFlowConservation();

	int createConsCapacityDetermination();

	int createConsDualFeasibility();

	Graph* g;	// Base graph
	Graph* s;	// solution graph
};
#endif
