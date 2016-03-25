#include "../include/Constraint.h"

#include <sstream>

int Constraint::consCounter = 0;

Constraint::Constraint(void)
	:Row(), vertex(Vertex()), s(Vertex()), t(Vertex()), arc(Arc())
{
	id = ++consCounter;
	consClass = -1;
	type = Constraint::C_UNKNOWN;

	return;
}

Constraint::Constraint(int _maxnnz, ROWSENSE _sense, double _rhs)
	:Row(_maxnnz, _sense, _rhs), vertex(Vertex()), s(Vertex()), t(Vertex()), arc(Arc())
{
	id = ++consCounter;
	consClass = -1;
	type = Constraint::C_UNKNOWN;

	return;
}

Constraint::Constraint(const Constraint & _constraint)
{
	++consCounter;
	
	type = _constraint.type;
	arc  = _constraint.arc;
	vertex = _constraint.vertex;
	s = _constraint.s;
	t = _constraint.t;
	consClass = _constraint.consClass;
	vertexSet = _constraint.vertexSet;
	
	return;
}

Constraint::~Constraint(void)
{
	return;
}

std::string Constraint::toString(void) const
{
	std::string consName;

	switch (type)
	{
	case Constraint::C_ASSIGNMENT:			consName = "C_ASSIGN";		break;
	case Constraint::C_FACILITY_ACTIVATION: consName = "C_F_ACTIV";		break;
	case Constraint::C_CONNECTEDNESS:		consName = "C_CONNECT";		break;
	case Constraint::C_TREE:				consName = "C_TREE";		break;
	case Constraint::C_EDGE_ACTIVATION:		consName = "C_E_ACTIV";		break;
	case Constraint::C_FLOW_CONSERVATION:	consName = "C_FLOW";		break;
	case Constraint::C_CAPACITY:			consName = "C_CAPAC";		break;
	case Constraint::C_DUAL_FEASIBILITY:	consName = "C_D_FEAS";		break;
	default:								consName = "C_UNKWN";		break;
	}

	if (consClass != -1)
	{
		std::stringstream aux;
		aux << consClass;

		consName += "_";
		consName += aux.str();
	}

	if (vertex.isValidVertex())
	{
		consName += "_";
		consName += vertex.toString();
	}

	if (s.isValidVertex())
	{
		consName += "_";
		consName += s.toString();
	}

	if (t.isValidVertex())
	{
		consName += "_";
		consName += t.toString();
	}

	if (arc.isValidArc())
	{
		consName += "_";
		consName += arc.toString();
	}

	if (vertexSet.begin() != vertexSet.end())
	{
		std::string aux("_v");

		for (std::set<int>::iterator it = vertexSet.begin(); it != vertexSet.end(); ++it)
			aux += *it;

		std::locale loc;
		const std::collate<char>& collateHasher = std::use_facet<std::collate<char> > (loc);
		long myHash = collateHasher.hash(aux.data(), aux.data() + aux.length());
		
		consName += aux;
	}

	return consName;
}
