#include "../include/Variable.h"

int Variable::varCounter = 0;

Variable::Variable(void)
	:Column(), v1(Vertex()), v2(Vertex()), arc(Arc())
{
	id = ++varCounter;
	category = 0;
	type = Variable::V_UNKNOWN;

	return;
}

Variable::Variable(COLTYPE _type, double _objcoef, double _lb, double _ub)
	:Column(_type, _objcoef, _lb, _ub), v1(Vertex()), v2(Vertex()), arc(Arc())
{
	id = ++varCounter;
	category = 0;
	type = Variable::V_UNKNOWN;

	return;
}

Variable::~Variable(void)
{
}

std::string Variable::toString(void) const
{
	std::string varName;

	switch (type)
	{
		case Variable::V_Y: varName = "Y"; break;
		case Variable::V_X: varName = "X"; break;
		case Variable::V_Z: varName = "Z"; break;
		case Variable::V_W: varName = "W"; break;
		case Variable::V_UNKNOWN: varName = "UNKN"; break;
	}

	if (category != '\0')
	{
		//varName += "_";
		varName += category;
	}

	if (v1.isValidVertex())
	{
		varName += "_";
		varName += v1.toString();
	}

	if (v2.isValidVertex())
	{
		varName += "_";
		varName += v2.toString();
	}

	if (arc.isValidArc())
	{
		varName += "_";
		varName += arc.toString();
	}

	return varName;
}
