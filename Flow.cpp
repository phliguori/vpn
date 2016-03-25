#include "../include/Flow.h"

#include <sstream>

Flow::Flow(Graph * _g)
{
	g = _g;
}

int Flow::createVariables()
{
	int nVars = 0, tmp = 0;
	clock_t start, end;

#ifdef VERBOSE
	std::cout << "\n";
	std::cout << "=======================================\n";
	std::cout << "-------------  Variables  -------------\n";
#endif

	start = clock();
	tmp = createVarX();
	end = clock();
	nVars += tmp;
#ifdef VERBOSE
	std::cout << "Variables X: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	start = clock();
	tmp = createVarY();
	end = clock();
	nVars += tmp;
#ifdef VERBOSE
	std::cout << "Variables Y: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	start = clock();
	tmp = createVarW();
	end = clock();
	nVars += tmp;
#ifdef VERBOSE
	std::cout << "Variables W: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	return nVars;
}

int Flow::createConstraints()
{
	int nCons = 0, tmp = 0;
	clock_t start, end;

#ifdef VERBOSE
	std::cout << "\n";
	std::cout << "========================================\n";
	std::cout << "-------------  Constraints  -------------\n";
#endif

	start = clock();
	tmp = createConsFlowConservation();
	end = clock();
	nCons += tmp;
#ifdef VERBOSE
	std::cout << "Constraint flow conservation: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	start = clock();
	tmp = createConsCapacityDetermination();
	end = clock();
	nCons += tmp;
#ifdef VERBOSE
	std::cout << "Constraint capacity on the edges: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	start = clock();
	tmp = createConsDualFeasibility();
	end = clock();
	nCons += tmp;
#ifdef VERBOSE
	std::cout << "Constraint dual feasibility: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	return nCons;
}

std::string  Flow::printXSol()
{
	std::ostringstream tmp;

	VariableHash::iterator vit;
	for (int varType = Variable::V_X; varType != Variable::V_UNKNOWN; varType += 10)
	{
		if (vHash.find((Variable::VARTYPE)varType) != vHash.end())
		{
			for (vit = vHash.at((Variable::VARTYPE)varType).begin(); vit != vHash.at((Variable::VARTYPE)varType).end(); ++vit)
			{
				int idx = vit->second;

				//if (pType == SolverMIP::PROBTYPE::INTEGER)
				//	xSol[idx] = (int) (xSol[idx] + 0.5);

				if (xSol[idx] > SOLVER_EPS)
					tmp << (vit->first).toString() << "; " << xSol[idx] << std::endl;
			}
		}
	}

	return tmp.str();
}

void Flow::builSolutionGraph()
{
	if (s == NULL)
		s = new Graph();

	int activeRouters = 0;
	int activeArcs = 0;

	s->nVertices = g->nVertices;
	s->nTerminals = g->nTerminals;
	s->nArcs = g->nArcs;
	s->reserve(g->nVertices + 1);

	for (int i = 0; i < g->nTerminals; ++i)
	{
		Vertex terminal = g->terminals[i];

		s->addTerminal(terminal);
		s->addVertex(terminal);
	}

	for (auto & vit : vHash[Variable::V_X])
	{
		int col = vit.second;
		Variable var = vit.first;
		Arc arc = var.getArc();
		Vertex h = var.getArc().getHead();
		Vertex t = var.getArc().getTail();

		if (!h.isTerminal() && xSol[col] > SOLVER_EPS)
		{
			++activeRouters;
			h.setColor("red");
			s->vertices[h.getCode()] = h;
		}

		if (!t.isTerminal() && xSol[col] > SOLVER_EPS)
		{
			++activeRouters;
			t.setColor("red");
			s->vertices[t.getCode()] = t;
		}

		if (xSol[col] > SOLVER_EPS)
		{
			++activeArcs;
			arc.setColor("red");
			arc.setPenwidht(3.5);
		}
		else
		{
			arc.setColor("black");
			arc.setStyle("dotted");
			arc.setPenwidht(1.0);
		}
		s->addArc(arc);
	}

#ifdef DEBUG
	s->toDot();
#endif

	return;

}

void Flow::solutionToDot(std::ostream & _output) const
{
	s->toDot(_output);

	return;
}

int Flow::createVarX()
{
	int nvars = 0;
	double coeff = 1.0;
	double lb = 0.;
	double ub = 1e20;
	VariableHash::iterator vit;
	Variable::VARTYPE varType = Variable::V_X;
	Column::COLTYPE colType = Column::COLTYPE::CONTINUOUS;

	for (int i = 1; i <= g->nVertices; ++i)
	{
		for (int j = 0; j < g->adjList[i].size(); ++j)
		{
			Arc arc = g->adjList[i][j];

			Variable x(colType, coeff, lb, ub);
			x.setType(varType);
			x.setArc(arc.toEdge());

			vit = vHash[varType].find(x);
			if (vit != vHash[varType].end())
				continue;

			bool isInserted = addCol(&x);
			if (isInserted)
			{
				x.setColIdx(getNCols() - 1);
				vHash[varType][x] = x.getColIdx();
				++nvars;
			}
		}
	}

	return nvars;
}

int Flow::createVarY()
{
	int nvars = 0;
	double coeff = 0.0;
	double lb = 0.;
	double ub = 1.;
	VariableHash::iterator vit;
	Variable::VARTYPE varType = Variable::V_Y;
	Column::COLTYPE colType = Column::COLTYPE::BINARY;

	for (int sIt = 0; sIt < g->nTerminals; ++sIt)
	{
		Vertex s = g->terminals[sIt];
		for (int tIt = 0; tIt < g->nTerminals; ++tIt)
		{
			Vertex t = g->terminals[tIt];

			if (s.getCode() == t.getCode())
				continue;

			for (int i = 1; i <= g->nVertices; ++i)
			{
				for (int j = 0; j < g->adjList[i].size(); ++j)
				{
					Arc arc = g->adjList[i][j];

					Variable y(colType, coeff, lb, ub);
					y.setType(varType);
					y.setVertex1(s);
					y.setVertex2(t);
					y.setArc(arc);

					vit = vHash[varType].find(y);
					if (vit != vHash[varType].end())
						continue;

					bool isInserted = addCol(&y);
					if (isInserted)
					{
						y.setColIdx(getNCols() - 1);
						vHash[varType][y] = y.getColIdx();
						++nvars;
					}
				}
			}
		}
	}

	return nvars;
}

int Flow::createVarW()
{
	int nvars = 0;

	double coeff = 0.0;
	double lb = 0.;
	double ub = 1e20;
	VariableHash::iterator vit;
	Variable::VARTYPE varType = Variable::V_W;
	Column::COLTYPE colType = Column::COLTYPE::CONTINUOUS;

	for (int sIt = 0; sIt < g->nTerminals; ++sIt)
	{
		Vertex s = g->terminals[sIt];

		for (int i = 1; i <= g->nVertices; ++i)
		{
			for (int j = 0; j < g->adjList[i].size(); ++j)
			{
				Arc arc = g->adjList[i][j];

				Variable w(colType, coeff, lb, ub);
				w.setType(varType);
				w.setCategory('p');
				w.setVertex1(s);
				w.setArc(arc.toEdge());

				vit = vHash[varType].find(w);
				if (vit != vHash[varType].end())
					continue;

				bool isInserted = addCol(&w);
				if (isInserted)
				{
					w.setColIdx(getNCols() - 1);
					vHash[varType][w] = w.getColIdx();
					++nvars;
				}

				w.setCategory('m');

				vit = vHash[varType].find(w);
				if (vit != vHash[varType].end())
					continue;

				isInserted = addCol(&w);
				if (isInserted)
				{
					w.setColIdx(getNCols() - 1);
					vHash[varType][w] = w.getColIdx();
					++nvars;
				}
			}
		}
	}

	return nvars;
}

int Flow::createConsFlowConservation()
{
	int ncons = 0;
	VariableHash::iterator vit;
	Row::ROWSENSE sense = Row::EQUAL;
	int maxnnz = 0;
	double rhs = 0.0;

	for (int iIt = 1; iIt <= g->nVertices; ++iIt)
	{
		Vertex i = g->vertices[iIt];

		for (int sIt = 0; sIt < g->nTerminals; ++sIt)
		{
			Vertex s = g->terminals[sIt];

			for (int tIt = 0; tIt < g->nTerminals; ++tIt)
			{
				Vertex t = g->terminals[tIt];

				if (s.getCode() == t.getCode())
					continue;

				if (i.getCode() == s.getCode())
					rhs = 1.0;
				else if (i.getCode() == t.getCode())
					rhs = -1.0;
				else
					rhs = 0.0;

				maxnnz = (2 * g->adjList[iIt].size()) + 1;

				Constraint c(maxnnz, sense, rhs);
				c.setType(Constraint::C_FLOW_CONSERVATION);
				c.setNode(i);
				c.setS(s);
				c.setT(t);

				for (int j = 0; j < g->adjList[iIt].size(); ++j)
				{
					Arc arc = g->adjList[iIt][j];

					Variable y;
					y.setType(Variable::VARTYPE::V_Y);
					y.setVertex1(s);
					y.setVertex2(t);
					y.setArc(arc);

					vit = vHash[Variable::VARTYPE::V_Y].find(y);
					if (vit != vHash[Variable::V_Y].end())
					{
						int col = vit->second;
						c.rowAddVar(col, 1.0);
					}

					y.setArc(arc.reverse());

					vit = vHash[Variable::VARTYPE::V_Y].find(y);
					if (vit != vHash[Variable::V_Y].end())
					{
						int col = vit->second;
						c.rowAddVar(col, -1.0);
					}
				}

				if (c.getRowNnz() > 0)
				{
					bool isInserted = addRow(&c);

					if (isInserted)
					{
						c.setRowIdx(getNRows() - 1);
						cHash[c] = c.getRowIdx();
						ncons++;
					}
				}
			}
		}
	}

	return ncons;
}

int Flow::createConsCapacityDetermination()
{
	int ncons = 0;
	VariableHash::iterator vit;
	Row::ROWSENSE sense = Row::LESS;
	int maxnnz = (2 * g->nTerminals) + 1;
	double rhs = 0.0;

	for (int i = 1; i <= g->nVertices; ++i)
	{
		for (int j = 0; j < g->adjList[i].size(); ++j)
		{
			Arc arc = g->adjList[i][j];

			Constraint c(maxnnz, sense, rhs);
			c.setType(Constraint::C_CAPACITY);
			c.setArc(arc.toEdge());

			for (int sIt = 0; sIt < g->nTerminals; ++sIt)
			{
				Vertex s = g->terminals[sIt];
				double bPlus = s.getEgree();
				double bMinus = s.getIngree();

				Variable w;
				w.setType(Variable::VARTYPE::V_W);
				w.setCategory('p');
				w.setVertex1(s);
				w.setArc(arc.toEdge());

				vit = vHash[Variable::VARTYPE::V_W].find(w);
				if (vit != vHash[Variable::V_W].end())
				{
					int col = vit->second;
					c.rowAddVar(col, bPlus);
				}

				w.setCategory('m');
				vit = vHash[Variable::VARTYPE::V_W].find(w);
				if (vit != vHash[Variable::V_W].end())
				{
					int col = vit->second;
					c.rowAddVar(col, bMinus);
				}
			}

			Variable x;
			x.setType(Variable::VARTYPE::V_X);
			x.setArc(arc.toEdge());

			vit = vHash[Variable::VARTYPE::V_X].find(x);
			if (vit != vHash[Variable::V_X].end())
			{
				int col = vit->second;
				c.rowAddVar(col, -1.0);
			}

			if (c.getRowNnz() > 0)
			{
				bool isInserted = addRow(&c);

				if (isInserted)
				{
					c.setRowIdx(getNRows() - 1);
					cHash[c] = c.getRowIdx();
					ncons++;
				}
			}
		}
	}

	return ncons;
}

int Flow::createConsDualFeasibility()
{
	int ncons = 0;
	VariableHash::iterator vit;
	Row::ROWSENSE sense = Row::GREATER;
	int maxnnz = 4;
	double rhs = 0.0;

	for (int sIt = 0; sIt < g->nTerminals; ++sIt)
	{
		Vertex s = g->terminals[sIt];

		for (int tIt = 0; tIt < g->nTerminals; ++tIt)
		{
			Vertex t = g->terminals[tIt];
			for (int i = 1; i <= g->nVertices; ++i)
			{
				for (int j = 0; j < g->adjList[i].size(); ++j)
				{
					Arc arc = g->adjList[i][j];

					Constraint c(maxnnz, sense, rhs);
					c.setType(Constraint::C_DUAL_FEASIBILITY);
					c.setS(s);
					c.setT(t);
					c.setArc(arc.toEdge());

					Variable wp;
					wp.setType(Variable::VARTYPE::V_W);
					wp.setCategory('p');
					wp.setVertex1(s);
					wp.setArc(arc.toEdge());

					vit = vHash[Variable::VARTYPE::V_W].find(wp);
					if (vit != vHash[Variable::V_W].end())
					{
						int col = vit->second;
						c.rowAddVar(col, 1.0);
					}

					Variable wm;
					wm.setType(Variable::VARTYPE::V_W);
					wm.setCategory('m');
					wm.setVertex1(s);
					wm.setArc(arc.toEdge());

					vit = vHash[Variable::VARTYPE::V_W].find(wm);
					if (vit != vHash[Variable::V_W].end())
					{
						int col = vit->second;
						c.rowAddVar(col, 1.0);
					}

					Variable y;
					y.setType(Variable::VARTYPE::V_Y);
					y.setVertex1(s);
					y.setVertex2(t);
					y.setArc(arc);

					vit = vHash[Variable::VARTYPE::V_Y].find(y);
					if (vit != vHash[Variable::V_Y].end())
					{
						int col = vit->second;
						c.rowAddVar(col, -1.0);
					}

					y.setArc(arc.reverse());

					vit = vHash[Variable::VARTYPE::V_Y].find(y);
					if (vit != vHash[Variable::V_Y].end())
					{
						int col = vit->second;
						c.rowAddVar(col, -1.0);
					}

					if (c.getRowNnz() > 0)
					{
						bool isInserted = addRow(&c);

						if (isInserted)
						{
							c.setRowIdx(getNRows() - 1);
							cHash[c] = c.getRowIdx();
							ncons++;
						}
					}
				}
			}
		}
	}

	return ncons;
}
