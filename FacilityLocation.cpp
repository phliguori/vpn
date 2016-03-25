#include "../include/FacilityLocation.h"

#include <set>
#include <list>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>

FacilityLocation::FacilityLocation(Graph * _g)
	:SolverMIP()
{
	g = _g;
	s = NULL;
	activeArcs = 0;
	activeRouters = 0;

	return;
}

FacilityLocation::~FacilityLocation()
{
	if (s != NULL)
		delete s;

	return;
}

void FacilityLocation::setSeparationRoutine(void * _cbdata, Separator * _separator)
{
	shutDownAutomaticCuts();

	shutDownAutomaticPreSolve();

	shutDownHeuristics();

	SolverMIP::setSeparationRoutine(_cbdata, _separator);

	return;
}

SolverMIP::SOLVERSTAT FacilityLocation::solveLRExtentedFormulations(MulltipleCutSetSeparation _cutSetSepFunc, double _tol)
{
	int pos = 0;
	int nbCuts = 0;
	int sentinel = 0;
	int MAX_ITER = 100;
	int tailOffCounter = 0;
	int tailOffTol = 3;
	int printInterval = 5;
	double currentLp = 0., lastLp = 0.;
	bool noViolatedCutFound = true;
	std::set<long> hashTable;
	TypeVariableHashPtr vHashPtr = &(vHash);
	std::vector<std::set<int>*> cutSets;
	SolverMIP::SOLVERSTAT ret = SolverMIP::SOLVERSTAT::SOLVERSTAT_UNKNOWN;

	if (xSol != NULL)
		delete[] xSol;

	xSol = new double[getNCols()];
	
	clock_t start = clock();
	do
	{
		lastLp = currentLp;
		status = SolverMIP::solve(SolverMIP::METHOD::METHOD_DUAL);

		if (status == SolverMIP::SOLVERSTAT::SOLVERSTAT_MIPOPTIMAL || status == SolverMIP::SOLVERSTAT::SOLVERSTAT_LPOPTIMAL ||
			status == SolverMIP::SOLVERSTAT::SOLVERSTAT_FEASIBLE)
		{
			ret = SolverMIP::SOLVERSTAT::SOLVERSTAT_FEASIBLE;
			currentLp = getObjVal();

			// Getting fractional node solution
			getX();

			if (sentinel % printInterval == 0)
			{
				printf("\n---- iter: %d\n", sentinel + 1);
				printf("OBJ_VAL = %lf\n", currentLp);
			}

#ifndef DEBUG
			//Printing xNode solution
			printf("\n\n---- iter: %d\n", sentinel + 1);
			for (int varType = Variable::V_X; varType != Variable::V_UNKNOWN; varType += 10)
			{
				VariableHash::iterator vit = vHashPtr->at((Variable::VARTYPE)varType).begin();
				for (; vit != vHashPtr->at((Variable::VARTYPE)varType).end(); ++vit)
				{
					int idx = vit->second;

					if (xSol[idx] > SOLVER_EPS)
						std::cout << (vit->first).toString() << "(" << idx << "); " << xSol[idx] << std::endl;
					printf("");
				}
				printf("");
			}
#endif
			// Verifying optimality conditions
			if (fabs(currentLp - lastLp) < _tol)
			{
				++tailOffCounter;
				if (tailOffCounter > tailOffTol)
				{
					ret = SolverMIP::SOLVERSTAT::SOLVERSTAT_LPOPTIMAL;
					break;
				}
			}
			else
				tailOffCounter = 0;

			// Calling the separation routine
			pos = cutSets.size();
			int cutSize = _cutSetSepFunc(g, vHashPtr, xSol, cutSets, true);

			// If a fractional cycle is found...
			if (cutSets.size() - pos > 0)
			{
				noViolatedCutFound = false;

				for (int i = pos; i < cutSets.size(); ++i)
				{
					std::set<int>* sPtr = cutSets[i];
					// Check whether the cut has already been generated
					unsigned long hashVal = hashFunc(*sPtr);
					std::set<long>::iterator it = hashTable.find(hashVal);
					if (it != hashTable.end())
					{
#ifdef DEBUG
						int warnCode = 990;
						std::string aux = convertSetToString(s);
						std::string msg = "The identified cut set was already separated " + convertSetToString(s);
						warningMsg(NULL, __func__, msg.data(), warnCode);
#endif
					}
					else
						hashTable.insert(hashVal);

					// ... we must find the cut set ...
					std::vector<Variable> cutSet;
					for (VariableHash::iterator vit = vHashPtr->at(Variable::V_Z).begin(); vit != vHashPtr->at(Variable::V_Z).end(); ++vit)
					{
						Variable z = vit->first;
						int head = z.getArc().getHead().getCode();
						int tail = z.getArc().getTail().getCode();

						std::set<int>::iterator hIt = sPtr->find(head);
						std::set<int>::iterator tIt = sPtr->find(tail);

						// ... which is composed by those arcs with one endpoint in S
						bool isHeadInS = hIt != sPtr->end();
						bool isTailInS = tIt != sPtr->end();
						if (!(isHeadInS && isTailInS) && (isHeadInS || isTailInS))
						{
							cutSet.push_back(z);
						}
					}

					// Identifying the y-variables involved in cut set constraints
					// And split them into two sets
					std::vector<Variable> sVec;
					std::vector<Variable> sCompVec;
					for (VariableHash::iterator vit = vHashPtr->at(Variable::V_Y).begin(); vit != vHashPtr->at(Variable::V_Y).end(); ++vit)
					{
						Variable v = vit->first;
						int nodeIdx = v.getVertex1().getCode();

						if (sPtr->find(nodeIdx) != sPtr->end())
						{
							sVec.push_back(v);
						}
						else
						{
							sCompVec.push_back(v);
						}
					}

					// Translating valid inequalities found into cplex/matrix representation
					int nzcnt = cutSet.size() + 2;
					std::vector<int> idx(nzcnt);
					std::vector<double> val(nzcnt);
					for (int i = 0; i < sVec.size(); ++i)
					{
						idx[0] = sVec[i].getColIdx();
						val[0] = -1.0;
						for (int j = 0; j < sCompVec.size(); ++j)
						{
							idx[1] = sCompVec[j].getColIdx();
							val[1] = -1.0;

							for (int k = 0; k < cutSet.size(); ++k)
							{
								idx[k + 2] = cutSet[k].getColIdx();
								val[k + 2] = 1.0;
							}

							// Adding user generated cut
							int nRows = 1;
							double rhs = -1;
							char sense = 'G';
							int rmatbeg = 0;
							int newColsAdded = 0;
							status = CPXaddrows(env, lp, newColsAdded, nRows, nzcnt, &rhs, &sense, &rmatbeg, &idx[0], &val[0], NULL, NULL);
							//status = CPXcutcallbackadd(_env, _cbdata, _wherefrom, nzcnt, -1, 'G', &idx[0], &val[0], CPX_USECUT_FORCE);
							if (status)
							{
								int warnCode = 999;
								std::string msg = "Failed to add integer cut.";
								warningMsg(NULL, __func__, msg.data(), warnCode);
							}
							else
								nbCuts++;

							printf("");
						}
					}
				}
#ifdef DEBUG
					// salva um arquivo .lp com o LP atual
					writeProbLP(".\\lpRelax");
#endif
			}
			else
			{
				// No violated cut was found
				noViolatedCutFound = true;
			}
		
			// If no violated cut was found
			if (noViolatedCutFound)
			{
				ret = SolverMIP::SOLVERSTAT::SOLVERSTAT_LPOPTIMAL;
				break;
			}
		}
		else
		{
			int warnCode = 201;
			std::string msg = "Model is infeasible";
			warningMsg(typeid(*this).name(), __func__, msg.data(), warnCode);

			// salva um arquivo .lp com o LP atual
			writeProbLP(".\\infeasible");

			ret = SolverMIP::SOLVERSTAT::SOLVERSTAT_INFEASIBLE;
			break;
		}

	} while (++sentinel < MAX_ITER);

	// Deallocating memory
	for (int i = 0; i < cutSets.size(); ++i)
	{
		if (cutSets[i] != NULL)
		{
			delete cutSets[i];
			cutSets[i] = NULL;
		}
	}

	clock_t end = clock();
	printf("\n-----");
	printf("\n----- iter: %d", sentinel + 1);
	printf("\n----- OBJ_VAL = %lf", currentLp);
	printf("\n----- Exectution time: %.4f", (end - start) / (double)CLOCKS_PER_SEC);

	return ret;
}

void findShortestPath(Vertex & _root, Vertex & _destination, Graph* _g, std::vector<Arc> & _path)
{
	int root = _root.getCode();
	bool isVertexFound = false;
	std::list<int> f;
	std::vector<bool> visited(_g->nVertices + 1, false);
	std::vector<Arc> pred(_g->nVertices + 1);


	// BFS from root vertex towards destination vertex
	f.push_back(root);
	while (!isVertexFound && f.size() != 0)
	{
		int v1 = *f.begin();
		visited[v1] = true;

		for (int j = 0; j < _g->adjList[v1].size(); ++j)
		{
			int v2 = _g->adjList[v1][j].getTail().getCode();

			if (!visited[v2])
			{
				pred[v2] = _g->adjList[v1][j];
				f.push_back(v2);
				visited[v2] = true;

				if (v2 == _destination.getCode())
				{
					isVertexFound = true;
					break;
				}
			}
		}
		f.pop_front();
	}

	// Retrieving the path from predecessors
	int predV = -1;
	int element = _destination.getCode();
	do
	{
		predV = pred[element].getHead().getCode();
		_path.push_back(pred[element]);
		element = predV;
	} while (element != _root.getCode() && element != -1);

	if (element == -1)
	{
		int warnCode = 970;
		std::string msg = "Problem retrieving the shortest path";
		warningMsg(NULL, __func__, msg.data(), warnCode);

		_path.clear();
	}

	printf("");
	return;
}

void FacilityLocation::builSolutionGraph()
{
	if (s == NULL)
		s = new Graph();

	activeRouters = 0;
	activeArcs = 0;

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

	for (auto& vit : vHash[Variable::V_Y])
	{
		int col = vit.second;
		Variable var = vit.first;
		Vertex v = var.getVertex1();
		int code = v.getCode();


		if (!v.isTerminal() && fabs(xSol[col] - 1.0) < SOLVER_EPS)
		{
			++activeRouters;
			v.setColor("red");
			s->vertices[code] = v;
		}
	}
	printf("");

	for (auto & vit : vHash[Variable::V_X])
	{
		Arc arc = vit.first.getArc();
		int col = vit.second;

		if (fabs(xSol[col] - 1.0) < SOLVER_EPS)
		{
			Vertex router = arc.getHead();
			Vertex terminal = arc.getTail();

			auto & element = g->arcSet[router.getCode()].find(arc);
			if (element != g->arcSet[router.getCode()].end())
			{
				arc.setColor("red");
				arc.setPenwidht(3.5);

				s->addArc(arc);
			}
			else
			{
				std::vector<Arc> path;
				findShortestPath(router, terminal, g, path);

				for (int i = 0; i < path.size(); ++i)
				{
					Arc aux = path[i];
					aux.setColor("red");
					aux.setPenwidht(1.5);

					s->addArc(aux);
					s->vertices[aux.getHead().getCode()].setColor("red");
				}
				activeArcs += path.size();
				printf("");
			}
		}

	}
	printf("");

	for (auto& vit : vHash[Variable::V_Z])
	{
		Variable v = vit.first;
		int col = vit.second;
		Arc arc = v.getArc();

		if (fabs(xSol[col] - 1.0) < SOLVER_EPS)
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
	printf("");
#ifdef DEBUG
	s->toDot();
#endif

	return;
}

int FacilityLocation::createVariables()
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
	std::cout << "Variables X: " << tmp << " (" << (end - start) / (double) CLOCKS_PER_SEC << ")\n";
#endif

	start = clock();
	tmp = createVarY();
	end = clock();
	nVars += tmp;
#ifdef VERBOSE
	std::cout << "Variables Y: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	start = clock();
	tmp = createVarZ();
	end = clock();
	nVars += tmp;
#ifdef VERBOSE
	std::cout << "Variables Z: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	return nVars;
}

int FacilityLocation::createConstraints()
{
	int nCons = 0, tmp = 0;
	clock_t start, end;

#ifdef VERBOSE
	std::cout << "\n";
	std::cout << "========================================\n";
	std::cout << "-------------  Constraints  -------------\n";
#endif

	start = clock();
	tmp = createConsAssignment();
	end = clock();
	nCons += tmp;
#ifdef VERBOSE
	std::cout << "Constraint assignment: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	start = clock();
	tmp = createConsFacilityActivation();
	end = clock();
	nCons += tmp;
#ifdef VERBOSE
	std::cout << "Constraint facility activation: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

//	start = clock();
//	tmp = createConsEdgeActivation();
//	end = clock();
//	nCons += tmp;
//#ifdef VERBOSE
//	std::cout << "Constraint edge activation: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
//#endif

	start = clock();
	tmp = createConsTree();
	end = clock();
	nCons += tmp;
#ifdef VERBOSE
	std::cout << "Constraint trees: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	start = clock();
	tmp = createConsRadiusDistance();
	end = clock();
	nCons += tmp;
#ifdef VERBOSE
	std::cout << "Constraint trees: " << tmp << " (" << (end - start) / (double)CLOCKS_PER_SEC << ")\n";
#endif

	return nCons;
}

void bfs(Graph* _g, std::vector<std::vector<int> >& _dist)
{
	//BFS
	for (int i = 1; i <= _g->nVertices; ++i)
	{
		int root = _g->vertices[i].getCode();

		std::list<int> f;
		std::vector<bool> visited(_g->nVertices + 1, false);
#ifdef DEBUG
		std::vector<int> distMat(_g->nVertices + 1, 0);
#endif
				
		f.push_back(root);
		while (f.size() != 0)
		{
			int v1 = *f.begin();
			visited[v1] = true;

			for (int j = 0; j < _g->adjList[v1].size(); ++j)
			{
				int v2 = _g->adjList[v1][j].getTail().getCode();

				if (!visited[v2])
				{
					visited[v2] = true;
					_dist[root][v2] = _dist[root][v1] + 1;
					f.push_back(v2);					
#ifdef DEBUG
					distMat[v2] = distMat[v1] + 1;
#endif
				}
			}

			f.pop_front();
		}
	}

	return;
}

Column::COLTYPE choseColType(int _type, double lb = 0., double ub = SOLVER_INF)
{
	Column::COLTYPE ret = Column::COLTYPE::UNKNOWN;

	switch (_type)
	{
	case SolverMIP::PROBTYPE::INTEGER:
		if (fabs(lb - 0.0) < SOLVER_EPS && fabs(ub - 1.0) < SOLVER_EPS)
			ret = Column::COLTYPE::BINARY;
		else
			ret = Column::COLTYPE::INTEGRAL;
		break;
	case SolverMIP::PROBTYPE::CONTINUOUS:
	default:
		ret = Column::COLTYPE::CONTINUOUS;
		break;
	}

	return ret;
}

int FacilityLocation::createVarX()
{
	int nVars = 0;
	double coeff = 0.0;
	double lb = 0.;
	double ub = 1.;
	VariableHash::iterator vit;
	Variable::VARTYPE varType = Variable::V_X;
	Column::COLTYPE colType = choseColType(pType, lb, ub);
	// @teste
	//colType = Column::CONTINUOUS;

	std::vector<std::vector<int> > dist((g->nVertices + 1), std::vector<int>((g->nVertices + 1), 0));

	clock_t start, end;

	start = clock();
	bfs(g, dist);
	end = clock();
#ifdef VERBOSE
	std::cout << "\nBFS_time: " << (end - start) / (double)CLOCKS_PER_SEC << "\n";
#endif

#ifdef DEBUG
	for (int i = 0; i <= g->nVertices; ++i)
	{
		for (int j = 0; j <= g->nVertices; ++j)
		{
			std::cout << i << " " << j << ": " << dist[i][j] << std::endl;
		}
	}
#endif

	// * Variable x_ip:
	// *
	// * Used to indicate whether the terminal vertex (client) p \in P is connected
	// * (or assigned) to the core vertex i \in I.

	double totalIngree = 0.;
	double totalEgree = 0.;
	for (int i = 1; i <= g->nVertices; ++i)
	{
		totalIngree += g->vertices[i].getIngree();
		totalEgree += g->vertices[i].getEgree();
	}

	for (int i = 1; i <= g->nVertices; ++i)
	{
		Vertex client = g->vertices[i];

		if (!client.isTerminal())
			continue;

		double tmp1 = std::min(g->vertices[i].getEgree(), totalIngree - g->vertices[i].getIngree());
		double tmp2 = std::min(g->vertices[i].getIngree(), totalEgree- g->vertices[i].getEgree());
		double b = tmp1 + tmp2;

		for (int j = 1; j <= g->nVertices; ++j)
		{
			Vertex router = g->vertices[j];
		
			// @annotation: comment next filter to resolve instances where
			// the vpn terminals may be considered internal nodes
			// 20160125
			if (client == router || router.isTerminal())
				continue;

			coeff = b * dist[client.getCode()][router.getCode()];

			Variable var(colType, coeff, lb, ub);
			var.setType(varType);
			var.setArc(Arc(router, client, 0.));

			vit = vHash[varType].find(var);
			if (vit != vHash[varType].end())
				continue;

			bool isInserted = addCol(&var);
			if (isInserted)
			{
				var.setColIdx(getNCols() - 1);
				vHash[varType][var] = var.getColIdx();
				++nVars;
			}
		}
	}

	return nVars;
}

int FacilityLocation::createVarY()
{
	int nVars = 0;
	double coeff = 0.0;
	double lb = 0.;
	double ub = 1.;
	VariableHash::iterator vit;
	Variable::VARTYPE varType = Variable::V_Y;
	Column::COLTYPE colType = choseColType(pType, lb, ub);
	// @teste
	colType = Column::CONTINUOUS;

	// * Variable y_i:
	// *
	// * Used to indicate whether the vertex i \in I is a core vertex (Steiner
	// * vertex) in the solution.
	for (int i = 1; i <= g->nVertices; ++i)
	{
		Vertex* v = &(g->vertices[i]);

		// @annotation: comment next filter to resolve instances where
		// the vpn terminals may be considered internal nodes
		// 20160125
		if (v->isTerminal())
			continue;

		Variable var(colType, coeff, lb, ub);
		var.setType(varType);
		var.setVertex1(*v);

		vit = vHash[varType].find(var);
		if (vit != vHash[varType].end())
			continue;

		bool isInserted = addCol(&var);
		if (isInserted)
		{
			var.setColIdx(getNCols() - 1);
			vHash[varType][var] = var.getColIdx();
			++nVars;
		}
	}

	return nVars;
}

int FacilityLocation::createVarZ()
{
	int nVars = 0;
	double coeff = 0.0;
	double lb = 0.;
	double ub = 1.;
	VariableHash::iterator vit;
	Variable::VARTYPE varType = Variable::V_Z;
	Column::COLTYPE colType = choseColType(pType, lb, ub);
	// @teste
	//colType = Column::CONTINUOUS;
	
	double bMinus = 0.;
	double bPlus = 0.;
	for (int i = 1; i <= g->nVertices; ++i)
	{
		bMinus += g->vertices[i].getIngree();
		bPlus += g->vertices[i].getEgree();
	}
	
	(bMinus < bPlus) ? coeff = bMinus : coeff = bPlus;

	// * Variable z_e:
	// *
	// * Used to indicate whether the edges e = (i, j) is used to connect the
	// * two core vertex i \in I and j \in I
	for (int i = 1; i <= g->nVertices; ++i)
	{
		Vertex* v = &(g->vertices[i]);

		// @annotation: comment next filter to resolve instances where
		// the vpn terminals may be considered internal nodes
		// 20160125
		if (v->isTerminal())
			continue;

		for (int j = 0; j < g->adjList[i].size(); ++j)
		{
			Vertex* u = &(g->adjList[i][j].getTail());
			Arc* arc = &(g->adjList[i][j]);

			// @annotation: comment next filter to resolve instances where
			// the vpn terminals may be considered internal nodes
			// 20160125
			if (arc->getHead().isTerminal() || arc->getTail().isTerminal())
				continue;

			Variable var(colType, coeff, lb, ub);
			var.setType(varType);
			var.setArc(arc->toEdge());

			vit = vHash[varType].find(var);
			if (vit != vHash[varType].end())
				continue;

			bool isInserted = addCol(&var);
			if (isInserted)
			{
				var.setColIdx(getNCols() - 1);
				vHash[varType][var] = var.getColIdx();
				++nVars;
			}
		}
	}

	return nVars;
}

int FacilityLocation::createConsAssignment()
{
	int numCons = 0;
	VariableHash::iterator vit;
	Row::ROWSENSE sense = Row::EQUAL;
	int maxnnz = (g->nVertices - g->nTerminals) + 1;
	double rhs = 1.0;

	for (int i = 1; i <= g->nVertices; ++i)
	{
		Vertex client = g->vertices[i];

		if (!client.isTerminal())
			continue;

		Constraint c(maxnnz, sense, rhs);
		c.setType(Constraint::C_ASSIGNMENT);
		c.setNode(client);

		Variable x;
		x.setType(Variable::V_X);

		for (int j = 1; j <= g->nVertices; ++j)
		{
			Vertex router = g->vertices[j];

			// @annotation: comment next filter to resolve instances where
			// the vpn terminals may be considered internal nodes
			// 20160125
			if (client == router || router.isTerminal())
				continue;

			x.setArc(Arc(router, client, 0.));

			vit = vHash[Variable::V_X].find(x);
			if (vit != vHash[Variable::V_X].end())
			{
				int colVarX = vit->second;
				c.rowAddVar(colVarX, 1.0);
			}
		}

		if (c.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c);

			if (isInserted)
			{
				c.setRowIdx(getNRows() - 1);
				cHash[c] = c.getRowIdx();
				numCons++;
			}
		}
	}

	return numCons;
}

int FacilityLocation::createConsFacilityActivation()
{
	int numCons = 0;
	VariableHash::iterator vit;
	Row::ROWSENSE sense = Row::LESS;
	int maxnnz = 2;
	double rhs = 0.0;
	int colVarX, colVarY, colVarZ;

	for (vit = vHash[Variable::V_X].begin(); vit != vHash[Variable::V_X].end(); ++vit)
	{
		colVarX = vit->second;
		Arc arc = vit->first.getArc();
		Vertex router = vit->first.getArc().getHead();
		Vertex client = vit->first.getArc().getTail();

		Constraint c(maxnnz, sense, rhs);
		c.setType(Constraint::C_FACILITY_ACTIVATION);
		c.setClass(1);
		c.setArc(arc);

		Variable y;
		y.setType(Variable::V_Y);
		y.setVertex1(router);
		VariableHash::iterator aux = vHash[Variable::V_Y].find(y);
		if (aux != vHash[Variable::V_Y].end())
		{
			colVarY = aux->second;
			
			c.rowAddVar(colVarX, 1.0);
			c.rowAddVar(colVarY, -1.0);
		}

		if (c.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c);

			if (isInserted)
			{
				c.setRowIdx(getNRows() - 1);
				cHash[c] = c.getRowIdx();
				numCons++;
			}
		}
	}

	for (vit = vHash[Variable::V_Z].begin(); vit != vHash[Variable::V_Z].end(); ++vit)
	{
		colVarZ = vit->second;
		Arc arc = vit->first.getArc();
		Vertex head = vit->first.getArc().getHead();
		Vertex tail = vit->first.getArc().getTail();

		Constraint c1(maxnnz, sense, rhs);
		c1.setType(Constraint::C_FACILITY_ACTIVATION);
		c1.setClass(2);
		c1.setArc(arc);

		Variable yHead;
		yHead.setType(Variable::V_Y);
		yHead.setVertex1(head);
		VariableHash::iterator aux = vHash[Variable::V_Y].find(yHead);
		if (aux != vHash[Variable::V_Y].end())
		{
			colVarY = aux->second;

			c1.rowAddVar(colVarZ, 1.0);
			c1.rowAddVar(colVarY, -1.0);
		}

		if (c1.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c1);

			if (isInserted)
			{
				c1.setRowIdx(getNRows() - 1);
				cHash[c1] = c1.getRowIdx();
				numCons++;
			}
		}

		Constraint c2(maxnnz, sense, rhs);
		c2.setType(Constraint::C_FACILITY_ACTIVATION);
		c2.setClass(2);
		c2.setArc(arc.reverse());

		Variable yTail;
		yTail.setType(Variable::V_Y);
		yTail.setVertex1(tail);
		aux = vHash[Variable::V_Y].find(yTail);
		if (aux != vHash[Variable::V_Y].end())
		{
			colVarY = aux->second;

			c2.rowAddVar(colVarZ, 1.0);
			c2.rowAddVar(colVarY, -1.0);
		}

		if (c2.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c2);

			if (isInserted)
			{
				c2.setRowIdx(getNRows() - 1);
				cHash[c2] = c2.getRowIdx();
				numCons++;
			}
		}
	}
	
	return numCons;
}

int FacilityLocation::createConsFacilityActivation(int)
{
	int numCons = 0;
	VariableHash::iterator vit;
	Row::ROWSENSE sense = Row::LESS;
	int maxnnz = 2;
	double rhs = 0.0;
	int colVarX, colVarY, colVarZ;

	for (vit = vHash[Variable::V_X].begin(); vit != vHash[Variable::V_X].end(); ++vit)
	{
		colVarX = vit->second;
		Arc arc = vit->first.getArc();
		Vertex router = vit->first.getArc().getHead();
		Vertex client = vit->first.getArc().getTail();

		Constraint c(maxnnz, sense, rhs);
		c.setType(Constraint::C_FACILITY_ACTIVATION);
		c.setClass(1);
		c.setArc(arc);

		Variable y;
		y.setType(Variable::V_Y);
		y.setVertex1(router);
		VariableHash::iterator aux = vHash[Variable::V_Y].find(y);
		if (aux != vHash[Variable::V_Y].end())
		{
			colVarY = aux->second;

			c.rowAddVar(colVarX, 1.0);
			c.rowAddVar(colVarY, -1.0);
		}

		if (c.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c);

			if (isInserted)
			{
				c.setRowIdx(getNRows() - 1);
				cHash[c] = c.getRowIdx();
				numCons++;
			}
		}
	}

	for (vit = vHash[Variable::V_Z].begin(); vit != vHash[Variable::V_Z].end(); ++vit)
	{
		colVarZ = vit->second;
		Arc arc = vit->first.getArc();
		Vertex head = vit->first.getArc().getHead();
		Vertex tail = vit->first.getArc().getTail();

		Constraint c1(maxnnz, sense, rhs);
		c1.setType(Constraint::C_FACILITY_ACTIVATION);
		c1.setClass(2);
		c1.setArc(arc);

		Variable yHead;
		yHead.setType(Variable::V_Y);
		yHead.setVertex1(head);
		VariableHash::iterator aux = vHash[Variable::V_Y].find(yHead);
		if (aux != vHash[Variable::V_Y].end())
		{
			colVarY = aux->second;

			c1.rowAddVar(colVarZ, 1.0);
			c1.rowAddVar(colVarY, -1.0);
		}

		if (c1.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c1);

			if (isInserted)
			{
				c1.setRowIdx(getNRows() - 1);
				cHash[c1] = c1.getRowIdx();
				numCons++;
			}
		}

		Constraint c2(maxnnz, sense, rhs);
		c2.setType(Constraint::C_FACILITY_ACTIVATION);
		c2.setClass(2);
		c2.setArc(arc.reverse());

		Variable yTail;
		yTail.setType(Variable::V_Y);
		yTail.setVertex1(tail);
		aux = vHash[Variable::V_Y].find(yTail);
		if (aux != vHash[Variable::V_Y].end())
		{
			colVarY = aux->second;

			c2.rowAddVar(colVarZ, 1.0);
			c2.rowAddVar(colVarY, -1.0);
		}

		if (c2.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c2);

			if (isInserted)
			{
				c2.setRowIdx(getNRows() - 1);
				cHash[c2] = c2.getRowIdx();
				numCons++;
			}
		}
	}

	return numCons;
}

int FacilityLocation::createConsTree()
{
	int numCons = 0;
	VariableHash::iterator vit;
	int maxnnz = vHash[Variable::V_Y].size() + vHash[Variable::V_Z].size();
	Row::ROWSENSE sense = Row::EQUAL;
	double rhs = -1;
	
	Constraint c(maxnnz, sense, rhs);
	c.setType(Constraint::C_TREE);

	if (cHash.find(c) == cHash.end())
	{
		for (vit = vHash[Variable::V_Z].begin(); vit != vHash[Variable::V_Z].end(); ++vit)
		{
			int colidx = vit->second;
			c.rowAddVar(colidx, 1.0);
		}

		for (vit = vHash[Variable::V_Y].begin(); vit != vHash[Variable::V_Y].end(); ++vit)
		{
			int colidx = vit->second;
			c.rowAddVar(colidx, -1.0);
		}

		if (c.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c);

			if (isInserted)
			{
				c.setRowIdx(getNRows() - 1);
				cHash[c] = c.getRowIdx();
				numCons++;
			}
		}
	}

	return numCons;
}

int FacilityLocation::createConsRadiusDistance()
{
	// TODO
	// Implementação do corte a seguir: os vértices roteadores da árvore VPN não poderão ser
	// folhas na solução. Isso implica que deverá existir pelo menos dois arcos com extremidada
	// em cada vértice roteador:
	//
	// 2 y_i \leq \sum\limits_{e \in \delta(i)} z_e + \sum\limits_{p \in P} x_{ip}, \forall i \in I
	int numCons = 0;

	std::vector<std::vector<int> > dist((g->nVertices + 1), std::vector<int>((g->nVertices + 1), 0));

	clock_t start, end;

	start = clock();
	bfs(g, dist);
	end = clock();
#ifdef VERBOSE
	std::cout << "\nBFS_time: " << (end - start) / (double)CLOCKS_PER_SEC << "\n";
#endif

	for (auto & terminal : g->terminals)
	{
		// Identifyng all possible distances from node terminal
		std::set<int> distSet;
		for (auto & d : dist[terminal.getCode()])
		{
			if (d > 1)
				distSet.insert(d);
		}
		printf("");

		std::map<int, std::vector<Vertex>> auxp;
		for (auto & router : g->vertices)
		{
			// @annotation: comment next filter to resolve instances where
			// the vpn terminals may be considered internal nodes
			// 20160125
			if (terminal == router || router.isTerminal())
				continue;

			// The pair (terminal, router) corresponds to a variable x_ip
			// and now we know the distance (number of hops) between these
			// two nodes, given by d.
			int dtr = dist[terminal.getCode()][router.getCode()];

			std::vector<Arc> path;
			findShortestPath(router, terminal, g, path);

			if (path.size() > 0)
			{
				Vertex element = (*path.rbegin()).getHead();
				auxp[dtr].push_back(element);
				printf("");
			}
		}
		printf("");
		for (auto & eleDist : auxp)
		{
			int dist = eleDist.first;
			std::vector<Vertex> & elements = eleDist.second;
			std::vector<Vertex> & elementsConstrained = auxp[dist - 1];

			for (int j = 0; j < elementsConstrained.size(); ++j)
			{
				Vertex vj = elementsConstrained[j];
				Variable yj;
				yj.setType(Variable::V_Y);
				yj.setVertex1(vj);

				VariableHash::iterator yIt = vHash[Variable::V_Y].find(yj);
				if (yIt == vHash[Variable::V_Y].end())
					continue;

				Constraint c(elements.size() + 2, Row::LESS, 1.0);
				c.setType(Constraint::C_UNKNOWN);
				c.setClass(dist);
				c.setNode(vj);

				c.rowAddVar(yIt->second, 1.0);

				//std::cout << c.toString() << ": " <<  yj.toString();

				for (int i = 0; i < elements.size(); ++i)
				{
					Vertex router = elements[i];
					Arc a = Arc(router, terminal, 0.);
					Variable xrt;
					xrt.setType(Variable::V_X);
					xrt.setArc(a);

					VariableHash::iterator xIt = vHash[Variable::V_X].find(xrt);
					if (xIt == vHash[Variable::V_X].end())
						continue;

					c.rowAddVar(xIt->second, 1.0);

					//std::cout << " + " << xrt.toString();
				}
				if (c.getRowNnz() > 1)
				{
					bool isInserted = addRow(&c);

					if (isInserted)
					{
						c.setRowIdx(getNRows() - 1);
						cHash[c] = c.getRowIdx();
						numCons++;
					}
				}
				//std::cout << " <= " << 1 << std::endl;
			}
		}
		printf("");
	}
	printf("");
	return numCons;
}

int FacilityLocation::createConsEdgeActivation()
{
	// TODO
	// Implementação do corte a seguir: os vértices roteadores da árvore VPN não poderão ser
	// folhas na solução. Isso implica que deverá existir pelo menos dois arcos com extremidada
	// em cada vértice roteador:
	//
	// 2 y_i \leq \sum\limits_{e \in \delta(i)} z_e + \sum\limits_{p \in P} x_{ip}, \forall i \in I

	int numCons = 0;
	VariableHash::iterator vit;
	Row::ROWSENSE sense = Row::LESS;
	int maxnnz = vHash[Variable::V_X].size() + vHash[Variable::V_Z].size() + 1;
	double rhs = 0.0;
	int colVarX, colVarZ, colVarY;

	for (auto & varY : vHash[Variable::V_Y])
	{
		int colVarY = varY.second;
		Vertex vertex = varY.first.getVertex1();

		Constraint c(maxnnz, sense, rhs);
		c.setType(Constraint::C_EDGE_ACTIVATION);
		c.setNode(vertex);

		c.rowAddVar(colVarY, 2.0);

		for (auto & varZ : vHash[Variable::V_Z])
		{
			colVarZ = varZ.second;
			Arc arc = varZ.first.getArc();

			//if (arc.getHead() == vertex.getCode() || arc.getTail() == vertex.getCode())
			if (arc.isIncident(vertex))
			{
				c.rowAddVar(colVarZ, -1.0);
			}
		}

		for (auto & varX : vHash[Variable::V_X])
		{
			colVarX = varX.second;
			Arc arc = varX.first.getArc();

			//if (arc.getHead() == vertex.getCode() || arc.getTail() == vertex.getCode())
			if (arc.isIncident(vertex))
			{
				c.rowAddVar(colVarX, -1.0);
			}
		}

		if (c.getRowNnz() > 0)
		{
			bool isInserted = addRow(&c);

			if (isInserted)
			{
				c.setRowIdx(getNRows() - 1);
				cHash[c] = c.getRowIdx();
				numCons++;
			}
		}
	}

	return numCons;
}

std::string FacilityLocation::printXSol()
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

void FacilityLocation::solutionToDot(std::ostream & _output) const
{
	s->toDot(_output);

	return;
}

