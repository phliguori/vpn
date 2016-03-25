#include "../include/CutSetSeparation.h"

#include "../include/Graph.h"
#include "../include/Variable.h"
#include "../include/CallbackData.h"

#include <set>
#include <map>
#include <vector>
#include <utility>

//#include <lemon/bfs.h>
//#include <lemon/core.h>
//#include <lemon/hao_orlin.h>	//double lemon::Tolerance<double>::def_epsilon = 1e-6;
#include <lemon/preflow.h>		//double lemon::Tolerance<double>::def_epsilon = 1e-6;
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>

double lemon::Tolerance<double>::def_epsilon = 1e-6;	
static double EPS = 1e-6;

// Some useful typedefs
typedef lemon::ListDigraph Digraph;
typedef Digraph::Arc DArc;
typedef Digraph::Node DNode;

void printSolution(TypeVariableHashPtr _vHashPtr, double* _xNode)
{
	for (int varType = Variable::V_X; varType != Variable::V_UNKNOWN; varType += 10)
	{
		VariableHash::iterator vit = _vHashPtr->at((Variable::VARTYPE)varType).begin();
		for (; vit != _vHashPtr->at((Variable::VARTYPE)varType).end(); ++vit)
		{
			int idx = vit->second;

			if (_xNode[idx] > EPS)
				std::cout << (vit->first).toString() << "(" << idx << "); " << _xNode[idx] << std::endl;
			printf("");
		}

		printf("");
	}

	return;
}

int __stdcall integerCutSet(Graph* _g, TypeVariableHashPtr _vHashPtr, double * _x, std::set<int> & _s)
{
	//Build a directed lemon graph from the integer solution
	Digraph graphSuportLemon;
	std::map<int, DNode> nodesMapLemon;
	std::vector<DArc> arcsVecLemon;
	Digraph::NodeMap<int> connectedComponent(graphSuportLemon);

	// Building nodes list
	for (VariableHash::iterator it = _vHashPtr->at(Variable::V_Y).begin(); it != _vHashPtr->at(Variable::V_Y).end(); ++it)
	{
		int col = it->second;
		
		if (_x[col] > 0.5) // This is a integer separation
		{
			int nodeIdx = it->first.getVertex1().getCode();;
			nodesMapLemon[nodeIdx] = graphSuportLemon.addNode();
		}
	}
	printf("");

	for (VariableHash::iterator it = _vHashPtr->at(Variable::V_Z).begin(); it != _vHashPtr->at(Variable::V_Z).end(); ++it)
	{
		int col = it->second;

		if (_x[col] < 0.5) continue;

		Variable z = it->first;

		int head = z.getArc().getHead().getCode();
		int tail = z.getArc().getTail().getCode();

		//Insert arcs into lemon graph
		std::map<int, DNode>::iterator hIt = nodesMapLemon.find(head);
		std::map<int, DNode>::iterator tIt = nodesMapLemon.find(tail);

		if (hIt != nodesMapLemon.end() && tIt != nodesMapLemon.end())
		{
			DNode headNode = hIt->second;
			DNode tailNode = tIt->second;
			DArc arc = graphSuportLemon.addArc(headNode, tailNode);
			DArc rarc = graphSuportLemon.addArc(tailNode, headNode);

			arcsVecLemon.push_back(arc);
			arcsVecLemon.push_back(rarc);
		}
	}

	// Call to Lemon Strongly Connected Component Routine
	int nconnected = lemon::stronglyConnectedComponents(graphSuportLemon, connectedComponent);

	// If there is more than one connected component, so there will exist some cut arcs
	if (nconnected > 1)
	{
		// Determining the size of each connected component
		std::map<int, int> componentSize;
		for (std::map<int, DNode>::iterator it = nodesMapLemon.begin(); it != nodesMapLemon.end(); ++it)
		{
			int nodeIdx = it->first;
			DNode node = it->second;

			componentSize[connectedComponent[node]]++;
		}

		// Chosing which component will play set S
		int maxElements = 0;
		int chosenComponent = 0;
		for (auto& c : componentSize)
		{
			if (c.second > maxElements)
			{
				chosenComponent = c.first;
				maxElements = c.second;
			}
		}

		// Building set S
		for (std::map<int, DNode>::iterator it = nodesMapLemon.begin(); it != nodesMapLemon.end(); ++it)
		{
			int nodeIdx = it->first;
			DNode node = it->second;

			//if (connectedComponent[node] == chosenComponent)
			if (connectedComponent[node] == 0)
			{
				_s.insert(nodeIdx);
			}
		}

		printf("");
	}

	return nconnected;
}

int
CPXPUBLIC separateIntegerCutSet(
	CPXCENVptr _env,
	void       *_cbdata,
	int        _wherefrom,
	void       *_cbhandle,
	int        *_useraction_p)
{
#ifdef DEBUG
	printf(__func__);
	printf("\n");
#endif

	*_useraction_p = CPX_CALLBACK_DEFAULT;

	CBData* cbHandlerPtr = (CBData*)_cbhandle;
	int status = 0;
	int numcols = cbHandlerPtr->numCols;
	Graph * g = cbHandlerPtr->graph;
	TypeVariableHashPtr vHashPtr = cbHandlerPtr->vHashPtr;
	std::set<long>* hashTable = cbHandlerPtr->hashTable;

	double* xNode = &(cbHandlerPtr->nodeSolX[0]);
	//std::vector<double> xNode(numcols, 0.0);
	
	// Getting integer node solution
	status = CPXgetcallbacknodex(_env, _cbdata, _wherefrom, xNode, 0, numcols - 1);
	if (status)
	{
		int warnCode = 901;
		std::string msg = "Failed to get node solution";
		warningMsg(NULL, __func__, msg.data(), warnCode);
		return status;
	}

#ifdef DEBUG
	//Printing xNode solution
	printSolution(vHashPtr, &xNode[0]);
#endif

	// Calling the separation routine
	std::set<int> s;
	int nComponents = integerCutSet(g, vHashPtr, &xNode[0], s);

	// If a cycle is found...
	if (nComponents > 1)
	{
		// ... check whether the cut has already been generated ...
		unsigned long hashVal = hashFunc(s);
		std::set<long>::iterator it = hashTable->find(hashVal);
		if (it != hashTable->end())
		{
#ifdef DEBUG
			int warnCode = 990;
			std::string aux = convertSetToString(s);
			std::string msg = "The identified cut set was already separated " + convertSetToString(s);
			warningMsg(NULL, __func__, msg.data(), warnCode);
#endif
			// ... in the affirmative case, stop callback execution.
			return 0;
		}
		else
			hashTable->insert(hashVal);

		// ... in the negative case, we must find the cut set ...
		std::vector<Variable> cutSet;
		for (VariableHash::iterator vit = vHashPtr->at(Variable::V_Z).begin(); vit != vHashPtr->at(Variable::V_Z).end(); ++vit)
		{
			Variable z = vit->first;
			int head = z.getArc().getHead().getCode();
			int tail = z.getArc().getTail().getCode();

			std::set<int>::iterator hIt = s.find(head);
			std::set<int>::iterator tIt = s.find(tail);

			// ... which is composed by those arcs with one endpoint in S
			bool isHeadInS = hIt != s.end();
			bool isTailInS = tIt != s.end();
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

			if (s.find(nodeIdx) != s.end())
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
				status = CPXcutcallbackadd(_env, _cbdata, _wherefrom, nzcnt, -1, 'G', &idx[0], &val[0], CPX_USECUT_FORCE);
				if (status)
				{
					int warnCode = 999;
					std::string msg = "Failed to add integer cut.";
					warningMsg(NULL, __func__, msg.data(), warnCode);
				}
				else
				{
					cbHandlerPtr->totalNbCuts++;
					cbHandlerPtr->totalNbIntCuts++;
				}

				printf("");
			}
		}

		printf("");

		// Telling CPLEX that cuts have been created
		*_useraction_p = CPX_CALLBACK_SET;

#ifdef DEBUG
		// salva um arquivo .lp com o LP atual
		CPXLPptr nodelp = NULL;
		status = CPXgetcallbacknodelp(_env, _cbdata, _wherefrom, &nodelp);
		if (status)
		{
			fprintf(stderr, "Failed to get node lp pointer.\n");
			return status;
		}

		status = CPXwriteprob(_env, nodelp, ".\\lazycutset.lp", "LP");
		if (status)
		{
			fprintf(stderr, "Failed to print lp file.\n");
			return status;
		}
#endif

	}

	return status;
}

int __stdcall fractionalCutSet(Graph* _g, TypeVariableHashPtr _vHashPtr, double * _x, std::set<int> & _s)
{
	//Build a directed lemon graph from the fractional solution
	Digraph graphSuportLemon;
	std::map<int, DNode> nodesMapLemon;
	std::vector<DArc> arcsVecLemon;
	Digraph::ArcMap<double> arcsCapMap(graphSuportLemon);
	//Digraph::NodeMap<int> connectedComponent(graphSuportLemon);

	// Building nodes list
	for (VariableHash::iterator it = _vHashPtr->at(Variable::V_Y).begin(); it != _vHashPtr->at(Variable::V_Y).end(); ++it)
	{
		int col = it->second;

		if (_x[col] > EPS)
		{
			int nodeIdx = it->first.getVertex1().getCode();
			nodesMapLemon[nodeIdx] = graphSuportLemon.addNode();
		}
	}

	bool haveFractionalCandidate = false;
	double maxRhs = 0;
	std::map<int, std::map<int, double> > rhs;
	for (auto & it : _vHashPtr->at(Variable::V_Y))
	{
		std::string aux1 = it.first.toString();
		int nodei = it.first.getVertex1().getCode();
		int coli = it.second;

		if (_x[coli] > EPS)
		{
			for (auto & jt : _vHashPtr->at(Variable::V_Y))
			{
				std::string aux2 = jt.first.toString();
				int nodej = jt.first.getVertex1().getCode();
				int colj = jt.second;

				if (coli != colj && _x[colj] > EPS)
				{
					double rhsij = _x[coli] + _x[colj] - 1.0;
					if (rhsij > EPS)
					{
						haveFractionalCandidate = true;
						rhs[nodei][nodej] = _x[coli] + _x[colj] - 1.0;

						if (rhsij > maxRhs)
							maxRhs = rhsij;
					}
				}
			}
		}
	}

	if (!haveFractionalCandidate)
		return 0;

	int chosenNode;
	int maxPossibilities = 0;
	for (auto & nodeMap : rhs)
	{
		int npossibilities = nodeMap.second.size();
		if (npossibilities > maxPossibilities)
		{
			maxPossibilities = npossibilities;
			chosenNode = nodeMap.first;
		}
	}
	printf("");

	for (VariableHash::iterator it = _vHashPtr->at(Variable::V_Z).begin(); it != _vHashPtr->at(Variable::V_Z).end(); ++it)
	{
		int col = it->second;

		if (_x[col] < EPS) continue;

		Variable z = it->first;

		int head = z.getArc().getHead().getCode();
		int tail = z.getArc().getTail().getCode();

		//Insert arcs into lemon graph
		std::map<int, DNode>::iterator hIt = nodesMapLemon.find(head);
		std::map<int, DNode>::iterator tIt = nodesMapLemon.find(tail);

		if (hIt != nodesMapLemon.end() && tIt != nodesMapLemon.end())
		{
			DNode headNode = hIt->second;
			DNode tailNode = tIt->second;
			DArc arc = graphSuportLemon.addArc(headNode, tailNode);
			DArc rarc = graphSuportLemon.addArc(tailNode, headNode);

			arcsVecLemon.push_back(arc);
			arcsVecLemon.push_back(rarc);

			arcsCapMap[arc] = _x[col];
			arcsCapMap[rarc] = _x[col];
		}
	}
	printf("");

	// Initializing and calling preflow algorithm
	int add_cuts = 0;
	DNode root = nodesMapLemon[chosenNode];
	std::set<std::set<int> > myCuts;

	for (std::map<int, double>::iterator st = rhs[chosenNode].begin(); st != rhs[chosenNode].end(); ++st)
	{
		int vSink = st->first;
		double rhs = st->second;
		DNode sink = nodesMapLemon[vSink];

		Digraph::NodeMap<bool> minCut(graphSuportLemon);
		lemon::Preflow<Digraph, Digraph::ArcMap<double> > preflow(graphSuportLemon, arcsCapMap, root, sink);
		preflow.run();

		double flow = preflow.flowValue();
		if (std::fabs(flow - rhs) > EPS)
		{
			preflow.minCutMap(minCut);

			// Building set S
			for (std::map<int, DNode>::iterator it = nodesMapLemon.begin(); it != nodesMapLemon.end(); ++it)
			{
				int nodeIdx = it->first;
				DNode node = it->second;

				if (minCut[node])
				{
					_s.insert(nodeIdx);
				}
			}
			break;
		}
	}

	return _s.size();
}

int __stdcall fractionalCutSet(Graph* _g, TypeVariableHashPtr _vHashPtr, double * _x, std::vector<std::set<int>* > & _cutSets, bool _findAllCutsSets)
{
	//Build a directed lemon graph from the fractional solution
	Digraph graphSuportLemon;
	std::map<int, DNode> nodesMapLemon;
	std::vector<DArc> arcsVecLemon;
	Digraph::ArcMap<double> arcsCapMap(graphSuportLemon);
	//Digraph::NodeMap<int> connectedComponent(graphSuportLemon);

	// Building nodes list
	for (VariableHash::iterator it = _vHashPtr->at(Variable::V_Y).begin(); it != _vHashPtr->at(Variable::V_Y).end(); ++it)
	{
		int col = it->second;

		if (_x[col] > EPS)
		{
			int nodeIdx = it->first.getVertex1().getCode();
			nodesMapLemon[nodeIdx] = graphSuportLemon.addNode();
		}
	}

	bool haveFractionalCandidate = false;
	double maxRhs = 0;
	std::map<int, std::map<int, double> > rhs;
	for (auto & it : _vHashPtr->at(Variable::V_Y))
	{
		std::string aux1 = it.first.toString();
		int nodei = it.first.getVertex1().getCode();
		int coli = it.second;

		if (_x[coli] > EPS)
		{
			for (auto & jt : _vHashPtr->at(Variable::V_Y))
			{
				std::string aux2 = jt.first.toString();
				int nodej = jt.first.getVertex1().getCode();
				int colj = jt.second;

				if (coli != colj && _x[colj] > EPS)
				{
					double rhsij = _x[coli] + _x[colj] - 1.0;
					if (rhsij > EPS)
					{
						haveFractionalCandidate = true;
						rhs[nodei][nodej] = _x[coli] + _x[colj] - 1.0;

						if (rhsij > maxRhs)
							maxRhs = rhsij;
					}
				}
			}
		}
	}

	if (!haveFractionalCandidate)
		return 0;

	int chosenNode;
	int maxPossibilities = 0;
	for (auto & nodeMap : rhs)
	{
		int npossibilities = nodeMap.second.size();
		if (npossibilities > maxPossibilities)
		{
			maxPossibilities = npossibilities;
			chosenNode = nodeMap.first;
		}
	}
	printf("");

	for (VariableHash::iterator it = _vHashPtr->at(Variable::V_Z).begin(); it != _vHashPtr->at(Variable::V_Z).end(); ++it)
	{
		int col = it->second;

		if (_x[col] < EPS) continue;

		Variable z = it->first;

		int head = z.getArc().getHead().getCode();
		int tail = z.getArc().getTail().getCode();

		//Insert arcs into lemon graph
		std::map<int, DNode>::iterator hIt = nodesMapLemon.find(head);
		std::map<int, DNode>::iterator tIt = nodesMapLemon.find(tail);

		if (hIt != nodesMapLemon.end() && tIt != nodesMapLemon.end())
		{
			DNode headNode = hIt->second;
			DNode tailNode = tIt->second;
			DArc arc = graphSuportLemon.addArc(headNode, tailNode);
			DArc rarc = graphSuportLemon.addArc(tailNode, headNode);

			arcsVecLemon.push_back(arc);
			arcsVecLemon.push_back(rarc);

			arcsCapMap[arc] = _x[col];
			arcsCapMap[rarc] = _x[col];
		}
	}
	printf("");

	// Initializing and calling preflow algorithm
	int add_cuts = 0;
	DNode root = nodesMapLemon[chosenNode];
	std::set<std::set<int> > myCuts;

	for (std::map<int, double>::iterator st = rhs[chosenNode].begin(); st != rhs[chosenNode].end(); ++st)
	{
		int vSink = st->first;
		double rhsVal = st->second;
		DNode sink = nodesMapLemon[vSink];

		Digraph::NodeMap<bool> minCut(graphSuportLemon);
		lemon::Preflow<Digraph, Digraph::ArcMap<double> > preflow(graphSuportLemon, arcsCapMap, root, sink);
		preflow.run();

		double flow = preflow.flowValue();
		if (std::fabs(flow - rhsVal) > EPS)
		{
			preflow.minCutMap(minCut);

			// Building set S
			std::set<int>* s = new std::set<int>();
			for (std::map<int, DNode>::iterator it = nodesMapLemon.begin(); it != nodesMapLemon.end(); ++it)
			{
				int nodeIdx = it->first;
				DNode node = it->second;

				if (minCut[node])
				{
					s->insert(nodeIdx);
				}
			}

			_cutSets.push_back(s);

			if (!_findAllCutsSets)
				break;
		}
	}

	return _cutSets.size();
}

int
CPXPUBLIC separateFractionalCutSet(
	CPXCENVptr _env,
	void       *_cbdata,
	int        _wherefrom,
	void       *_cbhandle,
	int        *_useraction_p)
{
#ifdef DEBUG
	printf(__func__);
	printf("\n");
#endif

	*_useraction_p = CPX_CALLBACK_DEFAULT;

	CBData* cbHandlerPtr = (CBData*)_cbhandle;
	int status = 0;
	int numcols = cbHandlerPtr->numCols;
	Graph * g = cbHandlerPtr->graph;
	TypeVariableHashPtr vHashPtr = cbHandlerPtr->vHashPtr;
	std::set<long>* hashTable = cbHandlerPtr->hashTable;

	double* xNode = &(cbHandlerPtr->nodeSolX[0]);
	//std::vector<double> xNode(numcols, 0.0);


	// Getting fractional node solution
	status = CPXgetcallbacknodex(_env, _cbdata, _wherefrom, xNode, 0, numcols - 1);
	if (status)
	{
		int warnCode = 902;
		std::string msg = "Failed to get node solution";
		warningMsg(NULL, __func__, msg.data(), warnCode);
		return status;
	}

#ifdef DEBUG
	//Printing xNode solution
	printSolution(vHashPtr, &xNode[0]);
#endif

	// Calling the separation routine
	std::set<int> s;
	int cutSize = fractionalCutSet(g, vHashPtr, &xNode[0], s);
	
	// If a fractional cycle is found...
	if (cutSize > 0)
	{
		// Check whether the cut has already been generated
		unsigned long hashVal = hashFunc(s);
		std::set<long>::iterator it = hashTable->find(hashVal);
		if (it != hashTable->end())
		{
#ifdef DEBUG
			int warnCode = 990;
			std::string aux = convertSetToString(s);
			std::string msg = "The identified cut set was already separated " + convertSetToString(s);
			warningMsg(NULL, __func__, msg.data(), warnCode);
#endif
			return 0;
		}
		else
			hashTable->insert(hashVal);

		// ... we must find the cut set ...
		std::vector<Variable> cutSet;
		for (VariableHash::iterator vit = vHashPtr->at(Variable::V_Z).begin(); vit != vHashPtr->at(Variable::V_Z).end(); ++vit)
		{
			Variable z = vit->first;
			int head = z.getArc().getHead().getCode();
			int tail = z.getArc().getTail().getCode();

			std::set<int>::iterator hIt = s.find(head);
			std::set<int>::iterator tIt = s.find(tail);

			// ... which is composed by those arcs with one endpoint in S
			bool isHeadInS = hIt != s.end();
			bool isTailInS = tIt != s.end();
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

			if (s.find(nodeIdx) != s.end())
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
				status = CPXcutcallbackadd(_env, _cbdata, _wherefrom, nzcnt, -1, 'G', &idx[0], &val[0], CPX_USECUT_FORCE);
				if (status)
				{
					int warnCode = 999;
					std::string msg = "Failed to add integer cut.";
					warningMsg(NULL, __func__, msg.data(), warnCode);
				}
				else
				{
					cbHandlerPtr->totalNbCuts++;
					cbHandlerPtr->totalNbFracCuts++;
				}

				printf("");
			}
		}

		printf("");

		// Telling CPLEX that cuts have been created
		*_useraction_p = CPX_CALLBACK_SET;

#ifdef DEBUG
		// salva um arquivo .lp com o LP atual
		CPXLPptr nodelp = NULL;
		status = CPXgetcallbacknodelp(_env, _cbdata, _wherefrom, &nodelp);
		if (status)
		{
			fprintf(stderr, "Failed to get node lp pointer.\n");
			return status;
		}

		status = CPXwriteprob(_env, nodelp, ".\\fraccutset-dp.lp", "LP");
		if (status)
		{
			fprintf(stderr, "Failed to print lp file.\n");
			return status;
		}
#endif

	}

	return status;
}
