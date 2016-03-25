#include "../include/ProblemDataLoader.h"

#include "Util.h"
#include <sstream>
#include <fstream>
#include <iostream>

ProblemDataLoader::ProblemDataLoader(std::string fileName, Graph * aGraph)
{
	g = aGraph;
	file = fileName;

	return;
}

ProblemDataLoader::~ProblemDataLoader(void)
{
	return;
}

void ProblemDataLoader::load(void)
{
	int nNodes, nTerminals, nArcs;

	std::ifstream fin(file, std::ios::in);

	if (!fin.good())
	{
		std::string msg = "Cannot open " + file;
		errorMsg(typeid(*this).name(), __func__, msg.c_str(), 1);
		
		exit(EXIT_FAILURE);
	}

#ifdef VERBOSE
		std::cout << "Reading file: " << file << "\n";;
#endif

	//reading general instance information
	fin >> nNodes >> nTerminals >> nArcs;

#ifdef VERBOSE
	std::cout << "\nNodes(" << nNodes << "), Terminals("<< nTerminals << "), Arcs(" << nArcs << ")\n";
#endif

	g->nVertices = nNodes;
	g->nTerminals = nTerminals;
	g->nArcs = nArcs;

	g->reserve(nNodes + 1);

	//reading terminal nodes information
	for (int i = 0; i < nTerminals; ++i)
	{
		int code;
		double ingree, egree;

		fin >> code >> ingree >> egree;

#ifdef VERBOSE
		std::cout << "\nterminal " << code << ", ingree " << ingree << ", egree " << egree;
#endif

		Vertex terminal(code, ingree, egree, true);
		g->addVertex(terminal);
		g->addTerminal(terminal);
	}

	//reading graph edges (network links)
	for (int i = 0; i < nArcs; ++i)
	{
		int head, tail;
		double capacity;

		fin >> head >> tail >> capacity;

		Arc aux(g->vertices[head], g->vertices[tail], capacity);

		if (head == tail)
		{
			int warnCode = 010;
			std::string msg = "Self arc : (" + aux.toString() + ")";
			warningMsg(typeid(*this).name(), __func__, msg.data(), warnCode);
			continue;
		}

		g->addArc(aux);
#ifdef VERBOSE
		std::cout << "\nArc(" << aux.toString() << ")";
#endif
	}

	return;
}
