#ifndef GRAPH_H
#define GRAPH_H

#include "Util.h"

#include <set>
#include <string>
#include <vector>
#include <iostream>

#include "Vertex.h"
#include "Edge.h"

class Graph
{
public:
	Graph();
	~Graph();

	void toDot(std::ostream & _output = std::cout) const;

	void reserve(int nodes);
	
	void addTerminal(Vertex _t);
	void addVertex(Vertex _v);
	void addArc(Arc _arc);

	int nVertices;
	int nTerminals;
	int nArcs;

	std::vector<Vertex> terminals;
	std::vector<Vertex> vertices;
	std::vector< std::vector<Arc> > adjList;
	std::vector< std::set<Arc>> arcSet;
};

#endif
