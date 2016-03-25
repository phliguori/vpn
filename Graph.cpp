#include "../include/Graph.h"
#include "Util.h"

#include <typeinfo>
#include <sstream>

Graph::Graph():
	nVertices(0), nTerminals(0)
{
	return;
}

Graph::~Graph()
{
	return;
}

void Graph::toDot(std::ostream & _output) const
{
	std::string indentation("   ");
	_output << "graph " << "v" << nVertices << "_t" << nTerminals << "_a" << nArcs <<"\n";
	_output << "{\n";
	_output << indentation << "label = \"<some_label>\"\n";

	for (int i = 1; i <= nVertices; ++i)
	{
		Vertex v = vertices[i];
		_output << indentation << v.toDot();
	}

	for (int i = 1; i <= nVertices; ++i)
	{
		for (int j = 0; j < adjList[i].size(); ++j)
		{
			Arc a = adjList[i][j];

			//if (a.getHead().getCode() < a.getTail().getCode())
			if (a.isCanonicalForm())
				_output << indentation << a.toDot();
		}
	}

	_output << "}" << std::endl;

	return;
}

void Graph::reserve(int _nodes)
{
	vertices.clear();
	vertices.resize(_nodes, Vertex(0));

	for (int i = 0; i < _nodes; ++i)
	{
		vertices[i] = Vertex(i);
	}

	adjList.clear();
	adjList.resize(_nodes, std::vector<Arc>());

	arcSet.clear();
	arcSet.resize(_nodes, std::set<Arc>());

	return;
}

void Graph::addTerminal(Vertex _t)
{
	terminals.push_back(_t);

	return;
}

void Graph::addVertex(Vertex _v)
{
	if (_v.getCode() < vertices.size())
	{
		vertices[_v.getCode()] = _v;
	}
	else
	{
		std::stringstream warn;
		warn << "Vertex out of bound [" << _v.getCode() << ", " << vertices.size() << "]";
		warningMsg(typeid(*this).name(), __func__, warn.str().c_str(), 201);
	}

	return;
}

void Graph::addArc(Arc _arc)
{
	if (_arc.getHead() < adjList.size() && _arc.getTail() < adjList.size())
	{
		adjList[_arc.getHead().getCode()].push_back(_arc);
		adjList[_arc.getTail().getCode()].push_back(_arc.reverse());

		arcSet[_arc.getHead().getCode()].insert(_arc);
		arcSet[_arc.getTail().getCode()].insert(_arc.reverse());
	}
	else
	{
		std::stringstream warn;
		warn << "Arc out of bound [" << nVertices << ", (" << _arc.toString() << ")]";
		warningMsg(typeid(*this).name(), __func__, warn.str().c_str(), 202);
	}

	return;
}

