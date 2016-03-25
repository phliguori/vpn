#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Util.h"
#include "Row.h"
#include "../include/Graph.h"

#include <set>
#include <map>
#include <string>
#include <locale>

class Constraint : public Row
{
public:
	typedef enum
	{
		C_ASSIGNMENT = 100,
		C_FACILITY_ACTIVATION = 200,
		C_TREE = 300,
		C_CONNECTEDNESS = 400,
		C_EDGE_ACTIVATION = 500,
		C_FLOW_CONSERVATION = 600,
		C_CAPACITY = 700,
		C_DUAL_FEASIBILITY = 800,
		C_UNKNOWN = 999
	} CONTYPE;

	class Comparator
	{
	public:
		bool operator()(const Constraint & _left, const Constraint & _right) const
		{
			if (_left.type != _right.type)
				return _left.type < _right.type;
			
			if (_left.vertex != _right.vertex)
				return _left.vertex < _right.vertex;

			if (_left.s != _right.s)
				return _left.s < _right.s;

			if (_left.t != _right.t)
				return _left.t < _right.t;

			if (_left.arc != _right.arc)
				return _left.arc < _right.arc;

			if (_left.vertexSet.size() != _right.vertexSet.size())
				return _left.vertexSet.size() < _right.vertexSet.size();

			long cleft = hashFunc(_left.toString());
			long cright = hashFunc(_right.toString());

			return cleft < cright;
		}
	};

	// Constructors & destructors
	Constraint(void);
	Constraint(int _maxnnz, ROWSENSE _sense, double _rhs);
	Constraint(const Constraint & _constraint);
	virtual ~Constraint(void);

	// Getters
	int getClass() const { return consClass; }
	int getType(void) const { return type; }
	Vertex getS(void) const { return s; }
	Vertex getT(void) const { return t; }
	Vertex getNode(void) const { return vertex; }
	Arc getArc(void) const { return arc; }

	// Setters
	void setClass(int _class) { consClass = _class; }
	void setType(CONTYPE _type) { type = _type; }
	void setNode(Vertex _vertex) { vertex = _vertex; }
	void setS(Vertex _vertex) { s = _vertex; }
	void setT(Vertex _vertex) { t = _vertex; }
	void setArc(Arc _arc) { arc = _arc; }

	// Util
	std::string toString(void) const;
	void addToSet(int vertex_id) { vertexSet.insert(vertex_id); }

private:
	static int consCounter;

	int id;

	int consClass;
	CONTYPE type;
	Arc arc;
	Vertex vertex;
	Vertex s, t;

	std::set<int> vertexSet;

	long int hashValue;
};

long hashFunc(const Constraint & _constraint);

typedef std::map<Constraint, int, Constraint::Comparator> ConstraintHash;

#endif