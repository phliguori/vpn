#ifndef VARIABLE_H
#define VARIABLE_H

#include "Column.h"
#include "../include/Graph.h"

#include <map>

class Variable:public Column
{
public:
	typedef enum type
	{
		V_X = 10,
		V_Y = 20,
		V_Z = 30,
		V_W = 40,
		V_UNKNOWN = 50
	} VARTYPE;

	class Comparator
	{
	public:
		bool operator()(const Variable & _left, const Variable & _right) const
		{
			return (_left < _right);
		}
	};

	class PtrComparator
	{
	public:
		bool operator()(const Variable* _left, const Variable* _right) const
		{
			Variable::Comparator compare;
			return compare(*_left, *_right);
		}
	};

	bool operator< (const Variable & _right) const 
	{
		if (this->type != _right.type)
			return (this->type < _right.type);

		if (this->category != _right.category)
			return (this->category < _right.category);

		if (this->v1 != _right.v1)
			return (this->v1 < _right.v1);

		if (this->v2 != _right.v2)
			return (this->v2 < _right.v2);

		return (this->arc < _right.arc);
	}

	bool operator!= (const Variable & _right) const
	{
		return (*this < _right) || (_right < *this);
	}

	bool operator== (const Variable & _right) const
	{
		return !(*this != _right);
	}

	bool operator> (const Variable & _right) const
	{
		if (*this != _right)
			return _right < *this;
		
		return false;
	}

	Variable(void);
	Variable(COLTYPE _type, double _objcoef, double _lb, double _ub);
	virtual ~Variable(void);

	std::string toString (void) const;

	// Getters
	VARTYPE getType(void) const { return type; }
	Arc getArc(void) const { return arc; }
	Vertex getVertex1(void) const { return v1; }
	Vertex getVertex2(void) const { return v2; }
	char getCategory(void) const { return category; }

	// Setters
	void setType(VARTYPE _type) { type = _type; }
	void setArc(Arc _arc) { arc = _arc; }
	void setVertex1(Vertex _vertex) { v1 = _vertex; }
	void setVertex2(Vertex _vertex) { v2 = _vertex; }
	void setCategory(char _category) { category = _category; }

protected:
	static int varCounter;

	int id;

	char category;
	VARTYPE type;
	Arc arc;
	Vertex v1, v2;
};

typedef std::map<Variable, int, Variable::Comparator> VariableHash;
typedef std::map<Variable*, int, Variable::PtrComparator> PtrVariableHash;
typedef std::map<Variable::VARTYPE, VariableHash> TypeVariableHash;
typedef std::map<Variable::VARTYPE, VariableHash>* TypeVariableHashPtr;

#endif