#ifndef PROBLEM_DATA_LOADER_H
#define PROBLEM_DATA_LOADER_H

#include "../include/Graph.h"

#include <string>

class ProblemDataLoader
{
public:
	ProblemDataLoader(std::string fileName, Graph* aGraph);
	~ProblemDataLoader(void);

	void load(void);

private:
	ProblemDataLoader(void);

	std::string file;
	Graph* g;
};


#endif
