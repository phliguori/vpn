#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <sys/stat.h>

#include <map>
#include <random>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "Util.h"
#include "Separator.h"
#include "SolutionWriter.h"
#include "../include/Graph.h"
#include "../include/ProblemDataLoader.h"
#include "../include/Flow.h"
#include "../include/FacilityLocation.h"
#include "../include/CutSetSeparation.h"

#ifndef PATH_SEPARATOR
#ifdef WIN32
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif
#endif

// * Command line syntax
// * -i <input_file>: indicates that input_file should be read 
// * 
// * -o <output_file>: indicates the name of file in which solution should be written
// * 
// * -d <directory>: specifies the directory where the files should de read/written
// * 
// * -s <option>: used to specify the separation routine, whether (f)ractional separation
// * should be carry over, or (i)nteger separation only, or (b)oth. Let the parameter
// * empty to not separate constraintes
// * 
// * -e: used to indicate wheter the algorithm of extended relaxation should
// * be used to calculate de lp relax
// *
// * -t <time_limit>: to set the time limit value
// *
// * -n <node_limit>: to set the maximum number of nodes to be explored in the branch-and-bound
// *
int main(int _argc, char* _argv[])
{
	int nodesLim = -1;
	double timeLim = -1.;
	bool solveExtendedLR = false;
	SolutionWriter* writer;
	Separator* separator=NULL;
	std::string inputFile;
	std::string fullInputFile;
	std::string outputFile;
	std::ostream & outputStream = std::cout;
	std::ofstream ofile;

	if (cmdOptionExists(_argv, _argv + _argc, "-i"))
	{
		inputFile += getCmdOption(_argv, _argv + _argc, "-i");
		fullInputFile += (cmdOptionExists(_argv, _argv + _argc, "-d")) ? getCmdOption(_argv, _argv + _argc, "-d") : ".";
		fullInputFile += PATH_SEPARATOR;
		fullInputFile += inputFile;
	}
	else
	{
		helpUsage();
		errorMsg(NULL, __func__, "Not valid parameters", -1);
	}

	if (cmdOptionExists(_argv, _argv + _argc, "-o"))
	{
		outputFile += (cmdOptionExists(_argv, _argv + _argc, "-d")) ? getCmdOption(_argv, _argv + _argc, "-d") : ".";
		outputFile += PATH_SEPARATOR;
		outputFile += getCmdOption(_argv, _argv + _argc, "-o");
		outputFile += ".sol";

		ofile.open(outputFile.c_str(), std::ios::app);
		if (!ofile.good())
		{
			//TODO: print header
			printf("");
		}
		else
			printf("");
	}

	std::string tmp;
	for (int i = 0; i < _argc; ++i)
	{
		char aux[512];
		strcpy_s(aux, _argv[i]);

		tmp += aux;
		tmp += " ";
	}

	if (cmdOptionExists(_argv, _argv + _argc, "-s"))
	{
		int f = 'f';
		int i = 'i';
		int b = 'b';
		int separationChoice = *(getCmdOption(_argv, _argv + _argc, "-s"));
		switch (separationChoice)
		{
		case 102:
			separator = new Separator(NULL, separateFractionalCutSet);
			break;
		case 105:
			separator = new Separator(separateIntegerCutSet, NULL);
			break;
		case 98:
			separator = new Separator(separateIntegerCutSet, separateFractionalCutSet);
			break;
		default:
			separator = new Separator(NULL, NULL);
			break;
		}
	}

	if (cmdOptionExists(_argv, _argv + _argc, "-e"))
		solveExtendedLR = true;

	if (cmdOptionExists(_argv, _argv + _argc, "-n"))
	{
		nodesLim = atoi(getCmdOption(_argv, _argv + _argc, "-n"));
	}

	if (cmdOptionExists(_argv, _argv + _argc, "-t"))
	{
		timeLim = atof(getCmdOption(_argv, _argv + _argc, "-t"));
	}
	

	std::cout << "---- Instance name: "<< fullInputFile << std::endl;


	Graph* g = new Graph();
	
	ProblemDataLoader* loader = new ProblemDataLoader(fullInputFile, g);
	loader->load();
	delete loader;

	SolverMIP::SOLVERSTAT status = SolverMIP::SOLVERSTAT::SOLVERSTAT_UNKNOWN;

	// * --------------------------------------------------
	// * SOLVING THE LINEAR RELAXATION
	// * --------------------------------------------------
	FacilityLocation* lpFL = new FacilityLocation(g);

	Flow* mipflow = new Flow(g);
	mipflow->loadProblem("flow", SolverMIP::FOSENSE::MINIMIZATION, SolverMIP::PROBTYPE::INTEGER);
	status = mipflow->solve(SolverMIP::METHOD::METHOD_MIP, true);
	printf("");
	double flow = 0.;
	if (status == SolverMIP::SOLVERSTAT::SOLVERSTAT_MIPOPTIMAL || status == SolverMIP::SOLVERSTAT::SOLVERSTAT_LPOPTIMAL ||
		status == SolverMIP::SOLVERSTAT::SOLVERSTAT_FEASIBLE)
	{
		flow = mipflow->getObjVal();
		mipflow->getX();
		mipflow->builSolutionGraph();

		std::cout << std::endl << "------";
		std::cout << std::endl;
		std::cout << mipflow->printXSol();

		std::ofstream dotFile(fullInputFile + ".dot");
		mipflow->solutionToDot(dotFile);
		dotFile.close();
	}

	// Instantiation solution writer
	writer = new SolutionWriter(lpFL);

	lpFL->loadProblem("facilityLocation", SolverMIP::FOSENSE::MINIMIZATION, SolverMIP::PROBTYPE::CONTINUOUS);
	if (solveExtendedLR)
		status = lpFL->solveLRExtentedFormulations(fractionalCutSet);
	else
		status = lpFL->solve(SolverMIP::METHOD::METHOD_PRIMAL);

	double lpRelax = -1.;
	if (status == SolverMIP::SOLVERSTAT::SOLVERSTAT_MIPOPTIMAL || status == SolverMIP::SOLVERSTAT::SOLVERSTAT_LPOPTIMAL ||
		status == SolverMIP::SOLVERSTAT::SOLVERSTAT_FEASIBLE)
	{
		lpFL->getX();
		lpRelax = lpFL->getObjVal();

		std::cout << std::endl << "------";
		std::cout << std::endl;
		std::cout << lpFL->printXSol();
		std::cout << writer->printStatus();
	}

	// Writing output solution
	if (ofile.good())
	{
		ofile << inputFile << ";";
		ofile << g->nVertices << ";";
		ofile << g->nTerminals << ";";
		ofile << g->nArcs << ";";
		ofile << writer->toString();
		ofile.close();
	}

	delete writer;
	delete lpFL;

	return 0;
	// * --------------------------------------------------
	// * SOLVING THE INTEGER PROBLEM
	// * --------------------------------------------------
	FacilityLocation* mipFL = new FacilityLocation(g);
	mipFL->loadProblem("facilityLocation", SolverMIP::FOSENSE::MINIMIZATION, SolverMIP::PROBTYPE::INTEGER);

	// Instantiation solution writer
	writer = new SolutionWriter(mipFL);

	CBData* cbData = new CBData(g, &(mipFL->vHash), mipFL->getNCols());
	//Separator* cutSet = new Separator(NULL, NULL);
	//Separator* cutSet = new Separator(separateIntegerCutSet, NULL);
	//Separator* cutSet = new Separator(NULL, separateFractionalCutSet);
	//Separator* cutSet = new Separator(separateIntegerCutSet, separateFractionalCutSet);
	
	mipFL->setNodesLim(nodesLim);
	mipFL->setTimeLim(timeLim);
	mipFL->setExternalLpRelax(lpRelax);
	mipFL->setSeparationRoutine(cbData, separator);
	status = mipFL->solve(SolverMIP::METHOD::METHOD_MIP);

	double mip = 0.;
	if (status == SolverMIP::SOLVERSTAT::SOLVERSTAT_MIPOPTIMAL || status == SolverMIP::SOLVERSTAT::SOLVERSTAT_LPOPTIMAL ||
		status == SolverMIP::SOLVERSTAT::SOLVERSTAT_FEASIBLE)
	{
		mipFL->getX();
		mipFL->builSolutionGraph();

		std::cout << std::endl << "------";
		std::cout << std::endl;
		std::cout << mipFL->printXSol();

		std::ofstream dotFile(fullInputFile + ".dot");
		mipFL->solutionToDot(dotFile);
		dotFile.close();
	}

	std::cout << writer->printStatus();
	std::cout << std::endl << "CallBack Data - IntCuts: " << cbData->totalNbIntCuts;
	std::cout << std::endl << "CallBack Data - FracCuts: " << cbData->totalNbFracCuts;
	std::cout << std::endl << "CallBack Data - TotalCuts: " << cbData->totalNbCuts;
	std::cout << std::endl;

	// Writing output solution
	if (ofile.good())
	{
		ofile << inputFile << ";";
		ofile << g->nVertices << ";";
		ofile << g->nTerminals << ";";
		ofile << g->nArcs << ";";
		ofile << mipFL->getActiveRouters() << ";";
		ofile << mipFL->getActiveArcs() << ";";
		ofile << writer->toString();
		ofile.close();
	}

	delete writer;
	delete separator;
	delete cbData;

	delete g;
	delete mipFL;

	std::cout << std::endl;
	warningMsg(NULL, __func__, "End of execution", 0);
	
	return 0;
}