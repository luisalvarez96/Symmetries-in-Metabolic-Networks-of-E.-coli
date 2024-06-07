/* Author: Ian Leifer <ianleifer93@gmail.com> */

#ifndef PROCESSOR_H
#define PROCESSOR_H
#include <vector>
#include <string>
#include "node.h"
#include "blocks.h"
using namespace std;

//#define inputFile "adjacency.txt"

class Processor {
private:
	Processor(int parallelId);
	static Processor *p_Processor;
public:
	static Processor* getProcessor(int parallelId);
private:
	int numberOfNodes;
	vector<Node> nodes;
	int numberOfConnections;
	void addConnection(int source, int destination);
	bool directed;
	bool weighted;
	int numberOfWeights;
	void readConfigFile();
    string inputFileName;
    string fiberFileName;
    string blocksFileName;
    void createFileNames(int parId);
private:
	// create 2D vector array to store all connections
	vector< vector<int> > connections;
	void readConnectionsFile();
	void printGroupoids(vector<int> groupoidIds);
private:
	vector<int> noInputNodes;
	vector<int> onlyLoopbackInputNodes;
	void findNoInputNodes();
	void calculateVectors(vector<int> nodeColors, vector< vector<int> > &vectors);
	int classifyNodes(vector< vector<int> > vectors, vector<int> &nodeColors);
	int findGroupoids(int numberOfNodes, vector< vector<int> > groupoidConnections, int numberOfColors, vector<int> &nodeColors);
private:
	vector<BuildingBlock> blocks;
	void prepareColors(vector<int> &nodeColors, int numberOfColors);
public:
	void run();
};

#endif
