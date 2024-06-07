/* Author: Ian Leifer <ianleifer93@gmail.com> */

#include "processor.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <algorithm>

Processor* Processor::p_Processor = 0;

Processor* Processor::getProcessor(int parallelId) {
	if(!p_Processor) {
		p_Processor = new Processor(parallelId);
	}
	return p_Processor;
}

Processor::Processor(int parallelId) {
    createFileNames(parallelId);
	readConfigFile();

	connections.resize(numberOfConnections);
	for(int i = 0; i < numberOfConnections; i++)
		connections[i].resize(weighted?3:2);

	for(int i = 0; i < numberOfNodes; i++)
		nodes.push_back(Node(i, 0));

	readConnectionsFile();
	findNoInputNodes();
}

void Processor::run() {
	/* Prepare colors for input */
	// different colors here just mean different type of nodes
	// defining all colors initially to be same
	// all nodes with no input will have different colors in directed case, because they are not synchronized
	// same goes for nodes with only connection to itself, but we want to allow this nodes to be part of groupoids,
	// cause they can synchronize
	vector<int> nodeColors(numberOfNodes, 0);
	int numberOfColors = noInputNodes.size() + onlyLoopbackInputNodes.size() + 1;
	for(int i = 0; i < noInputNodes.size(); i++) {
		nodeColors[noInputNodes[i]] = i;
	}

	for(int i = 0; i < onlyLoopbackInputNodes.size(); i++) {
		nodeColors[onlyLoopbackInputNodes[i]] = noInputNodes.size() + i;
	}

	numberOfColors = findGroupoids(numberOfNodes, connections, numberOfColors, nodeColors);

	printGroupoids(nodeColors);
	/*
	for(int i = 0; i < numberOfNodes; i++) {
		cout << i << "\t" << nodeColors[i] << endl;
	}*/
	
	/* Here we find building blocks */
	prepareColors(nodeColors, numberOfColors);

	/* First we create arrays of different colors */
	vector<vector<Node*>> colorSets(numberOfColors);

	for(int i = 0; i < numberOfNodes; i++) {
		int c = nodes[i].getColor();
		if(c != -1) {colorSets[c].push_back(&nodes[i]);}
	}

	/*for(int i = 0; i < numberOfColors; i++) {
		cout << "Color[" << i << "] = ";
		for(int j = 0; j < colorSets[i].size(); j++) {
			cout << colorSets[i][j]->getId() << "[" << colorSets[i][j]->getColor() << "]\t";
		}
		cout << endl;
	}*/

	for(int i = 0; i < numberOfColors; i++) {
		/* We use size of blocks as id, cause we want block ids to be from 0 to n, where n is number of building blocks */
		BuildingBlock bb(blocks.size());
		stack<Node*> *toAdd = new stack<Node*>;
		for(int j = 0; j < colorSets[i].size(); j++) {
			bb.addNode(colorSets[i][j]->getId());
			for(int k = 0; k < colorSets[i][j]->getNumberOfInputs(); k++) {
				if(bb.addNode(colorSets[i][j]->getInput(k)->getId(), colorSets[i][j]->getInput(k)->getColor())) {
					toAdd->push(colorSets[i][j]->getInput(k));
				}
			}
		}
		while(1) {
			bool added = false;
			stack<Node*> *newToAdd = new stack<Node*>;
			//bb.print();
			while(!toAdd->empty()) {
				Node* newNode = toAdd->top();
				toAdd->pop();
				//newNode->print();
				if(bb.isAddableColor(newNode->getColor())) {
					for(int k = 0; k < newNode->getNumberOfInputs(); k++) {
						if(bb.addNode(newNode->getInput(k)->getId(), newNode->getInput(k)->getColor())) {
							added = true;
							newToAdd->push(newNode->getInput(k));
							//cout << "Node " << newNode->getInput(k)->getId() << " with color " << newNode->getInput(k)->getColor() << " added" << endl;
							//bb.print();
						}
					}
				}
			}
			if(added == 0) {break;}
			toAdd = newToAdd;
		}
		delete toAdd;
		if(bb.getNumberOfNodes() != 0) {blocks.push_back(bb);}
	}
	
	for(int i = 0; i < blocks.size(); i++) {
		blocks[i].print(blocksFileName);
	}
}

void Processor::createFileNames(int parId) {
    if(parId != -1) {
        inputFileName = "parallel/adjacency" + to_string(parId) + ".txt";
        fiberFileName = "parallel/fibers" + to_string(parId) + ".txt";
        blocksFileName = "parallel/buildingBlocks" + to_string(parId) + ".txt";
    } else {
        inputFileName = "adjacency.txt";
        fiberFileName = "fibers.txt";
        blocksFileName = "buildingBlocks.txt";
    }
}

void Processor::prepareColors(vector<int> &nodeColors, int numberOfColors) {
	vector<int> colorDistribution(numberOfColors, 0);
	for(int i = 0; i < numberOfNodes; i++) {
		nodes[i].setColor(nodeColors[i]);
		colorDistribution[nodeColors[i]]++;
	}

	/* Now we put all elements of unique colors together to the color -1*/
	for(int i = 0; i < numberOfNodes; i++) {
		if(colorDistribution[nodes[i].getColor()] == 1) {nodes[i].setColor(-1);}
	}
}

void Processor::addConnection(int source, int destination) {
	if(source < 0 || source >= numberOfNodes || destination < 0 || destination >= numberOfNodes) {
		cout << "Error: Trying to add connection out of bound. Number of nodes = " << numberOfNodes << ", source = " << source << ", destination = " << destination << endl;
	}
	nodes[source]     .addOutput(&nodes[destination]);
	nodes[destination].addInput (&nodes[source]);
}

void Processor::readConfigFile() {
// First line is number of nodes
// Second is if graph is directed(true or false)
// Third is weighted(true or false)
// Forth is number of different weights
/* Here the small remark included. We assume that we know exactly the amount of different possible
n weights and they are from 0..n. From the perspective of the algorithm there is no difference if
the weights are numbers or names, but the person who creates input has to take care of it being exactly n in form 0..n. */
// then connections follow up
	string line;
	ifstream config;

	config.open(inputFileName, ifstream::in);
	cout << inputFileName << endl;

	//cout << "estoy aqui1" << endl;
	getline(config, line, '\n');
	numberOfNodes = stoi(line);

	//cout << "estoy aqui2" << endl;
	getline(config, line, '\n');
	directed = stoi(line);

	//cout << "estoy aqui3" << endl;
	getline(config, line, '\n');
	weighted = stoi(line);

	//cout << "estoy aqui4" << endl;
	getline(config, line, '\n');
	if(weighted) {
		numberOfWeights = stoi(line);
	} else {
		numberOfWeights = 0;
	}

	numberOfConnections = 0;
	while(1) {
		if(!std::getline(config, line, '\n')) {break;}
		numberOfConnections++;
	}
	config.close();
}

void Processor::readConnectionsFile() {
	string line;
	ifstream config;

	config.open(inputFileName, ifstream::in);

	// skip 3 lines with config data
	getline(config, line, '\n');
	getline(config, line, '\n');
	getline(config, line, '\n');
	getline(config, line, '\n');

	int i = 0;
	while(1) {
		if(!getline(config, line, '\t')) {break;}
		connections[i][0] = stoi(line);
		getline(config, line, weighted?'\t':'\n');
		connections[i][1] = stoi(line);
		if(weighted) {
			getline(config, line, '\n');
			connections[i][2] = stoi(line);
		}
		addConnection(connections[i][0], connections[i][1]);
		i++;
	}
	config.close();
	// print connections
	/*for(int i = 0; i < numberOfConnections; i++) {
		cout << "Connection " << i << ": " << connections[i][0] + 1 << " -> " << connections[i][1] + 1;
		if(weighted) {cout << ". Weight = " << connections[i][2];}
		cout << endl;
	}*/
}

void Processor::printGroupoids(vector<int> groupoidIds) {
	ofstream fiberFile;
	fiberFile.open(fiberFileName);
	for(int i = 0; i < groupoidIds.size(); i++) {
		fiberFile << i << "\t" << groupoidIds[i] << endl;
	}
	fiberFile.close();
}

void Processor::findNoInputNodes() {
	if(numberOfNodes == 0 || directed == 0) {return;}
	vector<bool> hasInput(numberOfNodes, 0);
	for(int i = 0; i < numberOfConnections; i++) {
		if(connections[i][0] == connections[i][1]) {continue;}
		hasInput[connections[i][1]] = 1;
	}
	for(int i = 0; i < numberOfConnections; i++) {
		if(hasInput[connections[i][0]] == 0 && connections[i][0] == connections[i][1]) {
			onlyLoopbackInputNodes.push_back(connections[i][0]);
			hasInput[connections[i][0]] = 1;
		}
	}
	for(int i = 0; i < numberOfNodes; i++) {
		if(hasInput[i] == 0) {noInputNodes.push_back(i);}
	}
}

void Processor::calculateVectors(vector<int> nodeColors, vector< vector<int> > &vectors) {
	for(int i = 0; i < connections.size(); i++) {
		if(directed == false) {
			int pos = 0;
			if(numberOfWeights == 0) {
				pos = nodeColors[connections[i][1]];
			} else {
				pos = nodeColors[connections[i][1]] * numberOfWeights + connections[i][2];
			}
			vectors[connections[i][0]][pos]++;
		}
		int pos = 0;
		if(numberOfWeights == 0) {
			pos = nodeColors[connections[i][0]];
		} else {
			pos = nodeColors[connections[i][0]] * numberOfWeights + connections[i][2];
		}
		vectors[connections[i][1]][pos]++;
	}

	// output vector value
	/*for(int i = 0; i < vectors.size(); i++) {
		for(int j = 0; j < vectors[0].size(); j++) {
			cout << "vectors[" << i << "][" << j << "] = " << vectors[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;*/
}

// returns new number of colors
int Processor::classifyNodes(vector< vector<int> > vectors, vector<int> &nodeColors) {
	// first let`s find how many unique types of vectors are out there
	vector< vector<int> > vectorTypes;
	vectorTypes.push_back(vectors[0]);

	for(int i = 1; i < nodeColors.size(); i++) {
		bool add = 1;
		for(int j = 0; j < vectorTypes.size(); j++) {
			if(vectors[i] == vectorTypes[j]) {add = 0;}
		}
		if(add == 1) {vectorTypes.push_back(vectors[i]);}
	}
/*
	// output unique vector types
	for(int i = 0; i < vectorTypes.size(); i++) {
		for(int j = 0; j < vectorTypes[0].size(); j++) {
			cout << "vectorTypes[" << i << "][" << j << "] = " << vectorTypes[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
*/
	// now let`s reshuffle node types for next step
	/* this for loop was i = 1 ... i++ before. I decided that it was starting from 1 because I copied this loop
	from one above, I don`t see any reason for 1 here, so I put 0, cause it causes problem with new nodes without
	input, cause they need color 0.
	If the reason was something else, don`t forget to put it here in commentary */
	for(int i = 0; i < nodeColors.size(); i++) {
		for(int j = 0; j < vectorTypes.size(); j++) {
			if(vectors[i] == vectorTypes[j]) {nodeColors[i] = j + noInputNodes.size();}
		}
	}

	if(directed != 0) {
		for(int i = 0; i < noInputNodes.size(); i++) {
			nodeColors[noInputNodes[i]] = i;
		}
	}
	return vectorTypes.size() + noInputNodes.size();
}

int Processor::findGroupoids(int numberOfNodes, vector< vector<int> > groupoidConnections, int numberOfColors, vector<int> &nodeColors) {
	while(1) {
		// create 2D vector array to store all vectors belonging to each node
		/* Explanation why array is of size numberOfNodes x (numberOfColors * numberOfWeights). There are two ways how to do it.
		Either it can be done as a 3D array and then we will need two realisations for weighted and non-weighted design.
		Or vectors themselves can be formed in a bit weird way, but we will classify nodes comparing vectors not worrying about their structure.
		It improves readability and simpliness only paying with the strange enumeration of array */
		vector< vector<int> > vectors(numberOfNodes);
		for(int i = 0; i < numberOfNodes; i++) {
			vectors[i].resize(numberOfColors * (weighted?numberOfWeights:1));
		}

		calculateVectors(nodeColors, vectors);
		int nOC = classifyNodes(vectors, nodeColors);
		if(nOC == numberOfColors) {break;}
		else {numberOfColors = nOC;}
	}
	return numberOfColors;
}
