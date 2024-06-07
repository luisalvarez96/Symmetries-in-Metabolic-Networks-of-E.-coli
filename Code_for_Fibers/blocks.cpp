/* Author: Ian Leifer <ianleifer93@gmail.com> */

#include "blocks.h"
#include <iostream>
#include <fstream>

bool BuildingBlock::addNode(int id, int color) {
	for(int i = 0; i < nodes.size(); i++) {
		if(nodes[i] == id) {return 0;}
	}
	nodes.push_back(id);
	if(color != -1) {
		if(color >= colors.size()) {
			colors.resize(color + 1);
		}
		colors[color]++;
	}
	return 1;
}

bool BuildingBlock::isAddableColor(int color) {
	if(color >= colors.size()) {
		return 0;
	}
	if(colors[color] > 1) {
		return 1;
	} else {
		return 0;
	}
}

void BuildingBlock::print(string filename) {
	ofstream buildingBlockFile;
	buildingBlockFile.open(filename, ofstream::out | ofstream::app);
	/*for(int i = 0; i < groupoidIds.size(); i++) {
		fiberFile << i << "\t" << groupoidIds[i] << endl;
	}*/
	buildingBlockFile << id << ":\t";
	for(int i = 0; i < nodes.size() - 1; i++) {
		buildingBlockFile << nodes[i] << ", ";
	}
	buildingBlockFile << nodes[nodes.size() - 1] << endl;
	buildingBlockFile.close();
	/*cout << "Colors:" << endl;
	for(int i = 0; i < colors.size(); i++) {
		cout << i << ":\t" << colors[i] << endl;
	}*/
}
