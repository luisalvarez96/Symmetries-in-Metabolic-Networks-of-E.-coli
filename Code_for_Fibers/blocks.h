/* Author: Ian Leifer <ianleifer93@gmail.com> */

#ifndef BUILDINGBLOCK_H
#define BUILDINGBLOCK_H
#include <vector>
#include <string>
using namespace std;

class BuildingBlock {
private:
	int id;
	vector<int> nodes;
	vector<int> colors;
public:
	BuildingBlock(int iD) {id = iD; colors.resize(1, 0);};
	bool addNode(int id, int color = -1);
	bool isAddableColor(int color);
	void print(string filename);
	int getNumberOfNodes() {return nodes.size();}
};

#endif
