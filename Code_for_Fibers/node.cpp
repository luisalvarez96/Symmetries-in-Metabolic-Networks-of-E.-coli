/* Author: Ian Leifer <ianleifer93@gmail.com> */

#include "node.h"
#include <iostream>

//TODO: check for duplicated connections
void Node::addInput(Node* entry) {
	input.push_back(entry);
}

void Node::addOutput(Node* entry) {
	output.push_back(entry);
}

void Node::print() {
	cout << "Printing node " << id << ". Color = " << color << endl;
	cout << "Number of connections (i/o) = " << input.size() << "/" << output.size() << endl;
}
