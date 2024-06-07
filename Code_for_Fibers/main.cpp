/* Author: Ian Leifer <ianleifer93@gmail.com> */

#include "processor.h"
#include <cstdlib>

int main(int argc, char *argv[]) {
  Processor *processor;
  int processId = -1;
  if(argc == 2) {
    processId = atoi(argv[1]);
  }
  processor = processor->getProcessor(processId);
	processor->run();
	return 0;
}
