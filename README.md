# Codes for paper *Symmetries in metabolic networks of E. coli*

The flow for analyzing each of the 4 metabolic networks is the following:
1. Get fibers for the network using the script **fiber.R** from the directory **Codes_for_fibers**.
2. Run *newBlocks(...)* from **changing_blocks_class.R** to run complete analysis of the network Building Blocks (BBs) obtained from previous step.
3. *plot_collapsed(...)* can be used to plot individual BBs if wanted.
4. **similarity.R** is used to analyze the similarity of the GO terms of different clusters of enzymes.

## 1. Codes_for_fibers

This is a modification of an older version of the codes that were made available on [this](https://github.com/makselab/fibrationSymmetries) repository. Most of the changes made were so that it would be easier to work with the rest of the codes for the paper. (Its important to note that all the scripts within this directory need to be on the same directory as **fiber.R**)

From **fiber.R** run the function *main(...)* to obtain the coloring of the nodes, or fibers, of the network (by internally calling all of the scripts in C++).
This calls the functions on the **functions.R** script, firstly for preparing the network files to be run on the C++ codes and after the C++ codes are run, to retrieve the obtained colors. After the colors are obtained the script calls **classifier.R** to classify the BBs. 

*main(...)* outputs a list with 2 entries. The first one, called *Nodes*, corresponds to a three column data frame. With the columns *Id*, *Label*, *FiberId*, corresponding to the internal id number of the node, its name or label, and the fiber id of the fiber that the node belongs to. The second entry, called *Blocks*, corresponds to a data frame listing the BBs. This list returns the same number of BBs and fibers, as it constructs a BB for each individual fiber.

## 2. *newBlocks(...)* from **changing_blocks_class.R**

This function fixes a couple of issues with the BBs output from *main(...)* when dealing with *Fibonacci* BBs. Firstly, it includes loops between the fibers and regulators that cross other nodes that do not directly regulate the fibers. Secondly, it removes repetetion of *Composite Feedback Fibonacci* BBs. This type of BBs occur when there is a feedback loop between multiple fibers, so by constructing the BB for each individual fiber you are essentially reconstructing the same BB, which leads to repetition in the intial list of BBs output from *main(...)*. For ease of (re-)classifying and anaylzing, it adds a column (*FiberIds*) listing all the FiberIds for all the fibers included in the BB and separate the fibers that are working as regulators in the BB. Column *FiberRegIds* lists the fiber ids of the regulating fibers in the BB, and *FiberRegs* lists the name of this nodes. *Fibers* includes the name of all nodes in fibers regulated by the BB and *Regulators* includes the name of external regulators (nodes that don't belong to a fiber). Lastly, column *nl* corresponds to the (re-)classification of the BB, column *r* to the branching ratio of the BB and *NumbFibs* how many fibers are regulated in the BB.

The output of *newBlocks(...)* corresponds to a list of lenght 3. First entry *Blocks* is the list of re-classified and completed BBs, second entry *counts* counts the number of fibers that appear in more than one BB, their frequency and the classes of each of these BBs, last entry *Fibers* lists all fibers and their classification.
