# electronicsdesignproject
Summer Term Electronics Design Project - Circuit simulator
Team Members : Marco, Weihan , James

Start Date: 11/5/2021
End Date: 13/6/2021

Report Link: (Will be uploaded)
Video Link: (Will be uploaded)

Abstract: 

The objective of this project was to implement a circuit simulator that performs AC analysis on a given circuit, 
returning the magnitude and phase of a transfer function over a range of stepped frequencies to form bode plots. 
Starting from a simplified netlist of sources, linear and non-linear devices, a conductance matrix system was 
used to represent the currents, voltages, and impedances in the network. From this, voltages of all nodes over 
a range of frequencies were obtained to derive the gain and phase of the circuit. Several tests were then conducted 
to determine the accuracy and efficiency of the programme. 

When using the software, the user will input a netlist of components as well as the simulation command. The input file
is then parsed to store the components and parameters in a data structure. From the data stored, the relevant nodal 
analysis equations formed by applying KCL can be represented in a conductance matrix system. The conductance matrix,
contains the conductances of every component. The position of a particular conductance in the matrix depends on the 
nodal connections of the impedance device.

The conductance matrix system in the form of G • V = I is analogous to the linear algebra system A • X = B. Thus, 
the linear system can be manipulated to give V by applying LU decomposition then solving, allowing the voltages at 
each node to be derived.

As reactive components such as inductors have different conductances for different frequencies, the conductance matrix
system needs to be solved for every given frequency defined by the simulation command during AC analysis, giving us the
necessary output data to form the bode plots.

Details:

Over the past month we’ve built a circuit simulator that can perform AC analysis on circuits containing linear and non-linear devices. 
Our aim was to design a circuit simulator that conducts an AC simulation with an accuracy of at least 95% compared to LTSpice within 
a reasonable duration.

To expedite the development process, our software was broken down into blocks that fulfil the main functional requirements. 

First the parsing block is used to parse a reduced SPICE format netlist containing the circuit information, whilst storing the 
components in a suitable data structure and then initialising the small-signal conversion of the circuit.

Secondly, the conductance matrix construction block uses the data stored to build conductance matrices.The Multisource Analysis block 
is then used to solve them to find nodal voltages.
If non-linear devices are present, the DC iteration block applies Newton-Raphson method to derive the DC operating point and find 
small signal parameters.

Finally, the AC analysis block should output the magnitude and phase of the transfer function at each 
frequency in a CSV file format to produce the bode plots. 
 
The multisource analysis block is particularly versatile because the algorithm implemented in it is also applied in the DC Iteration
and AC Analysis blocks to process equivalent circuits that contain multiple sources, making it a pivotal aspect of our simulator. 

Taking a deeper look into the algorithm implemented, the multisource analysis block can be broken down into 4 separate source analysis 
sub-blocks. The order in which these sub-blocks are arranged is key towards the optimization of the algorithm. Current sources and VCCSes 
are considered first because they add more variables to the conductance matrix, increasing its complexity. Whilst voltage sources are 
considered afterwards because they overwrite the conductance rows of the matrix system, adding more zeros in the process. The addition 
of zeroes simplifies the matrix system by reducing the number of operations required to perform LU decomposition, shortening the time 
taken to solve the system for nodal voltages.




Files Details:

Main.cpp - This contains the circuit simulator software. To test it, first write a .txt file containing a reduced netlist format.
Then, run the program, input the file name (i.e. netlist_bjt.txt), node number (i.e. 3), and source (i.e. V1).
The execution time will be printed. A file named output.txt will be created which contains the stepped frequency, 
Magnitude, and phase in CSV Format. The file can be plotted using matlab with the file plotsim.m

quicktestallcases.cpp - automated time test programme that first creates a specified number of netlist with increasing node number. 
After the creation of a file, it will read it and simulate the netlist that it contains and output the execution time. The purpose
of this is to find out the relationship between execution time against the number of nodes contained in a circuit. 

netlist_xxx - netlist files that contains the reduced spice format of a circuit.














Github Naming Convention:

Stable - stable	- Accepts merges from Working and Hotfixes
Working - master - Accepts merges from Features/Issues and Hotfixes

Features/Bug topic-* - Always branch off HEAD of Working
Hotfix	hotfix-*	- Always branch off Stable

If you are working on a new feature i.e. Solving inverse of conductance matrix, start a new branch named feature/solveconductancematrix 

Side note:
1. Head is a pointer that points to the latest commit of a branch. Always create a new branch of the HEAD of master when trying to add a new feature

2. DO NOT use hotfix for non-urgent bugs. Hotfix is only used to fix bugs off the Stable version of program. Use branch (bug)

If unsure of the convention, always refer to the link below!
https://gist.github.com/digitaljhelms/4287848
(There are a lot of naming conventions out there. This is the one we are using for this project)



Commit Messages convention:

<type>[optional scope]:<description>
etc. fix: solve inverse matrix
which means the feature of solve inverse matrix has a bug and is fixed

info - used for info and comments
fix - used for bug fixes
feat - used for feature
hotfix - used for bug fixes off stable


some useful Git Commands:

git fetch --prune (Delete merged branches in vscode, restart vscode after command)

git status (List which files are staged, unstaged, and untracked.)

git log (Shows commits history)

Git commands cheatsheet:
https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet


