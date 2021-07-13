# strategicvoting
Strategic Voting under Plurality and RCV
This repository contains the replication code for Eggers and Nowacki (n.d.).

# Scripts and running order

Due to the demandingness of computations, we ran the computations described in the paper on Stanford's [Sherlock](https://www.sherlock.stanford.edu/) computing cluster. We ran the scripts in the following order:

	* 'server.sh' runs the main algorithm for both Plurality and RCV, and saves the algorithm results for every case.
		* The batch/queueing command takes two arguments (e.g., `sbatch server.sh #1 #2'), where #1 is the parameter setting for the precision of beliefs (s), and #2 is the speed of learning in the algorithm (lambda).
	* 'all_gather.sh' takes the algorithm output and produces summary statistics.
	* 'random_iter.sh' runs the algorithm with random starting points (for RCV only) to demonstrate quasi-convergence irrespective of initial beliefs.
	* 'random_iter_process.sh' processes the algorithm output.
	* 'conjecture.sh' runs additional tests to show the mechanism conjectures.
	* 'server_four.sh' runs the algorithm adapted to the four-party setting (limited number of iterations)

# Figures reference

These will be updated shortly.