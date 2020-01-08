# strategicvoting
Strategic Voting under Plurality and RCV

# Scripts and running order

Server-side script. Due to the demandingness of computations, we recommend that the following script be executed on a computing cluster or another machine with adequate memory (we ran it with the following configuration: ....)

when running the script from the console, the command takes two arguments appended at the end that define the parameter values:

RScript replication_server.R A B

where 'A' decides the precision of beliefs (s) and 'B' decides the speed of learning in the iterative polling algorithm (lambda). We ran the script for the following values:

1 1 (s = 10, lambda = 0.05)
1 2 (s = 10, lambda = 0.01)
1 3 (s = 10, lambda = 0.10)
4 1 (s = 55, lambda = 0.05)
4 2 (s = 55, lambda = 0.01)
4 3 (s = 55, lambda = 0.10)
6 1 (s = 85, lambda = 0.05)
6 2 (s = 85, lambda = 0.01)
6 3 (s = 85, lambda = 0.10)

The scripts generate data (saved in output/files) and figures (saved in output/figs_v2).

Subsequently, we run the following scripts to generate the remaining figures:
* gather_v_vecs_v2.R
* random_starting_iterations_v2.R
* summary_new.R

# Figures reference

Figures in .pdf correspond to:
Figure 1 -- output/figs_v2/1/85/v_vec_path_v2.pdf 	[server.R]
Figure 2 -- output/figs_v2/1/85/dist_br_lin.pdf 	[server.R]
Figure 3 -- output/figs_v2/random_select.pdf 		[random_start.R]
Figure 4 -- output/figs_v2/iterated_complete.pdf 	[summary.R]
Figure 5 -- output/figs_v2/conj1.pdf 				[server.R]
		 -- output/figs_v2/conj2a.pdf 				[server.R]
Figure 6 -- 
Figure 7 --
Figure 8 -- output/figs_v2/1/85/dist_br.pdf 		[server.R]
Figure 9 -- output/figs_v2/1/10/dist_br.pdf 		[server.R]
Figure 10 -- output/figs_v2/1/55/dist_br.pdf 		[server.R]
Figure 11 -- output/figs_v2/3/85/dist_br.pdf 		[server.R]
Figure 12 -- output/figs_v2/3/10/dist_br.pdf 		[server.R]
Figure 13 -- output/figs_v2/3/55/dist_br.pdf 		[server.R]
Figure 14 -- output/figs_v2/2/85/dist_br.pdf 		[server.R]
Figure 15 -- output/figs_v2/2/55/dist_br.pdf 		[server.R]
Figure 16 -- output/figs_v2/2/10/dist_br.pdf 		[server.R]
Figure 17 -- output/figs_v2/1/85/distbr_lag_lin.pdf [server.R]
Figure 18 -- output/figs_v2/distance_to_baseline.pdf [v_vecs.R]
Figure 19 -- output/figs_v2/distance_to_baseline_by_s.pdf [v_vecs.R]
Figure 20 -- output/figs_v2/quantile_summary.pdf 	[random_start.R]
Figure 21 -- output/figs_v2/1/85/main_results.pdf 	[server.R]
Figure 22 -- output/figs_v2/1/10/main_results.pdf 	[server.R]
		  -- output/figs_v2/1/55/main_results.pdf 	[server.R]
Figure 23 -- output/figs_v2/3/85/main_results.pdf 	[server.R]
Figure 24 -- output/figs_v2/3/10/main_results.pdf 	[server.R]
	      -- output/figs_v2/3/55/main_results.pdf 	[server.R]
Figure 25 -- output/figs_v2/2/85/main_results.pdf 	[server.R]
Figure 26 -- output/figs_v2/2/10/main_results.pdf 	[server.R]
	      -- output/figs_v2/2/55/main_results.pdf 	[server.R]
Figure 27 -- output/figs_v2/1/85/nc_CHE_2011.pdf 	[server.R]
	      -- output/figs_v2/1/85/nc_GRC_2009.pdf 	[server.R]
	      -- output/figs_v2/1/85/nc_HRV_2007.pdf 	[server.R]