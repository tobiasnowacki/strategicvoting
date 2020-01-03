Figure 1 -- v_vec share paths (s = 85, lambda = 0.05)
Figure 2 -- v_vec distance (poll-to-poll) (s = 85, lambda = 0.05)
Figure 3 -- random starting points, select cases (s = 85, lambda = 0.05)
Figure 4 -- Expected benefit etc. summaries (s = 10, s = 55, s = 85; lambda = 0.05, first 60 iters)
Figure 5 -- conjectures (s = 85, lambda = 0.05; first 60 iters)
Figure 6 -- Analytical approach [Andy]
Figure 7 -- RMSE figure [Andy]
Figure 8 to 16 -- v_vec distance (poll-to-poll) all parameter combs
Figure 17 -- v_vec distance (lagged polls, s = 85, lambda = 0.05)
Figure 18 -- IRV convergence (jth poll w parameter combs vs 250th poll baseline case)
Figure 19 -- IRV convergence (jth poll w parameter combs vs 250th poll case with same s and lambda = 0.05)
Figure 20 -- IRV convergence (random starting points)
Figure 21 to 26 -- expected benefit etc. summaries (parameter combos)
Figure 27 -- case-specific voting behaviour (baseline)

common to all parameter combos:
(1) v_vec distances (although w/ linear scale for baseline)
(2) v_vec comparisons to baseline
(3) expected benefit summaries

special to baseline case:
(1) linear v_vec distances
(2) shorter iterations scale for main body
(3) conjecture tests
(4) random starting points

Write server-side script (feed params in):
	for each param combo
		- save v_vec trajectory
		- save expected benefit summaries
	if baseline
		- run conjecture tests 5(?)
		- case-specific voting behaviour

	for baseline case
		- save random starting point trajectories

Write client-side script:
	for each param combo
		- load v_vec trajectory
		- create distance plots 8 to 16
		- load expected benefit summaries
		- create exp summary plots 21 to 26

	jointly...
		- create IRV convergence comparisons 18, 19

	for baseline case:
		- create special distance plots (trajectory, linear) 1, 2
		- create special expected benefit summaries (truncated) 4
		- create lagged distance plots 17
		- create random starting point figures 3, 20
