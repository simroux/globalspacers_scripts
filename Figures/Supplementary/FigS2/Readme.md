# Step 1 - get metrics from the simulation results
`./get_metrics.pl run`
Generate the summary file from the simulations used to plot the figures, needs as input the tsv file "Compare_spacer_recovery_sub200.tsv", see "Analyses/Extra/Tools_benchmark/" for how to generate.
# Step 2 - use R to generate the figure
`plot_fig_benchmarks.R`
To plot the figures in R