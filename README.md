Use `build_TreeSim.sh` to use `TreeSim_2.4.tar.gz` to access the **sim.rateshift.age** function for simulating trees. 

The `simulate_trees.R` produces 200 trees with one skyline interval (so uses constant rates for simulations).

`templates/ssRanges_inference_template_one_int.xml` provides the basic XML used to generate analysis files for each simulated tree.

The rates used in each simulation is saved in `true_rates.csv`. The trees are saved in `simulated_trees/200_trees.out`.

For each tree index from 0 to 199, there is a subfolder in `inf` which contains the generated XML, and resulting log and tree files.

The XML files are run using `slurm-inferences.sh`. `ssr.jar` is used to access the ssRanges BEAST2 package.

