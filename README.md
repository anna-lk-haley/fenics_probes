# fenics_probes
A short module for using probes with Fenics

Instantiate the probes with the command:
p = Probe(x, V, probe_id)
where x is a numpy array, V is the function space, and probe_id is a numeric id which indicates which probe is printed in each file.

To update, use:
p(u, t, probe_folder, V)
where u is a 3-component function defined on a function space, t is the time, probe_folder is the folder you want the probes to be printed to, and V is the function space.

