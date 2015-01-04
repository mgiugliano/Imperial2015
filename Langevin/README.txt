
This directory contains files to reproduce results presented in the
paper:

Linaro, D., Storace, M., and Giugliano, M.  "Accurate and fast
simulation of channel noise in conductance-based model neurons by
diffusion approximation".

The subdirectories contain the following items.

1. C++ contains the source files for simulating the open-close
kinetics of sodium and potassium ion channels. The programs allow to
simulate the models discussed in the paper in the condition of voltage
clamp.

2. mod-files contains the mod-files that can be used in NEURON to
simulate the full models.

3. python contains a script, HHneuron.py, that allows to simulate the
models implemented through mod-files.

4. matlab contains some auxiliary Matlab scripts that can be used to
read data files and produce raster plots.

5. supplem contains the Matlab scripts employed to generate the figures and analysis,
included as Supplemental Material.

For any question, feel free to contact: daniele.linaro@unige.it, marco.storace@unige.it, and michele.giugliano@ua.ac.be
