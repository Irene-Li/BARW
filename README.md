# Branching annihilating random walk

This is adapted from the code used in [1].  

### Installation 
- Clone the repository.
- Inside the repository, type `make`. 
  
### Running the code 
To run for a specific set of parameters, type:

    ./bws -L 2049 -N 1000 -h 0.308 -p 0.34 -q 0.33 > BARW/data_L2049_N1000_p0.34_q0.33_h0.308.dat

where L is the length of the simulation box, N is the number of realisations and h is the ratio between the branching rate and diffusion. You can replace `data_h_0.5.txt` with any filename you'd like. 

Then you can use the python notebook `plot_data.ipynb` in the folder BARW to plot the videos, which will be saved in the same directory. For the moment, you still need to manually change the L parameter to agree with the one used in the simulation. 

### References 
[1] : Bordeu, I., Amarteifio, S., Garcia-Millan, R., Walter, B., Wei, N. and Pruessner, G., 2019. Volume explored by a branching random walk on general graphs. Scientific reports, 9(1), pp.1-9.
