# Branching annihilating random walk


### Installation 
- Clone the repository.
- Inside the repository, type `make`. 
  
### Running the code 
To run for a specific set of parameters, type:

    ./bws -L 255 -N 10 -h 0.4 -p 0.8 -q 0.1 > BARW/data_h_0.4_p_0.8_q_0.1.txt

where L is the length of the simulation box, N is the number of realisations and h is the ratio between the branching rate and diffusion. You can replace `data_h_0.5.txt` with any filename you'd like. 

Then you can use the python notebook `plot_data.ipynb` in the folder BARW to plot the videos, which will be saved in the same directory. For the moment, you still need to manually change the L parameter to agree with the one used in the simulation. 