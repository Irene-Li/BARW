# Branching annihilating random walk


### Installation 
- Clone the repository.
- Inside the repository, type `make`. 
  
### Running the code 
To run for a specific set of parameters, type:

    ./bws -L 63 -N 10 -h 0.5 > BARW/data_h_0.5.txt

where L is the length of the simulation box, N is the number of realisations and h is the ratio between the branching rate and diffusion. You can replace `data_h_0.5.txt` with any filename you'd like. 

Then you can use the python notebook `plot_data.ipynb` in the folder BARW to plot the videos, which will be saved in the same directory. For the moment, you still need to manually change the L parameter to agree with the one used in the simulation. 

### To reproduce the cache problem 

Go into `lattice_core.c` and set `capacity_exceed = 0`. Then go through the process of running the code again and use `plot_data.ipynb` to plot movies of the realisations. You should be able to see that `movie_4.mp4` retains a patch from `movie_3.mp4`. 
