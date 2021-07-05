# Branching annihilating random walk


### Installation 
- Clone the repository.
- Inside the repository, type `make`. 
  
### Running the code 
To run for a specific set of parameters, type:

    ./bws -L 63 -N 5 -h 0.5

where L is the length of the simulation box, N is the number of realisations and h is the ratio between the branching rate and diffusion. 

Then you can use the python notebook in the folder BARW to plot the videos, which will be saved in the same directory. For the moment, you still need to manually change the L parameter to agree with the one used in the simulation. 
