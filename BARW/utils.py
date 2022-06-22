import numpy as np


def read_file(file, verbose=False): 

	realisation_flag = False 
	moments_flag = False

	realisations = [] 
	moments = [] 
	n = 0 

	line = file.readline() 
	while line: 
		if line[0] == '#':
			if verbose: 
				print(line)
			if line.startswith('# Starting'):
				realisation_flag = True 
			if line.startswith('#Writing'): 
				moments_flag = True 
		else: 
			if moments_flag: 
				moments.append(line)
			elif realisation_flag: 
				new_realisation = [line]
				realisations.append(new_realisation)
				realisation_flag = False
			else:
				new_realisation.append(line)   
		line = file.readline()

	return realisations, moments 


def extract_lite(realisation, N): 
    edge_reach = [] 
    msd = [] 
    time = []  

    for (i, line) in enumerate(realisation):  
        if line.startswith('time'):
            time.append(float(line[5:-1]))
        elif line.startswith('msd'): 
            msd.append(float(line[4:-1]))
        elif line.startswith('edge reach'):
            edge_reach = np.fromstring(line[11:-1], sep=',')
        elif line.startswith('coarse_grain_tracer'):
            f = lambda x: np.fromstring(x, sep=',')
            coarse_grain_moments = np.vstack(list(map(f, realisation[i+1:i+N+1])))
    return time, msd, edge_reach, coarse_grain_moments

def extract_moments(moments): 
    active_moments = [] 
    tracer_moments = [] 
    active_flag = False 
    tracer_flag = False 
    
    for line in moments: 
        if line.startswith(' active moments'):
            active_flag = True 
        elif line.startswith('tracer moments'):
            active_flag = False 
            tracer_flag = True 
        else: 
            if active_flag: 
                active_moments.append(np.fromstring(line, sep='\t'))
            if tracer_flag: 
                tracer_moments.append(np.fromstring(line, sep='\t'))
    return np.vstack(active_moments), np.vstack(tracer_moments)
        
def extract_evolution(realisation, L): 
    tracers = []
    active_particles = [] 
    edges = [] 
    msds = [] 
    total_active = [] 
    total_tracer = [] 
    time = [] 
    passive_flag = False 
    active_flag = False
    for line in realisation:  
        if line.startswith('time'):
            time.append(float(line[5:-1]))
        elif line.startswith('total active'):
            total_active.append(float(line[13:-1]))
        elif line.startswith('total tracer'):
            total_tracer.append(float(line[13:-1]))
        elif line.startswith('passive'):
            passive_indices = np.fromstring(line[8:-1], sep=',')
            passive_indices = np.unravel_index(passive_indices.astype('int'), (L, L))
            passive_flag = True 
        elif line.startswith('active'):
            active_indices = np.fromstring(line[7:-1], sep=',')
            active_indices = np.unravel_index(active_indices.astype('int'), (L, L))
            active_flag = True 
        elif line.startswith('msd'):
            msds.append(float(line[4:-1]))
        elif line.startswith('edge'): 
            edges.append(np.fromstring(line[5:-1], sep=','))
        else: 
            tracer_snapshot = np.zeros((L, L)) # clear the snapshot for the next time slice 
            active_snapshot = np.zeros((L, L)) 
            
        if passive_flag and active_flag: 
            tracer_snapshot[passive_indices] = 1 # set occupied sites to -1 
            
            for (i,j) in zip(*active_indices):
                active_snapshot[i, j] += 1 
            tracers.append(tracer_snapshot)
            active_particles.append(active_snapshot)
            passive_flag = False 
            active_flag = False 
            
    return np.array(edges), np.array(msds), np.array(counts), np.array(tracers), np.array(active_particles)

def make_movies(evolution, label):
    fig = plt.figure(figsize=(20, 20))
    ims = []
    plt.axis('off')
    for xy in evolution:
        im = plt.imshow(xy, animated=True, vmin=0, vmax=2, cmap='Blues', origin='lower') 
        plt.axis('off')
        ims.append([im])
        
    ani = am.ArtistAnimation(fig, ims, interval=100, blit=True,
                                    repeat_delay=1000)
    mywriter = am.FFMpegWriter()
    ani.save("{}_movie_{}.mp4".format(label, i), writer=mywriter)
    plt.close()  
    
    
def pad(nested_lists): 
    max_length = max(map(len, nested_lists))
    for l in nested_lists: 
        length = len(l)
        if length == 0: 
            raise Exception('no time series data recorded')
        elif length < max_length: 
            l.extend([l[-1]]*(max_length-length))
    return np.vstack(nested_lists)

def select_longest(times): 
    return max(times, key=lambda x:len(x))
    
    


