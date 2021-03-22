import numpy as np

BETA = 0.1
EPSILON = 0.05

class Ackley:
    def __init__(self, N=30, C1:=20, C2=0.2, C3=2*np.pi):
        self.N = N
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3

    def eval(self, x: list):
        first_summation = -self.C2 * np.sqrt((1/self.N) * sum(xi**2 for xi in x))
        second_summation = (1/self.N) * sum([np.cos(self.c3 * xi) for xi in x])
        np.random.normal()
        return -self.C1 * np.exp(first_summation) - np.exp(second_summation) + self.c1 + 1

class Chromosome:
    def __init__(self, N, multiple_sigmas = False, alphas = False):
        self.values = [0] * N
        if multiple_sigmas:
            self.sigmas = [0] * N
        else:
            self.sigma = 0
        if alphas:
            K = N * (N-1) / 2
            self.alphas = [0] * K
        #self.mutations = 0
        #self.successes = 0
    
    def matrix(i, j):
        pass
    
    def ratio():
        return self.successes/self.mutations
        
class Optimizer:
    def __init__(self,
                generations=10000,
                popsize=30,
                children=200,
                mutation="ancy",
                recombination=("ancy", "fancy"),
                survival="pansy",
                auto=False):
        self.generations = generations
        self.popsize = popsize
        self.children = children
        self.set_mutation_and_chromosomes(mutation)
        self.set_survivor_selection(survival)
        self.recombination = recombination
        self.auto = auto
        self.pop = []

    # mutation functions

    def mutation_simple(tau, N, fun):
        for chromosome in pop:
            #chromosome.mutations+=1
            new_sigma = choromosome.sigma * np.exp(tau*np.random.normal(0, 1))
            mutation_vector = np.random.normal(0, new_sigma, N)
            new_chromosome = [x + y for x, y in zip(mutation_vector, chromosome.values)]
            if fun(new_chromosome) < fun(chromosome.values):
                chromosome.values = new_chromosome.copy()
                chromosome.sigma = new_sigma
                #chromosome.successes += 1
            
        
    def mutation_noncorrelated(tau_global, tau_fine, N, fun):
        for chromosome in pop:
            mutation_global = np.random.normal(0,1)
            mutation_vector = np.random.normal(0,1,N)
            new_sigmas = [sigmai * np.exp(tau_global*mutation_global + tau_fine*mutation_vector[i]) for i, sigmai in enumarate(chromosome.sigmas)]
            mutation_sigma = [sigmai * mutationi for sigmai, mutationi in zip(new_sigmas, mutation_vector)]
            new_chromosome = [x+y for x, y in zip(mutation_sigma, chromosome.values)]
            if fun(new_chromosome) < fun(chromosome.values):
                chromosome.values = new_chromosome.copy()
                chromosome.sigma = new_sigmas.copy()

    def mutation_correlated(tau_global, tau_fine, N):
        pass

    # survival functions

    def survival_common(self, candidates, fun):
        candidates.sort(key=lambda c: fun(c.values))
        self.pop = candidates[0:self.popsize]

    def survival_comma(self, pop_children, fun):
        self.survival_common(pop_children, fun)

    def survival_plus(self, pop_children, fun):
        self.survival_common(self.pop + pop_children, fun)

    #optimization functions

    def pop_init(self)
        self.pop = []
        for _ in range(self.popsize):
            self.pop.append(self.chromosome_init())

    def optimize(self, fun, N):
        self.pop_init()
        for _ in range(self.generations):
            pop_children = [] # create children through mutation and recombination
            self.select_survivors(pop_children, fun)
        best = min(self.pop, key=lambda chrom: fun(chrom))
        return best.values

    #init functions

    def set_survivor_selection(self, survival_string):
        if survival_string == "ancy":
            self.select_survivors = self.survival_comma
        elif survival_string == "fancy":
            self.select_survivors = self.survival_plus

    def set_mutation_and_chromosomes(self, mutation_string):
        if mutation_string == "ancy":
            tau = 1/np.sqrt(N)
            self.mutation = lambda: mutation_simple(tau)
            self.chromosome_init = lambda: Chromosome(N)
        else
            tau_global = 1/np.sqrt(2*N)
            tau_fine   = 1/np.sqrt(2*np.sqrt(N))
            if mutation_string == "fancy":
                self.mutation = lambda: mutation_noncorrelated(tau_global, tau_fine)
                self.chromosome_init = lambda: Chromosome(N, multiple_sigmas = True)
            elif mutation_string == "pancy":
                self.mutation = lambda: mutation_correlated(tau_global, tau_fine)
                self.chromosome_init = Chromosome(N, multiple_sigmas = True, alphas = True)
            else:
                sys.exit("invalid mutation")