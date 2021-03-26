import numpy as np
import sys
import random
import pandas as pd
import time

BETA = np.pi * (5/180)
epsilon = 1
AUTO = 0.1

class Ackley:
    def __init__(self, N=30, C1=20, C2=0.2, C3=2*np.pi):
        self.N = N
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
    
    def eval(self, x: list):
        first_summation = (-self.C2) * np.sqrt((1/self.N) * sum([xi*xi for xi in x]))
        second_summation = (1/self.N) * sum([np.cos(self.C3 * xi) for xi in x])
        return ((-self.C1) * np.exp(first_summation)) + self.C1 + np.e - np.exp(second_summation)

class Chromosome:
    def __init__(self, N, fun, multiple_sigmas = False, alphas = False):
        self.values = np.random.uniform(-15, 15, size=N)
        if multiple_sigmas:
            self.sigmas = np.random.uniform(1, 3, size=N)
        else:
            self.sigma = np.random.uniform(1, 3)
        if alphas:
            K = N * (N-1) // 2
            #self.alphas = [0] * K
            self.alphas = np.random.normal(0, BETA, K)
        #self.mutations = 0
        #self.successes = 0
        self.fun = fun
        self.fitness = 0

    def updateFit(self):
        self.fitness = self.fun(self.values)
        return 'fitness updated'
    
    def ratio(self):
        pass
        #return self.successes/self.mutations
        
class Optimizer:
    def __init__(self,
                generations=10000,
                popsize=30,
                children=200,
                mutation="non-correlated",
                recombination="two-mean",
                survival="comma",
                auto=False):
        self.generations = generations
        self.popsize = popsize
        self.children = children
        self.set_mutation_and_chromosomes(mutation)
        self.set_survivor_selection(survival)
        self.set_recombination(recombination)
        self.auto = auto
        self.pop = []

    # recombination functions

    def two_mean(self):
        parents = random.sample(self.pop, 2)
        new_values = [(x+y)/2 for x, y in zip(parents[0].values, parents[1].values)]
        child = self.chromosome_init()
        child.values = new_values
        if hasattr(self.pop[0], "sigma"):
            child.sigma = (parents[0].sigma + parents[1].sigma) / 2
        else:
            child.sigmas = [(s1 + s2)/2 for (s1, s2) in zip(parents[0].sigmas, parents[1].sigmas)]
        if hasattr(self.pop[0], "alphas"):
            child.alphas = [(a1 + a2)/2 for (a1, a2) in zip(parents[0].alphas, parents[1].alphas)]
        return child
    
    def multiple_mean(self):
        new_values = []
        single_sigma = hasattr(self.pop[0], "sigma")
        uses_alphas = hasattr(self.pop[0], "alphas")
        new_alphas = []
        new_sigmas = []
        for i in range(self.N):
            parents = random.sample(self.pop, 2)
            mean = (parents[0].values[i] + parents[1].values[i])/2
            new_values.append(mean)
            if not single_sigma:
                sigma_mean = (parents[0].sigmas[i] + parents[1].sigmas[i])/2
                new_sigmas.append(sigma_mean)
            if uses_alphas:
                alpha_mean = (parents[0].alphas[i] + parents[1].alphas[i])/2
                new_alphas.append(alpha_mean)
        child = self.chromosome_init()
        if uses_alphas:
            for i in range(self.N, len(self.pop[0].alphas)):
                parents = random.sample(self.pop, 2)
                if uses_alphas:
                    alpha_mean = (parents[0].alphas[i] + parents[1].alphas[i])/2
                    new_alphas.append(alpha_mean)
            child.alphas = new_alphas
        if single_sigma:
            parents = random.sample(self.pop, 2)
            child.sigma = (parents[0].sigma + parents[1].sigma) / 2
        else:
            child.sigmas = new_sigmas
        child.values = new_values
        return child

    # mutation functions

    def bounded(self, val):
        bounded_val = max(val, self.min)
        bounded_val = min(bounded_val, self.max)
        return bounded_val

    def mutation_simple(self, tau, total_pop):
        for chromosome in total_pop:
            #chromosome.mutations+=1
            new_sigma = max(epsilon, chromosome.sigma * np.exp(tau*np.random.normal(0, 1)))
            mutation_vector = np.random.normal(0, new_sigma, self.N)
            new_chromosome = [self.bounded(x + y) for x, y in zip(mutation_vector, chromosome.values)]
            new_chromosome_fit = self.fun(new_chromosome)
            if new_chromosome_fit < chromosome.fitness:
                chromosome.values = new_chromosome#.copy()
                chromosome.sigma = new_sigma
                chromosome.fitness = new_chromosome_fit
                #chromosome.successes += 1
            
        
    def mutation_noncorrelated(self, tau_global, tau_fine, total_pop):
        for chromosome in total_pop:
            mutation_global = np.random.normal(0,1)
            mutation_vector = np.random.normal(0,1,self.N)
            new_sigmas = [max(epsilon, sigmai * np.exp(tau_global*mutation_global + tau_fine*mutation_vector[i])) for i, sigmai in enumerate(chromosome.sigmas)]
            mutation_sigma = [sigmai * mutationi for sigmai, mutationi in zip(new_sigmas, mutation_vector)]
            new_chromosome = [self.bounded(x + y) for x, y in zip(mutation_sigma, chromosome.values)]
            new_chromosome_fit = self.fun(new_chromosome)
            if new_chromosome_fit < chromosome.fitness:
                chromosome.values = new_chromosome#.copy()
                chromosome.sigmas = new_sigmas#.copy()
                chromosome.fitness = new_chromosome_fit

    def mutation_correlated(self, tau_global, tau_fine, total_pop):
        for chromosome in total_pop:
            mutation_global = np.random.normal(0,1)
            mutation_vector = np.random.normal(0,1,self.N)
            new_sigmas = [max(epsilon, sigmai * np.exp(tau_global*mutation_global + tau_fine*mutation_vector[i])) for i, sigmai in enumerate(chromosome.sigmas)]
            new_alphas = [self.alphaCheck(alpha+(BETA*mutation_global)) for alpha in chromosome.alphas]
            new_C = self.Cmatrix(new_sigmas, new_alphas)
            mutation_matrix = np.random.multivariate_normal(np.zeros(self.N), new_C)
            new_chromosome = [self.bounded(x + y) for x, y in zip(mutation_matrix, chromosome.values)]
            new_chromosome_fit = self.fun(new_chromosome)
            if new_chromosome_fit < chromosome.fitness:
                chromosome.values = new_chromosome#.copy()
                chromosome.sigmas = new_sigmas#.copy()
                chromosome.alphas = new_alphas#.copy()
                chromosome.fitness = new_chromosome_fit
    
    def alphaCheck(self, alpha):
        if abs(alpha) > np.pi:
            return alpha - (2*np.pi*np.sign(alpha))
        else:
            return alpha

    def Cmatrix(self, sigmas, alphas):
        C = np.zeros((self.N, self.N))
        count=0
        for i in range(self.N):
            for j in range(self.N):
                if i==j:
                    C[i][j] = sigmas[i]*sigmas[i]
                elif j>i:
                    C[i][j] = (((sigmas[i]*sigmas[i]) * sigmas[j]*sigmas[j])/2)*np.tan(2*alphas[count])
                    count+=1
                else:
                    C[i][j] = C[j][i]
        return C

    # survival functions

    def select_survivors(self, candidates):
        candidates.sort(key=lambda c: c.fitness)
        self.pop = candidates[0:self.popsize]

    #optimization functions

    def pop_init(self):
        self.pop = []
        for _ in range(self.popsize):
            self.pop.append(self.chromosome_init())

    def optimize(self, N, fun, min_val, max_val, recording_interval = 20):
        self.N = N
        self.fun = fun
        self.pop_init()
        self.min = min_val
        self.max = max_val
        global epsilon
        orig_epsilon = epsilon
        gen_reached_e15 = 0
        per_iteration_info = []
        last_best_i = 0
        last_best = float("inf")
        init_time = time.time()
        for i in range(self.generations):
            if gen_reached_e15 == 0 and self.pop[0].fitness <= 1e-14:
                gen_reached_e15 = i
            if i%recording_interval == 0:
                per_iteration_info.append((i, self.pop[0].fitness, sum([c.fitness for c in self.pop])/self.N))
            #self.printbest(i)
            last_best, last_best_i = self.autoadapt(i, last_best, last_best_i)
            pop_children = []
            for _ in range(self.children):
                child = self.recombination()
                pop_children.append(child)
            candidates = self.surviving_parents() + pop_children
            self.mutation(candidates)
            self.select_survivors(candidates)
        best = self.pop[0]
        epsilon = orig_epsilon
        return per_iteration_info, best.fitness, gen_reached_e15, time.time() - init_time

    def printbest(self, i):
        if i%5 != 0:
            return
        print(f"{i}-rdndth generation")
        best = self.pop[0]
        print(f"best fitness: {best.fitness}")
        if hasattr(best, "sigma"):
            print(f"best sigma: {best.sigma}")
        else:
            print(f"best sigmas mean: {sum(best.sigmas)/len(best.sigmas)}")


    def autoadapt(self, i, last_best, last_best_i):
        best = self.pop[0]
        if best.fitness < last_best * 0.9:
            last_best_i = i
            last_best = best.fitness
        else:
            if i - last_best_i >= 10:
                global epsilon
                epsilon *= AUTO
                last_best_i = i
        return last_best, last_best_i


    #init functions

    def set_survivor_selection(self, survival_string):
        if survival_string == "comma":
            self.surviving_parents = lambda: []
        elif survival_string == "plus":
            self.surviving_parents = lambda: self.pop
        else:
            sys.exit("invalid survivor selection")

    def set_mutation_and_chromosomes(self, mutation_string):
        def simple_init():
            c = Chromosome(self.N, self.fun)
            c.updateFit()
            return c
        
        def noncorrelated_init():
            c = Chromosome(self.N, self.fun, multiple_sigmas = True)
            c.updateFit()
            return c
        
        def correlated_init():
            c = Chromosome(self.N, self.fun, multiple_sigmas = True, alphas = True)
            c.updateFit()
            return c

        if mutation_string == "simple":
            tau = lambda: 1/np.sqrt(self.N)
            self.mutation = lambda total_pop: self.mutation_simple(tau(), total_pop)
            self.chromosome_init = simple_init
        else:
            tau_global = lambda: 1/np.sqrt(2*self.N)
            tau_fine   = lambda: 1/np.sqrt(2*np.sqrt(self.N))
            if mutation_string == "non-correlated":
                self.mutation = lambda total_pop: self.mutation_noncorrelated(tau_global(), tau_fine(), total_pop)
                self.chromosome_init = noncorrelated_init
            elif mutation_string == "correlated":
                self.mutation = lambda total_pop: self.mutation_correlated(tau_global(), tau_fine(), total_pop)
                self.chromosome_init = correlated_init
            else:
                sys.exit("invalid mutation")
    
    def set_recombination(self, recombination_string):
        if recombination_string == "two-mean":
            self.recombination = self.two_mean
        elif recombination_string == "multiple-mean":
            self.recombination = self.multiple_mean
        else:
            sys.exit("invalid recombination")

def dumbfun(listy):
    x = listy[0]
    y = listy[1]
    return np.cos(x/500) + np.cos(x) + np.cos(y/500) + np.cos(y)

mutations = ["simple", "non-correlated", "correlated"]
recombinations = ["two-mean", "multiple-mean"]
survivals = ["comma", "plus"]

def main():
    a = Ackley(N=30)
    for mutation in mutations:
        for recomb in recombinations:
            for survival in survivals:
                data = {'history': [],
                        'best_fitness': [],
                        'generation_n': [],
                        'runtime': []}
                for _ in range(30 if mutation != "correlated" else 3):
                    opti = Optimizer(mutation=mutation, recombination=recomb, survival=survival, generations=1000)
                    history, best_fitness, generation_n, runtime = opti.optimize(30, a.eval, -15, 15)
                    data['history'].append(history)
                    data['best_fitness'].append(best_fitness)
                    data['generation_n'].append(generation_n)
                    data['runtime'].append(runtime)
                file_name = f"{mutation}_{recomb}_{survival}.csv" 
                df = pd.DataFrame(data, columns = ['history', 'best_fitness', 'generation_n', 'runtime'])
                df.to_csv(file_name, sep='\t')
                print(f'saving file: {file_name}')

main()