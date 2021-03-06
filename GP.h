/*  Copyright © 2012 Mauro Castelli
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.
 
 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//!  \file            GP.h
//! \brief            file containing the definition of the classes used to represent a symbol, a node of the tree, a population of individuals and definition of the functions
//! \date            created on 01/09/2012

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <tuple>
#include <limits>
#include <time.h>
/// Macro used to generate a random number
#define frand() ((double) rand() / (RAND_MAX))


/// variable containing the numbers of terminal symbols (variables, not constant values)
int NUM_VARIABLE_SYMBOLS;
/// variable that stores the numbers of terminal symbols that contain constant values
int NUM_CONSTANT_SYMBOLS;
/// variable containing the numbers of functional symbols
int NUM_FUNCTIONAL_SYMBOLS;

/// struct used to store a single instance in memory
typedef struct Instance_{
	/// array containing the values of the independent variables
	double *vars;
	/// variable that stores the result of the evaluation of an individual on the particular instance represented by the values stored in the variable double *vars
	double res;
	/// target value
	double y_value;
} Instance;

/// variable used to store training and test instances in memory
Instance *set;

/// variable containing the numbers of rows (instances) of the training dataset
int nrow;
/// variable containing the numbers of columns (excluding the target) of the training dataset
int nvar;
/// variable containing the numbers of rows (instances) of the validation dataset (used to identify overfitting)
int nrow_val;
/// variable containing the numbers of rows (instances) of the test dataset
int nrow_test;
/// variable containing the numbers of columns (excluding the target) of the test dataset
int nvar_test;

using namespace std;

/// struct used to store the parameters of the configuration.ini file
typedef struct cfg_{
	/// size of the population: number of candidate solutions
	int population_size;
	/// number of iterations of the GP algorithm
	int max_number_generations;
	/// initialization method: 0 -> grow method, 1-> full method, 2-> ramped half and half method
	int init_type;
	/// crossover rate
	double p_crossover;
	/// mutation rate
	double p_mutation;
	/// maximum depth of a newly created individual
	int max_depth_creation;
	/// size of the tournament selection
	int tournament_size;
	/// variable that indicates if it is possible to accept single-node individual in the initial population
	int zero_depth;
	/// mutation step of the geometric semantic mutation
	double mutation_step;
	/// variable that indicates the number of constants to be inserted in the terminal set
	int num_random_constants;
	/// variable that indicates the minimum possible value for a random constant
	double min_random_constant;
	/// variable that indicates the maximum possible value for a random constant
	double max_random_constant;
	/// variable that indicates if the problem is a minimization problem (1) or a maximization problem (0)
	int minimization_problem;
	/// variable that defines the proportion of the validation set in percent of the training set
	double validation_set_size;
	/// variable that defines the minimum number of generations to run in the beginning without using the repulsors
	double repulsor_min_age;
	/// variable that defines the maximum number of repulsors to keep during execution
	int semantic_repulsor_max_number;
	/// variable that defines the numbers of individuals to keep as elite on the validation set
	int validation_elite_size;
	/// variable that indicates whether and how individuals will be checked for overfitting (always using the n best of the population)
	int use_best_as_rep_candidate;
	/// variable that indicates whether overfitting is determined by the median or the average of the validation elites fitness
	int overfit_by_median;
	/// variable that indicates whether the split of the training set into validation and training data should be shuffled (different every time)
	int shuffle_validation_split;
	/// variable that indicates whether to log semantics of individuals and repulsors to an output file
	int log_semantics;
	/// variable that indicates whether to use repulors as separate objectives or to aggregate them into one single objective
	int aggregate_repulsors;
	/// variable that indicates whether force recreation of individuals, if they are equal to any repulsor
	int force_avoid_repulsors;
	/// variable that indicates whether the fraction of distance, below which two individuals are considered equal
	double equality_delta;
	/// variable that indicates whether to randomly select from the best pareto front during tournament selection
	int true_pareto_selection;
	/// variable that indicates whether to exclude fitness from domination determination (as long as repuslers aren't aggregated)
	int domination_exclude_fitness;
	/// variable that indicates whether to merge repulsers rather than replace
	int merge_repulsors;
}cfg;

/// struct variable containing the values of the parameters specified in the configuration.ini file
cfg config;



/**
 * \class symbol
 *
 * \brief
 *
 * This class represents a symbol of the set T (terminal symbols) or F (functional symbols).
 *
 * \author Mauro Castelli
 *
 * \version 0.0.1
 *
 * \date 01/09/2012
 *
 */
class symbol{
public:
	/// boolean variable used to discriminate between functional and terminal symbols
	bool type;
	/// int variable that contains the number of arguments accepted by a symbol. It is 0 for a terminal symbol
	int arity;
	/// int variable that contains a unique identifier for the symbol
	int id;
	/// symbolic name of the symbol
	char name[30];
	/// variable that contains the current value of a terminal symbol
	double value;
	///default constructor
	symbol(){};
	///constructor that creates a symbols and sets the variable type, num_arguments, id and value
	symbol(bool p_type, int p_arity, int p_id, const char *p_name){
		type=p_type;
		arity=p_arity;
		id=p_id;
		strcpy(name,p_name);
	};
};

///array containing terminal and functional symbols
vector <symbol *> symbols;

/**
 * \class node
 *
 * \brief
 *
 * This class is used to represent a node of the tree.
 *
 * \author Mauro Castelli
 *
 * \version 0.0.1
 *
 * \date 01/09/2012
 *
 */
class node{
public:
	///symbol inside a node
	symbol* root;
	/// parent node
	node* parent;
	///pointers to children
	node **children;
	/// class destructor
	~node() {delete[] children;}
};

/// tuple containing a fitness measure and a non-dominated sorting rank
typedef tuple<double, int> fitness_data;
/// vector of fitness data
typedef vector< fitness_data > fitness_list;

/**
 * \class population
 *
 * \brief
 *
 * This class is used to represent a GP population.
 *
 * \author Mauro Castelli
 *
 * \version 0.0.1
 *
 * \date 01/09/2012
 *
 */
class population{
public:
	/// pointers to individuals
	node **individuals;
	/// int variable that contains the index of the best individual in the population
	int index_best;
	/// int variable that contains the number of individuals that are inside the population
	int num_ind;
	/// array of training fitness data
	fitness_data *fitness;
	/// array of validation fitness data
	fitness_data *fitness_val;
	/// array of test fitness data
	fitness_data *fitness_test;
	/// class constructor
	population(){
		individuals=new node* [config.population_size];
		num_ind=0;
		fitness=new fitness_data [config.population_size];
		fitness_val=new fitness_data [config.population_size];
		fitness_test=new fitness_data [config.population_size];
	};
	/// class destructor
	~population() { delete[] individuals; }
};

/// array where each element (that is a 3-tuple) contains the fitness, the nondominated sorting rank and the crowded distance measure on the training set at generation g
fitness_list fit_;
/// array where each element (that is a 3-tuple) contains the fitness, the nondominated sorting rank and the crowded distance measure on the validation set at generation g
fitness_list fit_val;
/// array where each element (that is a 3-tuple) contains the fitness, the nondominated sorting rank and the crowded distance measure on the test set at generation g
fitness_list fit_test;
/// array where each element (that is a 3-tuple) contains the fitness, the nondominated sorting rank and the crowded distance measure on the training set at generation g+1
fitness_list fit_new;
/// array where each element (that is a 3-tuple) contains the fitness, the nondominated sorting rank and the crowded distance measure on the validation set at generation g+1
fitness_list fit_new_val;
/// array where each element (that is a 3-tuple) contains the fitness, the nondominated sorting rank and the crowded distance measure on the test set at generation g+1
fitness_list fit_new_test;

/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the training set at generation g
vector < vector<double> > sem_train_cases;
/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the training set at generation g+1
vector < vector<double> > sem_train_cases_new;
/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the validation set at generation g
vector < vector<double> > sem_val_cases;
/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the validation set at generation g+1
vector < vector<double> > sem_val_cases_new;
/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the test set at generation g
vector < vector<double> > sem_test_cases;
/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the test set at generation g+1
vector < vector<double> > sem_test_cases_new;

/// array where each element (that is also an array) contains the semantics of a repulsor individual
vector < vector<double> > sem_repulsors;
/// array where each element (that is also an array) contains the semantics of a repulsor individual that will be additionally checked in the next generations
vector < vector<double> > sem_repulsors_new;
/// array where each element (that is also an array) contains the distances of an individual each repulsor at generation g
vector < vector<double> > repulsor_distances;
/// array where each element (that is also an array) contains the distances of an individual each repulsor at generation g+1
vector < vector<double> > repulsor_distances_new;
// array containing the severity of overfitting of the repulsors (based on the fitness difference to the validation elite)
vector <double> overfit_severity;
// array containing the severity of overfitting of the repulsors (based on the fitness difference to the validation elite), for the next generation
vector <double> overfit_severity_new;
// array containing a mapping of an id, to an index in the repulsor tables and a marker of the generation, to keep track of changing positions for logging
vector <tuple<int, int, int>> repulsor_map;

// array containing the semantics of the elite n individuals on the validation set (generation independent list of best individuals)
vector < vector<double> > sem_val_elite;
// array containing the fitness of the elite n individuals on the validation set
vector <double> fit_val_elite;
// index of the worst individual in the array of elite n individuals on the validation set
int val_elite_worse_idx;
// average fitness of the elite n individualson the validation set
double val_elite_avg_fit;

/// variable that stores the index of the best individual.
int index_best;
/// variable that stores the maximum distance between individuals in the population.
double population_max_distance;
double population_combined_max_distance;

/*!
 * \fn                 void read_config_file(cfg *config)
 * \brief             function that reads the configuration file
 * \param          cfg *config: pointer to the struct containing the variables needed to run the program
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void read_config_file(cfg *config, char *file);


/*!
 * \fn                void create_T_F()
 * \brief             function that creates the terminal and functional sets.
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void create_T_F();


/*!
 * \fn                 int choose_function()
 * \brief             function that randomly selects a functional symbol
 * \return           int: the ID of the chosen functional symbol
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
int choose_function();

/*!
 * \fn                 int choose_terminal()
 * \brief             function that randomly selects a terminal symbol. With probability 0.7 a variable is selected, while random constants have a probability of 0.3 to be selected. To change these probabilities just change their values in the function.
 * \return           int: the ID of the chosen terminal symbol
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
int choose_terminal();


/*!
 * \fn                void create_grow_pop(population **p)
 * \brief             function that creates a population using the grow method.
 * \param          population **p: pointer to an empty population
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void create_grow_pop(population **p);

/*!
 * \fn                void create_full_pop(population **p)
 * \brief             function that creates a population of full trees (each tree has a depth equal to the possible maximum length).
 * \param          population **p: pointer to an empty population
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void create_full_pop(population** p);

/*!
 * \fn                void create_ramped_pop(population **p)
 * \brief             function that creates a population with the ramped half and half algorithm.
 * \param          population **p: pointer to an empty population
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void create_ramped_pop(population **p);

/*!
 * \fn                void create_population(population **p, int i)
 * \brief             function that creates a population using the method specified by the parameter int i.
 * \param          population **p: pointer to an empty population
 * \param          int i: type of initialization method
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void create_population(population **p, int i);


/*!
 * \fn                void create_grow_tree(node **el, int depth, node *parent, int max_depth)
 * \brief             function that creates a random tree with depth in the range [0;max_depth]
 * \param          node **el: pointer to the node that must be added to the tree
 * \param          int depth: current depth of the tree
 * \param          node *parent: parent node
 * \param          int max_depth: maximum depth of the tree
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void create_grow_tree(node **el, int depth, node *parent, int max_depth);

/*!
 * \fn                void create_full_tree(node **el, int depth, node *parent, int max_depth)
 * \brief             function that creates a tree with depth equal to the ones specified by the parameter max_depth
 * \param          node **el: pointer to the node that must be added to the tree
 * \param          int depth: current depth of the tree
 * \param          node *parent: parent node
 * \param          int max_depth: maximum depth of the tree
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void create_full_tree(node **el, int depth, node *parent, int max_depth);


/*!
 * \fn                 double protected_division(double num, double den)
 * \brief             function that implements a protected division. If the denominator is equal to 0 the function returns 1 as a result of the division;
 * \param          double num: numerator
 * \param          double den: denominator
 * \return           double: the result of the division if denominator is different from 0; 1 if denominator is equal to 0
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
double protected_division(double num, double den);


/*!
 * \fn                 double eval(node *tree)
 * \brief             function that evaluates a tree.
 * \param          node *tree: radix of the tree to be evaluated
 * \return           double: the value of the evaluation
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
double eval(node *tree);


/*!
 * \fn                 void evaluate(population **p)
 * \brief             function that calculates the fitness of all the individuals and determines the best individual in the population
 * \param          population **p: pointer to the population containing the individuals
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void evaluate(population **p);


/*!
 * \fn                 double Myevaluate(node *el)
 * \brief             function that calculates the training fitness of an individual (representing as a tree)
 * \param          node *el: radix of the tree
 * \return           double: the training fitness of the individual
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
double Myevaluate(node *el);

/*!
 * \fn                 double Myevaluate_val(node *el)
 * \brief             function that calculates the validation fitness of an individual (representing as a tree)
 * \param          node *el: radix of the tree
 * \return           double: the validation fitness of the individual
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
double Myevaluate_val(node *el);


/*!
 * \fn                 double Myevaluate_test(node *el)
 * \brief             function that calculates the test fitness of an individual (representing as a tree)
 * \param          node *el: radix of the tree
 * \return           double: the test fitness of the individual
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
double Myevaluate_test(node *el);


/*!
 * \fn                void Myevaluate_random (node *el, vector <double> & sem)
 * \brief             function that calculates the semantics (considering training instances) of a randomly generated tree. The tree is used to perform the semantic geometric crossover or the geometric semantic mutation
 * \param          node* el: radix of the tree to be evaluated
 * \param          vector <double> & sem: reference to an empty vector that, at the end of the function, will contain the semantics of the individual calculated on the training instances
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file              GP.h
 */
void Myevaluate_random (node *el, vector <double> & sem);

/*!
 * \fn                void Myevaluate_random_val (node *el, vector <double> & sem)
 * \brief             function that calculates the semantics (considering validation instances) of a randomly generated tree. The tree is used to perform the semantic geometric crossover or the geometric semantic mutation
 * \param          node* el: radix of the tree to be evaluated
 * \param          vector <double> & sem: reference to an empty vector that, at the end of the function, will contain the semantics of the individual calculated on the validation instances
 * \return           void
 * \date             TODO add date
 * \author          Paul Englert
 * \file              GP.h
 */
void Myevaluate_random_val (node *el, vector <double> & sem);


/*!
 * \fn                 double Myevaluate_random_test(node *el, , vector <double> & sem)
 * \brief             function that calculates the semantics (considering test instances) of a randomly generated tree. The tree is used to perform the semantic geometric crossover or the geometric semantic mutation
 * \param          node* el: radix of the tree to be evaluated
 * \param          vector <double> & sem: reference to an empty vector that, at the end of the function, will contain the semantics of the individual calculated on the test instances
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file              GP.h
 */
void Myevaluate_random_test(node *el, vector <double> & sem);

/*!
 * \fn                 void nsga_II_sort(population **p)
 * \brief             function that assigns a non-domination rank and a crowded distance measure to the individuals fitness data of a population
 * \param          population **p: pointer to the population containing the individuals
 * \return           void
 * \date             TODO add date
 * \author          Paul Englert
 * \file              GP.h
 */
void nsga_II_sort(population **p);


/*!
 * \fn                 void perform_fast_non_domination_sort(population **p, vector<int> *d_front, int **d_counts, vector< vector<int> > *d_individuals)
 * \brief             function that identifies non-domination levels, domination counts, and dominated individuals as needed for NSGA-II
 * \param          population **p: pointer to the population containing the individuals
 * \param          vector<int> *d_front: pointer to a vector of indexes that will be filled with all indexes lying on the first pareto level
 * \param          int **d_counts: pointer to an array of counts, where each position corresponds to an individual in the population, it will be filled with the count that each individual has been dominated by another indivdual
 * \param          vector< vector<int> > *d_individuals: pointer to a vector corresponds to an individual in the population and each element is a vector of indexes, that the individual dominates
 * \return           void
 * \date             TODO add date
 * \author          Paul Englert
 * \file              GP.h
 */
void perform_fast_non_domination_sort(population **p, vector<int> *d_front, int **d_counts, vector< vector<int> > *d_individuals);

/*!
 * \fn                 vector <int>* extract_next_front(population **p, vector<int> *d_front, int **d_counts, vector< vector<int> > *d_individuals)
 * \brief             function that extracts the next deeper level of domination front
 * \param          int cur_front: index denoting the parameter *d_fronts level in the pareto hierarchie
 * \param          population **p: pointer to the population containing the individuals
 * \param          vector <int>* next_front: pointer to a vector of indexes representing the next front (one level deeper)
 * \param          vector<int> *d_front: pointer to a vector of the indexes of the parent front
 * \param          int **d_counts: pointer to an array of counts, where each position corresponds to an individual in the population and the value determines the number of dominating individuals
 * \param          vector< vector<int> > *d_individuals: pointer to a vector corresponds to an individual in the population and each element is a vector of indexes, that the individual dominates
 * \return           void
 * \date             TODO add date
 * \author          Paul Englert
 * \file              GP.h
 */
void extract_next_front(int cur_front, vector <int>* next_front, population **p, vector<int> *d_front, int **d_counts, vector< vector<int> > *d_individuals);

/*!
 * \fn                bool dominates(int i, int j)
 * \brief             function that determines whether individual i dominates individual j
 * \param            int i: individual i to check
 * \param            int j: individual j to check as opponent
 * \return           bool: whether i dominates j
 * \date             TODO add date
 * \author          Paul Englert
 * \file               GP.h
 */
bool dominates(int i, int j);
/*!
 * \fn                void update_terminal_symbols(int i)
 * \brief             function that updates the value of the terminal symbols in a tree.
 * \param          int i: line of the dataset containing the values of the terminal symbols
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void update_terminal_symbols(int i);


/*!
 * \fn                 void delete_individual(node * el)
 * \brief             function that deletes a tree and frees the the memory allocated to store the tree
 * \return           node* el: radix of the tree to be deleted
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void delete_individual(node * el);


/*!
 * \fn                int tournament_selection()
 * \brief             function that implements a tournament selection procedure
 * \return           int: index of the best individual among the ones that participate at the tournament
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
int tournament_selection();


/*!
 * \fn                void reproduction(int i)
 * \brief             function that copy an individual of the population at generation g-1 to the current population(generation g)
 * \param            int i: index of the individual in the current population
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void reproduction(int i);


/*!
 * \fn                void geometric_semantic_crossover(int i)
 * \brief             function that performs a geometric semantic crossover
 * \param            int i: index of the newly created individual in the new population
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void geometric_semantic_crossover(int i);

/*!
 * \fn                void geometric_semantic_mutation(int i)
 * \brief             function that performs a geometric semantic mutation
 * \param            int i: index of the mutated individual in the new population
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file             GP.h
 */
void geometric_semantic_mutation(int i);


/*!
 * \fn                void update_training_fitness(vector <double> semantic_values, bool crossover)
 * \brief             function that calculate the training fitness of an individual using the information stored in its semantic vector. The function updates the data structure that stores the training fitness of the individuals
 * \param            vector <double> semantic_values: vector that contains the semantics (calculated on the training set) of an individual
 * \param            bool crossover: variable that indicates if the function has been called by the geometric semantic crossover or by the geometric semantic mutation
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void update_training_fitness(vector <double> semantic_values, bool crossover);



/*!
 * \fn                void update_validation_fitness(vector <double> semantic_values, bool crossover)
 * \brief             function that calculate the validation fitness of an individual using the information stored in its semantic vector. The function updates the data structure that stores the validation fitness of the individuals
 * \param            vector <double> semantic_values: vector that contains the semantics (calculated on the validation set) of an individual
 * \param            bool crossover: variable that indicates if the function has been called by the geometric semantic crossover or by the geometric semantic mutation
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void update_validation_fitness(vector <double> semantic_values, bool crossover);


/*!
 * \fn                void update_test_fitness(vector <double> semantic_values, bool crossover)
 * \brief             function that calculate the test fitness of an individual using the information stored in its semantic vector. The function updates the data structure that stores the test fitness of the individuals
 * \param            vector <double> semantic_values: vector that contains the semantics (calculated on the test set) of an individual
 * \param            bool crossover: variable that indicates if the function has been called by the geometric semantic crossover or by the geometric semantic mutation
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void update_test_fitness(vector <double> semantic_values, bool crossover);


/*!
 * \fn                void update_repulsor_distances(vector <double> semantic_values, bool crossover)
 * \brief             function that calculates the semantic distance of an indivdual for all existing repulsors. The function updates the structure that stores the distances of each individual to each repulsor
 * \param            vector <double> semantic_values: vector that contains the semantics (calculated on the training set) of an individual
 * \param            bool crossover: variable that indicates if the function has been called by the geometric semantic crossover or by the geometric semantic mutation
 * \return           void
 * \date             TODO add date
 * \author          Paul Englert
 * \file               GP.h
 */
void update_repulsor_distances(vector <double> semantic_values, bool crossover);


/*!
 * \fn                void update_validation_elite(vector <double> semantic_values, double fitness)
 * \brief             function that checks the given semantics and adds them to the elite group on the validation set, if acceptable
 * \param            vector <double> semantic_values: vector that contains the semantics (calculated on the validation set) of an individual
 * \param            double fitness: variable that contains the fitness of the individual
 * \return           void
 * \date             TODO add date
 * \author          Paul Englert
 * \file               GP.h
 */
void update_validation_elite(vector <double> semantic_values, double fitness);

/*!
 * \fn               double is_overfitting(double fit)
 * \brief             function that calculates whether an individual with fitness fit is overfitting, based on the difference in validation and training data and the configured threshold.
 * \param            double fit: fitness of the individual to check
 * \return           bool: if the individual i is overfitting
 * \date             TODO add date
 * \author          Paul Englert
 * \file               GP.h
 */
bool is_overfitting(double fit);

/*!
 * \fn               double get_overfitting_severity(double fit)
 * \brief             function that returns the difference between the average validation elite fitness and the individuals fitness on the validation data
 * \param            double fit: fitness of the individual to check
 * \return           double: margin by which the individual is worse than the average validation elite fitness
 * \date             TODO add date
 * \author          Paul Englert
 * \file               GP.h
 */
double get_overfitting_severity(double fit);

/*!
 * \fn               void add_repulsor(vector <double> semantics, double validation_fitness)
 * \brief             function that adds the given semantics as a repulsor if it isn't already recorded
 * \param            vector <double> semantics: semantics of the repulsor
 * \param            double validation_fitness: fitness of the repulsor on the validation data
 * \return           void
 * \date             TODO add date
 * \author          Paul Englert
 * \file               GP.h
 */
void add_repulsor(vector <double> semantics, double validation_fitness);

/*!
 * \fn                int best_individual()
 * \brief             function that finds the best individual in the population
 * \return           int: the index of the best individual
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
int best_individual();


/*!
 * \fn               void update_tables()
 * \brief            function that updates the tables used to store fitness values and semantics of the individual. It is used at the end of each iteration of the algorithm
 * \return           void
 * \date             01/09/2012
 * \author           Mauro Castelli
 * \file             GP.h
 */
void update_tables();


/*!
 * \fn               int update_repulsors()
 * \brief            function that updates the tables used to store the semantics of the repulsors, as well as recalculates the distances of the individuals to the repulsors
 * \return           int: the number of repulsors lost due to size constaint from configuration
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
int update_repulsors();

/*!
 * \fn               int merge_repulsors(int count)
 * \brief            function that merges n repulsers in the currently active repulser table
 * \param 			 int count: number of repulsers to recursively merge
 * \return           int: index of the first freed slot in the repusler table
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
int merge_repulsors(int count);


/*!
 * \fn               void read_input_data(char *train_file, char *test_file)
 * \brief            function that reads the data from the training file and from the test file.
 * \return           void
 * \date             01/09/2012
 * \author           Mauro Castelli
 * \file             GP.h
 */
void read_input_data(char *train_file, char *test_file);


/*!
 * \fn               bool better (double f1, double f2)
 * \brief            function that compare the fitness of two solutions.
 * \param            double f1: fitness value of an individual
 * \param            double f2: fitness value of an individual
 * \return           bool: true if f1 is better than f2, false in the opposite case
 * \date             01/09/2012
 * \author           Mauro Castelli
 * \file             GP.h
 */
bool better (double f1, double f2);

/*!
 * \fn               bool nsga_II_better (tuple<double, int, double> i1, tuple<double, int, double> i2)
 * \brief            function that compares two solutions.
 * \param            tuple<double, int, double> i1: fitness, pareto rank and crowded distance value of an individual
 * \param            tuple<double, int, double> i2: fitness, pareto rank and crowded distance value of an individual
 * \return           bool: true if i1 is better than i2, false in the opposite case
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
bool nsga_II_better (fitness_data i1, fitness_data i2);

/*!
 * \fn               void calculateMaxDistance ()
 * \brief            function that calculates the maximum distance in the population.
 * \return           void
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
void calculateMaxDistance ();

/*!
 * \fn               bool isEqualToAnyRepulsor (int i)
 * \brief            function that compares an individuals semantics to all repulsors.
 * \param            int i: index of the individual to test
 * \return           bool: true if ind is equal to any repulsor (based on equality_delta)
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
bool isEqualToAnyRepulsor (int i);

/*!
 * \fn               double calculateSmallestRepulsorDistance (int i)
 * \brief            function that returns the distance to the closest repulsor.
 * \param            int i: index of the individual to test
 * \return           double: distance to closest repulsor
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
double calculateSmallestRepulsorDistance (int i);

/*!
 * \fn               void log_semantics (ofstream *csem, int num_gen)
 * \brief            function that logs all semantics of the current population and the repulsors.
 * \param 			 ofstream csem: the output stream to write the semantics to
 * \param 			 int num_gen: the current generation number
 * \return           void
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
void log_semantics (ofstream *csem, int num_gen);


void read_config_file(cfg *config, char *file){
	clog<<"\tReading Configuration:"<<endl;
	fstream f(file, ios::in);
	if (!f.is_open()) {
    		cerr<<"CONFIGURATION FILE NOT FOUND." << endl;
    		exit(-1);
	}
	while(!f.eof()){
		char in_str[100]="";
		char str1[100]="";
		char str2[100]="";
		int j=0;
		f.getline(in_str,100);
		string line = string(in_str);
		clog<<"\t\t"<<in_str<<endl;
		if(in_str[0]!='\0'){
			// secure the format
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\t'), line.end());
			strcpy(in_str, line.c_str());
			while(in_str[j]!='='){
				str1[j] = in_str[j];
				j++;
			}
			j++;
			int i=0;
			while(in_str[j]!='\0'){
				str2[i] = in_str[j];
				j++;
				i++;
			}
		}
		if(strcmp(str1, "population_size") == 0)
			config->population_size = atoi(str2);
		if(strcmp(str1, "max_number_generations") == 0)
			config->max_number_generations=atoi(str2); 
		if(strcmp(str1, "init_type") == 0)
			config->init_type=atoi(str2);
		if(strcmp(str1, "p_crossover") == 0)
			config->p_crossover=atof(str2);
		if(strcmp(str1, "p_mutation") == 0)
			config->p_mutation=atof(str2);
		if(strcmp(str1, "max_depth_creation") == 0)	
			config->max_depth_creation=atoi(str2);
		if(strcmp(str1, "tournament_size") == 0)	
			config->tournament_size=atoi(str2);
		if(strcmp(str1, "zero_depth") == 0)
			config->zero_depth=atoi(str2);
		if(strcmp(str1, "mutation_step") == 0)
			config->mutation_step=atof(str2);
		if(strcmp(str1, "num_random_constants") == 0){
			config->num_random_constants=atoi(str2);	
			NUM_CONSTANT_SYMBOLS=config->num_random_constants;
        }
        if(strcmp(str1, "min_random_constant") == 0)
			config->min_random_constant=atof(str2);
		if(strcmp(str1, "max_random_constant") == 0)
			config->max_random_constant=atof(str2);
		if(strcmp(str1, "minimization_problem") == 0)
			config->minimization_problem=atoi(str2);
		if(strcmp(str1, "log_semantics") == 0)
			config->log_semantics=atoi(str2);
		if(strcmp(str1, "validation_set_size") == 0)
			config->validation_set_size=atof(str2);
		if(strcmp(str1, "repulsor_min_age") == 0)
			config->repulsor_min_age=atoi(str2);
		if(strcmp(str1, "semantic_repulsor_max_number") == 0)
			config->semantic_repulsor_max_number=atoi(str2);
		if(strcmp(str1, "validation_elite_size") == 0)
			config->validation_elite_size=atoi(str2);
		if(strcmp(str1, "use_best_as_rep_candidate") == 0)
			config->use_best_as_rep_candidate=atoi(str2);
		if(strcmp(str1, "overfit_by_median") == 0)
			config->overfit_by_median=atoi(str2);
		if(strcmp(str1, "shuffle_validation_split") == 0)
			config->shuffle_validation_split=atoi(str2);
		if(strcmp(str1, "aggregate_repulsors") == 0)
			config->aggregate_repulsors=atoi(str2);
		if(strcmp(str1, "force_avoid_repulsors") == 0)
			config->force_avoid_repulsors=atoi(str2);
		if(strcmp(str1, "equality_delta") == 0)
			config->equality_delta=atof(str2);
		if(strcmp(str1, "true_pareto_selection") == 0)
			config->true_pareto_selection=atoi(str2);
		if(strcmp(str1, "domination_exclude_fitness") == 0)
			config->domination_exclude_fitness=atoi(str2);
		if(strcmp(str1, "merge_repulsors") == 0)
			config->merge_repulsors=atoi(str2);
	}
	f.close();
	if(config->p_crossover<0 || config->p_mutation<0 || config->p_crossover+config->p_mutation>1){
		cout<<"ERROR: CROSSOVER RATE AND MUTATION RATE MUST BE GREATER THAN (OR EQUAL TO) 0 AND THEIR SUM SMALLER THAN (OR EQUAL TO) 1.";
		exit(-1);
	}
	if(config->validation_set_size<0 || config->validation_set_size>=1){
		cout<<"ERROR: VALIDATION SET SIZE CANNOT BE SMALLER THAN 0 OR EQUAL (OR LARGER) TO 1.";
		exit(-1);
	}
}


void create_T_F(){
	NUM_VARIABLE_SYMBOLS=nvar;
	NUM_FUNCTIONAL_SYMBOLS=4;
	symbols.push_back(new symbol(1,2,1,"+"));
	symbols.push_back(new symbol(1,2,2,"-"));
	symbols.push_back(new symbol(1,2,3,"*"));
	symbols.push_back(new symbol(1,2,4,"/"));
	for(int i=NUM_FUNCTIONAL_SYMBOLS;i<NUM_VARIABLE_SYMBOLS+NUM_FUNCTIONAL_SYMBOLS;i++){
		char str[50] = "x";
		char buf[50]="";
		sprintf(buf, "%d", i-NUM_FUNCTIONAL_SYMBOLS);
		strcat( str, buf);
		symbols.push_back(new symbol(0,0,i,str));
	}
	for(int i=NUM_VARIABLE_SYMBOLS+NUM_FUNCTIONAL_SYMBOLS;i<NUM_VARIABLE_SYMBOLS+NUM_FUNCTIONAL_SYMBOLS+NUM_CONSTANT_SYMBOLS;i++){
		int a=config.min_random_constant+frand()*(config.max_random_constant-config.min_random_constant);
		char buf [50]="";
		stringstream s;
		s << a;
		string f;
		s>>f;
		strcpy(buf,f.c_str());
		symbols.push_back(new symbol(0,0,i,buf));
		symbols[symbols.size()-1]->value=a;
	}
}


int choose_function(){
	int index;
	index=int(frand()*(NUM_FUNCTIONAL_SYMBOLS-1));
	return index;
}

int choose_terminal(){
	int index;
	if(NUM_CONSTANT_SYMBOLS==0){
		index=int(NUM_FUNCTIONAL_SYMBOLS+frand()*(NUM_VARIABLE_SYMBOLS-1));
	}
	else{
		if(frand()<0.7)
			index=int(NUM_FUNCTIONAL_SYMBOLS+frand()*(NUM_VARIABLE_SYMBOLS-1));
		else
			index=int(NUM_FUNCTIONAL_SYMBOLS+NUM_VARIABLE_SYMBOLS+frand()*(NUM_CONSTANT_SYMBOLS-1));
	}
	return index;
}


void create_grow_pop(population **p){
	int i=(*p)->num_ind;
	while(i<config.population_size){
		(*p)->individuals[i]=new node;
		create_grow_tree((node**)&((*p)->individuals[i]),0, NULL, config.max_depth_creation);
		i++;
	}
}

void create_full_pop(population** p){
	int i=(*p)->num_ind;
	while(i<config.population_size){
		(*p)->individuals[i]=new node;
		create_full_tree((node**)&((*p)->individuals[i]),0, NULL, config.max_depth_creation);
		i++;
	}
}

void create_ramped_pop(population **p){
	int sub_pop;
	int r;
	int i=(*p)->num_ind;
	int min_depth;
	if(config.zero_depth==0){
		sub_pop=(config.population_size-(*p)->num_ind)/config.max_depth_creation;
		r=(config.population_size-(*p)->num_ind)%config.max_depth_creation;
		min_depth=1;
	}
	else{
		sub_pop=(config.population_size-(*p)->num_ind)/(config.max_depth_creation+1);
		r=(config.population_size-(*p)->num_ind)%(config.max_depth_creation+1);
		min_depth=0;
	}
	int j=config.max_depth_creation;
	while(j>=min_depth){
		if(j<config.max_depth_creation){
			for(int k=0; k<(int)(ceil((double)sub_pop/2)); k++){
				(*p)->individuals[i]=new node;
				create_full_tree((node**)&((*p)->individuals[i]),0, NULL, j);
				i++;
				(*p)->num_ind++;
			}
			for(int k=0; k<(int)(floor((double)sub_pop/2)); k++){
				(*p)->individuals[i]=new node;
				create_grow_tree((node**)&((*p)->individuals[i]),0, NULL ,j);
				i++;
				(*p)->num_ind++;
			}
		}
		else{
			for(int k=0; k<(int)(ceil((double)(sub_pop+r)/2)); k++){
				(*p)->individuals[i]=new node;
				create_full_tree((node**)&((*p)->individuals[i]),0, NULL, j);
				i++;
				(*p)->num_ind++;
			}
			for(int k=0; k<(int)(floor((double)(sub_pop+r)/2)); k++){
				(*p)->individuals[i]=new node;
				create_grow_tree((node**)&((*p)->individuals[i]),0, NULL, j);
				i++;
				(*p)->num_ind++;
			}
		}
		j--;
	}
}


void create_population(population **p, int i){
	if(i==0)
		create_grow_pop((population **)&(*p));
	if(i==1)
		create_full_pop((population **)&(*p));
	if(i==2)
		create_ramped_pop(p);
}


void create_grow_tree(node **el, int depth, node *parent, int max_depth){
	if(depth==0 && config.zero_depth==0){
		(*el)->root=symbols[choose_function()];
		(*el)->parent=NULL;
		(*el)->children=new node* [(*el)->root->arity];
		for (int i=0; i<(*el)->root->arity; i++){
			(*el)->children[i]=new node;
			create_grow_tree(((node **)&((*el)->children[i])), depth+1, *el, max_depth);
		}
		return;
	}
	if(depth==max_depth){
		(*el)->root=symbols[choose_terminal()];
		(*el)->parent=parent;
		(*el)->children=NULL;
		return;
	}
	if((depth>0 && depth<max_depth) || (depth==0 && config.zero_depth==1)){
		if(frand()>0.5){
			(*el)->root=symbols[choose_function()];
			(*el)->parent=parent;
			(*el)->children=new node* [(*el)->root->arity];
			for (int i=0; i<(*el)->root->arity; i++){
				(*el)->children[i]=new node;
				create_grow_tree(((node **)&((*el)->children[i])), depth+1, *el, max_depth);
			}
		}
		else{
			(*el)->root=symbols[choose_terminal()];
			(*el)->parent=parent;
			(*el)->children=NULL;
			return;
		}
	}
}

void create_full_tree(node **el, int depth, node *parent, int max_depth){
	if(depth==0 && depth<max_depth){
		(*el)->root=symbols[choose_function()];
		(*el)->parent=NULL;
		(*el)->children=new node* [(*el)->root->arity];
		for (int i=0; i<(*el)->root->arity; i++){
			(*el)->children[i]=new node;
			create_full_tree(((node **)&((*el)->children[i])), depth+1, *el, max_depth);
		}
		return;
	}
	if(depth==max_depth){
		(*el)->root=symbols[choose_terminal()];
		(*el)->parent=parent;
		(*el)->children=NULL;
		return;
	}
	if(depth>0 && depth<max_depth){
		(*el)->root=symbols[choose_function()];
		(*el)->parent=parent;
		(*el)->children=new node* [(*el)->root->arity];
		for (int i=0; i<(*el)->root->arity; i++){
			(*el)->children[i]=new node;
			create_full_tree(((node **)&((*el)->children[i])), depth+1, *el, max_depth);
		}
	}
}


double protected_division(double num, double den){
	if(den==0)
		return 1;
	else
		return  (num/den);
}


double eval(node *tree){
	if(tree->root->type==1){
		if(strcmp(tree->root->name,"+")==0){
			return (eval(tree->children[0])+eval(tree->children[1]));
		}
		if(strcmp(tree->root->name,"-")==0){
			return (eval(tree->children[0])-eval(tree->children[1]));
		}
		if(strcmp(tree->root->name,"*")==0){
			return (eval(tree->children[0])*eval(tree->children[1]));
		}
		if(strcmp(tree->root->name,"/")==0){
			return protected_division(eval(tree->children[0]),eval(tree->children[1]));
		}
	}
	else{
		return (tree->root->value);
	}
	cout<<"ERROR: UNDEFINED SYMBOL"<<endl;
	exit(-1);
}


void evaluate(population **p){
	(*p)->fitness[0]=make_tuple(Myevaluate((*p)->individuals[0]), 0);
	(*p)->fitness_val[0]=make_tuple(Myevaluate_val((*p)->individuals[0]), 0);
	(*p)->fitness_test[0]=make_tuple(Myevaluate_test((*p)->individuals[0]), 0);
	(*p)->index_best=0;
	fit_.push_back((*p)->fitness[0]);
	fit_val.push_back((*p)->fitness_val[0]);
	fit_test.push_back((*p)->fitness_test[0]);
	for(int i=1; i<config.population_size; i++){
		(*p)->fitness[i]=make_tuple(Myevaluate((*p)->individuals[i]), 0);
		(*p)->fitness_val[i]=make_tuple(Myevaluate_val((*p)->individuals[i]), 0);
		(*p)->fitness_test[i]=make_tuple(Myevaluate_test((*p)->individuals[i]), 0);
		fit_.push_back((*p)->fitness[i]);
		fit_val.push_back((*p)->fitness_val[i]);
		fit_test.push_back((*p)->fitness_test[i]);
		
		if(nsga_II_better((*p)->fitness[i],(*p)->fitness[(*p)->index_best])){
			(*p)->index_best=i;
		}
	}
	nsga_II_sort(p);
}


double Myevaluate (node *el) {
	double d=0;
	vector <double> val;
	for(int i=0;i<nrow;i++){
		update_terminal_symbols(i);
		set[i].res=eval(el);
		val.push_back(set[i].res);
		d+=(set[i].res-set[i].y_value)*(set[i].res-set[i].y_value);
	}
	sem_train_cases.push_back(val);
	d=sqrt(d/nrow);
	return d;
}

double Myevaluate_val (node *el) {
	double d=0;
	vector <double> val;
	for(int i=nrow;i<nrow+nrow_val;i++){
		update_terminal_symbols(i);
		set[i].res=eval(el);
		val.push_back(set[i].res);
		d+=(set[i].res-set[i].y_value)*(set[i].res-set[i].y_value);
	}
	sem_val_cases.push_back(val);
	d=sqrt(d/nrow_val);
	return d;
}


double Myevaluate_test (node *el) {
	double d=0;
	vector <double> val;
	for(int i=nrow+nrow_val;i<nrow+nrow_val+nrow_test;i++){
		update_terminal_symbols(i);
		set[i].res=eval(el);
		val.push_back(set[i].res);
		d+=(set[i].res-set[i].y_value)*(set[i].res-set[i].y_value);
	}
	sem_test_cases.push_back(val);
	d=sqrt(d/nrow_test);
	return d;
}


void Myevaluate_random (node *el, vector <double> & sem){
	for(int i=0;i<nrow;i++){
		update_terminal_symbols(i);
		set[i].res=eval(el);
		sem.push_back(set[i].res);
	}
}

void Myevaluate_random_val (node *el, vector <double> & sem){
	for(int i=nrow;i<nrow+nrow_val;i++){
		update_terminal_symbols(i);
		set[i].res=eval(el);
		sem.push_back(set[i].res);
	}
}

void Myevaluate_random_test(node *el, vector <double> & sem) {
	for(int i=nrow+nrow_val;i<nrow+nrow_val+nrow_test;i++){
		update_terminal_symbols(i);
		set[i].res=eval(el);
		sem.push_back(set[i].res);
	}
}

void nsga_II_sort(population **p) {
	
	if (sem_repulsors.size() == 0){
		clog<<"\t"<<"No semantic repulsors collected - skipping nsga_II_sort()"<<endl;
		return;
	}
	
	vector <int> domination_front;
	int *domination_counts = new int [config.population_size] {};
	vector< vector<int> > dominated_individuals;
	
	perform_fast_non_domination_sort(p, &domination_front, (int**)&domination_counts, &dominated_individuals);
	
	int front = 1;
	int count = 0;
	while (domination_front.size()!=0){
		count += domination_front.size();
		clog<<"\t"<<"Number of Individuals in front "<<front<<": "<<domination_front.size()<<endl;
		
		vector <int> next_front;
		extract_next_front(front, &next_front, p, &domination_front, (int**)&domination_counts, &dominated_individuals);
		
		// update iteration data
		front++;
		domination_front = next_front;
	}
	delete[] domination_counts;
	cout<<"Total individuals processed for fronts "<<count<<endl;
}

void perform_fast_non_domination_sort(population **p, vector<int> *d_front, int **d_counts, vector< vector<int> > *d_individuals){
	clog<<"\t"<<"Calculating domination count and set of dominated individuals for each individual in p"<<endl;
	for(int i=0; i<config.population_size; i++){
		vector<int> d_inds;
		for(int j=0; j<config.population_size; j++){
			if (i == j) continue; // no need to compare to oneself
			// determine domination of i over j, or vice versa based on fitness and all repulsor distances
			bool iDominatesJ = dominates(i, j);
			bool jDominatesI = dominates(j, i);
			if (iDominatesJ && jDominatesI){
				cout<<"WARNING: two individuals can't dominate each other.";
			}
			// update data structures
			if (iDominatesJ) // add j to set of dominated solutions of i
				d_inds.push_back(j);
			else if (jDominatesI)    // increment count of times that i has been dominated
				(*d_counts)[i]++;
		}
		(*d_individuals).push_back(d_inds);
		if ((*d_counts)[i] == 0){
			(*d_front).push_back(i);
			// add rank to data
			fit_new[i]=make_tuple(get<0>(fit_new[i]),1);
			(*p)->fitness[i]=make_tuple(get<0>(fit_new[i]),1);
		}
	}
}

void extract_next_front(int cur_front, vector <int>* next_front, population **p, vector<int> *d_front, int **d_counts, vector< vector<int> > *d_individuals){
	for (int i = 0; i < (*d_front).size(); i++){
		// extract next front and update ranks
		for (int d = 0; d < (*d_individuals)[(*d_front)[i]].size(); d++){
			int q = (*d_individuals)[(*d_front)[i]][d];
			if (find((*next_front).begin(), (*next_front).end(), q) == (*next_front).end()) { // check if the index has been added to the next front
				fit_new[q]=make_tuple(get<0>(fit_new[q]),cur_front+1);
				(*p)->fitness[q]=make_tuple(get<0>(fit_new[q]),cur_front+1);
				(*d_counts)[q]--;
				if ((*d_counts)[q] == 0)
					(*next_front).push_back(q);
			}
		}
	}
}

bool dominates(int i, int j){
	// determine domination of i over j
	bool iDominatesJ = better(get<0>(fit_new[i]), get<0>(fit_new[j]));
	double avg_dist_i = 0;
	double avg_dist_j = 0;
	bool iIsRepulsor = false;
	bool jIsRepulsor = false;
	// cout<<"Fitness: "<<get<0>(fit_new[i])<<" vs "<<get<0>(fit_new[j])<<""<<endl;
	for (int r = 0; r < sem_repulsors.size(); r++){
		// check if distance == 0 -> i/j is an repulsor
		if (repulsor_distances_new[i][r] == 0){
			iIsRepulsor = true;
		}
		if (repulsor_distances_new[j][r] == 0){
			jIsRepulsor = true;
		}
		avg_dist_i += repulsor_distances_new[i][r];
		avg_dist_j += repulsor_distances_new[j][r];
		if (config.domination_exclude_fitness == 1){
			if (r == 0){ // reset iDominatesJ (was initialized based on fitness)
				iDominatesJ = repulsor_distances_new[i][r] > repulsor_distances_new[j][r];
			} else{
				iDominatesJ = (iDominatesJ && repulsor_distances_new[i][r] > repulsor_distances_new[j][r]);
			}
		} else{
			iDominatesJ = (iDominatesJ && repulsor_distances_new[i][r] > repulsor_distances_new[j][r]);
		}
	}
	// check whether to aggregate the repulsers and reset domination (check if average distance is larger -> then it dominates)
	avg_dist_i = avg_dist_i/sem_repulsors.size();
	avg_dist_j = avg_dist_j/sem_repulsors.size();
	// cout<<i<<";"<<j<<endl;
	// cout<<get<0>(fit_new[i])<<";"<<get<0>(fit_new[j])<<endl;
	// cout<<avg_dist_i<<";"<<avg_dist_j<<endl;
	if (config.aggregate_repulsors==1){
		// cout<<"Avg. Distance: "<<avg_dist_i<<" vs "<<avg_dist_j<<""<<endl;
		iDominatesJ = better(get<0>(fit_new[i]), get<0>(fit_new[j])) && (avg_dist_i > avg_dist_j);
	}
	// force domination if i is repulsor (only if j is not a repulsor either)
	if (iIsRepulsor){
		iDominatesJ = false;
	} else if (!iIsRepulsor && jIsRepulsor){
		iDominatesJ = true;
	} else if (iIsRepulsor && jIsRepulsor){
		iDominatesJ = false;
	}
	return iDominatesJ;
}


void update_terminal_symbols(int i){
	for(int j=0; j<NUM_VARIABLE_SYMBOLS; j++){
		symbols[j+NUM_FUNCTIONAL_SYMBOLS]->value=set[i].vars[j];
	}
}


void delete_individual(node * el){
	if(el==NULL)
		return;
	if(el->children!=NULL){
		for(int i=0; i<el->root->arity; i++){
			delete_individual(el->children[i]);
		}
	}
	delete el;
}


int tournament_selection(){
	int *index=NULL;
	index=new int [config.tournament_size];
	// clog<<"\t[Tournament with "<<sem_repulsors.size()<<" repulsers of ";
	for(int i=0;i<config.tournament_size;i++){
		index[i]=int(frand()*(config.population_size-1));
		// clog<<index[i]<<",";
	}
	// clog<<"]";
	if (config.true_pareto_selection==0)
	{
		fitness_data best=fit_[index[0]];
		double best_traditional=get<0>(fit_[index[0]]);
		int best_index=index[0];
		int best_index_traditional=index[0];
		for(int j=1;j<config.tournament_size;j++){
			fitness_data opponent=fit_[index[j]];
			if (better(get<0>(opponent), get<0>(best))){
				best_traditional=get<0>(opponent);
				best_index_traditional=index[j];
			}
			if (nsga_II_better(opponent, best) && !isEqualToAnyRepulsor(best_index)){
				best=opponent;
				best_index=index[j];
			}
		}
		if (best_index_traditional != best_index)
			clog<<"[TOURNAMENT] Idx NSGA II: "<<best_index<<" vs. Idx Traditional: "<<best_index_traditional<<" | Ranks "<<get<1>(best)<<" vs. "<<get<1>(fit_[best_index_traditional])<<endl;
		delete[] index;
		return best_index;
	} else {
		// get best rank and apply frand() > 0.5
		vector<int> best_indices;
		int best_rank=get<1>(fit_[index[0]]);
		best_indices.push_back(0);
		for (int j = 1; j < config.tournament_size; j++){
			if (get<1>(fit_[index[j]]) < best_rank){
				// clear best indices and update best rank
				best_indices.clear();
				best_rank=get<1>(fit_[index[j]]);
				best_indices.push_back(j);
			} else if (get<1>(fit_[index[j]]) == best_rank){
				best_indices.push_back(j);
			} else {/*nothing*/}
		}
		int idx = best_indices[floor(frand()*best_indices.size())];
		delete[] index;
		return idx;
	}
}


void reproduction(int i){
	int p1 = i;
	bool allow_variation=true;
	if (i==index_best){
		allow_variation=false;
		if (calculateSmallestRepulsorDistance(i) == 0)// don't use eitism if index_best is repulsor
			allow_variation=true;
	}
	if(allow_variation){
		p1=tournament_selection();
	}
	// train
	sem_train_cases_new.push_back(sem_train_cases[p1]);
	fit_new.push_back(fit_[p1]);
	if (sem_repulsors.size()>0)
		repulsor_distances_new.push_back(repulsor_distances[p1]);
	// validation
	sem_val_cases_new.push_back(sem_val_cases[p1]);
	fit_new_val.push_back(fit_val[p1]);
	// test
	sem_test_cases_new.push_back(sem_test_cases[p1]);
	fit_new_test.push_back(fit_test[p1]);
}


void geometric_semantic_crossover(int i){
	bool allow_variation=true;
	if (i==index_best){
		allow_variation=false;
		if (calculateSmallestRepulsorDistance(i) == 0)// don't use eitism if index_best is repulsor
			allow_variation=true;
	}
	if(allow_variation){
		int p1=tournament_selection();
		int p2=tournament_selection();
		
		node* RT=new node;
		create_grow_tree((node**)&(RT),0, NULL, config.max_depth_creation);
		
		vector <double> sem_RT;
		vector <double> sem_RT_val;
		vector <double> sem_RT_test;
		
		Myevaluate_random(RT, sem_RT);
		Myevaluate_random_val(RT, sem_RT_val);
		Myevaluate_random_test(RT, sem_RT_test);
		
		delete_individual(RT);
		
		vector <double> val;
		vector <double> val_val;
		vector <double> val_test;
		
		// train
		for(int j=0;j<nrow;j++){
			double sigmoid=1.0/(1+exp(-(sem_RT[j])));
			val.push_back(sem_train_cases[p1][j]*(sigmoid)+sem_train_cases[p2][j]*(1-sigmoid));
		}
		sem_train_cases_new.push_back(val);
		update_training_fitness(val,1);

		// validation
		for(int j=0;j<nrow_val;j++){
			double sigmoid_val=1.0/(1+exp(-(sem_RT_val[j])));
			val_val.push_back(sem_val_cases[p1][j]*(sigmoid_val)+sem_val_cases[p2][j]*(1-sigmoid_val));
		}
		sem_val_cases_new.push_back(val_val);
		update_validation_fitness(val_val,1);
		update_validation_elite(sem_val_cases_new[i],get<0>(fit_new_val[i]));
		
		// test
		for(int j=0;j<nrow_test;j++){
			double sigmoid_test=1.0/(1+exp(-(sem_RT_test[j])));
			val_test.push_back(sem_test_cases[p1][j]*(sigmoid_test)+sem_test_cases[p2][j]*(1-sigmoid_test));
		}
		sem_test_cases_new.push_back(val_test);
		update_test_fitness(val_test,1);
		
		// repulsers
		vector <double> sems = val;
		sems.insert(sems.end(), val_val.begin(), val_val.end());
		update_repulsor_distances(sems,1);
	}
	
	else{
		// train
		sem_train_cases_new.push_back(sem_train_cases[i]);
		fit_new.push_back(fit_[i]);
		if (sem_repulsors.size()>0)
			repulsor_distances_new.push_back(repulsor_distances[i]);
		// validation
		sem_val_cases_new.push_back(sem_val_cases[i]);
		fit_new_val.push_back(fit_val[i]);
		// test
		sem_test_cases_new.push_back(sem_test_cases[i]);
		fit_new_test.push_back(fit_test[i]);
	}
}

void geometric_semantic_mutation(int i){
	bool allow_variation=true;
	if (i==index_best){
		allow_variation=false;
		if (calculateSmallestRepulsorDistance(i) == 0)// don't use eitism if index_best is repulsor
			allow_variation=true;
	}
	if(allow_variation){
		node* RT=new node;
		create_grow_tree((node**)&(RT),0, NULL, config.max_depth_creation);
		node* RT_2=new node;
		create_grow_tree((node**)&(RT_2),0, NULL, config.max_depth_creation);
		
		vector <double> sem_RT1;
		vector <double> sem_RT1_val;
		vector <double> sem_RT1_test;
		vector <double> sem_RT2;
		vector <double> sem_RT2_val;
		vector <double> sem_RT2_test;
		
		Myevaluate_random(RT,sem_RT1);
		Myevaluate_random_val(RT,sem_RT1_val);
		Myevaluate_random_test(RT,sem_RT1_test);
		Myevaluate_random(RT_2,sem_RT2);
		Myevaluate_random_val(RT_2,sem_RT2_val);
		Myevaluate_random_test(RT_2,sem_RT2_test);
		delete_individual(RT);
		delete_individual(RT_2);
		
		// train
		for(int j=0;j<nrow;j++){
			double sigmoid_1=1.0/(1+exp(-(sem_RT1[j])));
			double sigmoid_2=1.0/(1+exp(-(sem_RT2[j])));
			sem_train_cases_new[i][j]=sem_train_cases_new[i][j]+config.mutation_step*(sigmoid_1-sigmoid_2);
		}
		
		update_training_fitness(sem_train_cases_new[i],0);
		
		// validation
		for(int j=0;j<nrow_val;j++){
			double sigmoid_val_1=1.0/(1+exp(-(sem_RT1_val[j])));
			double sigmoid_val_2=1.0/(1+exp(-(sem_RT2_val[j])));
			sem_val_cases_new[i][j]=sem_val_cases_new[i][j]+config.mutation_step*(sigmoid_val_1-sigmoid_val_2);
		}
		// QUESTION: Is it okay to check overfitting before? It won't interfere right? if the individual would be added to the lsit, it wouldn't be overfitting...
		update_validation_fitness(sem_val_cases_new[i],0);
		update_validation_elite(sem_val_cases_new[i],get<0>(fit_new_val[i]));
		
		// test
		for(int j=0;j<nrow_test;j++){
			double sigmoid_test_1=1.0/(1+exp(-(sem_RT1_test[j])));
			double sigmoid_test_2=1.0/(1+exp(-(sem_RT2_test[j])));
			sem_test_cases_new[i][j]=sem_test_cases_new[i][j]+config.mutation_step*(sigmoid_test_1-sigmoid_test_2);
		}
		update_test_fitness(sem_test_cases_new[i],0);

		// repulsors
		vector <double> sems = sem_train_cases_new[i];
		sems.insert(sems.end(), sem_val_cases_new[i].begin(), sem_val_cases_new[i].end());
		update_repulsor_distances(sems,1);
	}
	else{
		long new_index = sem_train_cases_new.size()-1;
		// train
		sem_train_cases_new[new_index] = sem_train_cases[i];
		fit_new[new_index] = fit_[i];
		if (sem_repulsors.size()>0)
			repulsor_distances_new[new_index] = repulsor_distances[i];
		// validation
		sem_val_cases_new[new_index] = sem_val_cases[i];
		fit_new_val[new_index] = fit_val[i];
		// test
		sem_test_cases_new[new_index] = sem_test_cases[i];
		fit_new_test[new_index] = fit_test[i];
	}
}



void update_training_fitness(vector <double> semantic_values, bool crossover){
	double d=0;
	for(int j=0;j<nrow;j++){
		d+=(semantic_values[j]-set[j].y_value)*(semantic_values[j]-set[j].y_value);
	}
	if(crossover==1)
		fit_new.push_back(make_tuple(sqrt(d/nrow), 0));
	else
		fit_new[fit_new.size()-1]=make_tuple(sqrt(d/nrow), 0);
}

void update_validation_fitness(vector <double> semantic_values, bool crossover){
	double d=0;
	for(int j=nrow;j<nrow+nrow_val;j++){
		d+=(semantic_values[j-nrow]-set[j].y_value)*(semantic_values[j-nrow]-set[j].y_value);
	}
	if(crossover==1)
		fit_new_val.push_back(make_tuple(sqrt(d/nrow_val), 0));
	else
		fit_new_val[fit_new_val.size()-1]=make_tuple(sqrt(d/nrow_val), 0);
	
	
	// check overfitting, and add the individual to the repulsors
	int idx = fit_new_val.size()-1;
	if (config.use_best_as_rep_candidate == 0 && is_overfitting(get<0>(fit_new_val[idx]))){
		// add_repulsor(sem_train_cases_new[idx], get<0>(fit_new_val[idx]));
		vector<double> sems = sem_train_cases_new[idx];
		sems.insert(sems.end(), sem_val_cases_new[idx].begin(), sem_val_cases_new[idx].end());
		add_repulsor(sems, get<0>(fit_new_val[idx]));
	}
}


void update_test_fitness(vector <double> semantic_values, bool crossover){
	double d=0;
	for(int j=nrow+nrow_val;j<nrow+nrow_val+nrow_test;j++){
		d+=(semantic_values[j-nrow-nrow_val]-set[j].y_value)*(semantic_values[j-nrow-nrow_val]-set[j].y_value);
	}
	if(crossover == 1)
		fit_new_test.push_back(make_tuple(sqrt(d/nrow_test), 0));
	else
		fit_new_test[fit_new_test.size()-1]=make_tuple(sqrt(d/nrow_test), 0);
}

void update_repulsor_distances(vector <double> semantic_values, bool crossover){
	/// calculate distances to all repulsors and update repulsor_distances_new
	if (sem_repulsors.size()==0)
		return;
	
	vector <double> rds;
	double d = 0;
	for (int r=0; r < sem_repulsors.size(); r++){
		d = 0;
		for (int v=0; v < sem_repulsors[r].size(); v++){
			d += (sem_repulsors[r][v] - semantic_values[v])*(sem_repulsors[r][v] - semantic_values[v]);
		}
		rds.push_back(sqrt(d / sem_repulsors[r].size()));
	}
	if(crossover == 1)
		repulsor_distances_new.push_back(rds);
	else
		repulsor_distances_new[repulsor_distances_new.size()-1] = rds;
}

void update_validation_elite(vector <double> semantic_values, double fitness){
	if (fit_val_elite.size() < config.validation_elite_size || better(fitness, fit_val_elite[val_elite_worse_idx])){
		// exchange worst with the current semantics and update the data fields
		if (sem_val_elite.size()<config.validation_elite_size){
			// not enough individuals in list yet, so push back
			sem_val_elite.push_back(semantic_values);
			fit_val_elite.push_back(fitness);
		} else {
			// replace
			fit_val_elite[val_elite_worse_idx] = fitness;
			sem_val_elite[val_elite_worse_idx] = semantic_values;
		}
		if (config.overfit_by_median == 1){
			// set median as average
			std::sort(sem_val_elite.begin(), sem_val_elite.end());
			std::sort(fit_val_elite.begin(), fit_val_elite.end());
			long index = fit_val_elite.size()/2;
			if (fit_val_elite.size()%2 == 0){
				val_elite_avg_fit = (fit_val_elite[index-1] + fit_val_elite[index])/2;
			} else {
				val_elite_avg_fit = fit_val_elite[index];
			}
		} else {
			// update average
			val_elite_avg_fit = 0;
			for (int i  = 1; i<sem_val_elite.size(); i++){
				val_elite_avg_fit += fit_val_elite[i];
			}
			val_elite_avg_fit = val_elite_avg_fit/sem_val_elite.size();
		}
		// update index worse
		val_elite_worse_idx = 0;
		for (int i = 1; i<sem_val_elite.size(); i++){
			if (better(fit_val_elite[val_elite_worse_idx], fit_val_elite[i])){
				val_elite_worse_idx = i;
			}
		}
	}
}

int best_individual(){
	fitness_data best=fit_[0];
	int best_index1=0;
	for(int i=0;i<fit_.size();i++){
		// if(nsga_II_better(fit_[i],best)){
		// 	best=fit_[i];
		// 	best_index1=i;
		// }
		if(better(get<0>(fit_[i]),get<0>(best))){
			best=fit_[i];
			best_index1=i;
		}
	}
	
	// check if best individual is overfitting or not
	if (config.use_best_as_rep_candidate == 1 && is_overfitting(get<0>(fit_val[best_index1]))){
		clog<<"\tBest individual has been found to be overfitting."<<endl;
		vector<double> sems = sem_train_cases[best_index1];
		sems.insert(sems.end(), sem_val_cases[best_index1].begin(), sem_val_cases[best_index1].end());
		add_repulsor(sems, get<0>(fit_val[best_index1]));
		// add_repulsor(sem_train_cases[best_index1], get<0>(fit_val[best_index1]));
	}
	else if (config.use_best_as_rep_candidate > 1) {
		// TODO add sorting and checking of all n individuals
		clog<<"\tTesting "<<config.use_best_as_rep_candidate<<" for overfitting."<<endl;
		vector<int> sortedIndices;
		for (int i = 0; i < fit_.size(); i++){
			bool added=false;
			for (int v = 0; v < sortedIndices.size(); v++){
				if (fit_[sortedIndices[v]] <= fit_[i]) continue;
				sortedIndices.insert( sortedIndices.begin() + v, i);
				added=true;
				break;
			}
			if (!added){
				sortedIndices.push_back(i);
			}
		}
		for (int i = 0; i < config.use_best_as_rep_candidate; i++){
			// test best n individuals for overfitting
			if (is_overfitting(get<0>(fit_val[sortedIndices[i]]))){
				vector<double> sems = sem_train_cases[sortedIndices[i]];
				sems.insert(sems.end(), sem_val_cases[sortedIndices[i]].begin(), sem_val_cases[sortedIndices[i]].end());
				add_repulsor(sems, get<0>(fit_val[sortedIndices[i]]));
			}
		}
		clog<<"\tBest "<<sem_repulsors_new.size()<<" found to overfit."<<endl;
	}
	return best_index1;
}

bool is_overfitting(double fit){
	// compare to average fitness of validation elite set
	if (sem_val_elite.size() < config.validation_elite_size)
		return false;
	return (better(val_elite_avg_fit, fit));
}

double get_overfitting_severity(double fit){
	clog<<"\t\tOverfit Margin "<<(val_elite_avg_fit-fit)<<endl;
	return (val_elite_avg_fit-fit);
}

void add_repulsor(vector <double> semantics, double validation_fitness){
	bool add = true;
	// check in repulsor list
	for (int r = 0; r < sem_repulsors.size(); r++){
		bool same = true; 
		for (int d = 0; d < sem_repulsors[r].size(); d++){
			if (semantics[d] != sem_repulsors[r][d]){
				same = false;
				break;
			}
		}
		if (same){
			add = false;
			break;
		}
	}
	// check in pending repulsor list (if not already found in repulsor list)
	if (add){
		for (int r = 0; r < sem_repulsors_new.size(); r++){
			bool same = true; 
			for (int d = 0; d < sem_repulsors_new[r].size(); d++){
				if (semantics[d] != sem_repulsors_new[r][d]){
					same = false;
					break;
				}
			}
			if (same){
				add = false;
				break;
			}
		}
	}
	if (add){
		sem_repulsors_new.push_back(semantics);
		overfit_severity_new.push_back(get_overfitting_severity(validation_fitness));
	}
	else{
		clog<<"\tRepulsor already recorded"<<endl;
	}
}

void update_tables(){
	// repulsors
	repulsor_distances.clear();
	repulsor_distances.assign(repulsor_distances_new.begin(), repulsor_distances_new.end());
	repulsor_distances_new.clear();
	// training set
	fit_.clear();
	fit_.assign(fit_new.begin(),fit_new.end());
	fit_new.clear();
	sem_train_cases.clear();
	sem_train_cases.assign(sem_train_cases_new.begin(),sem_train_cases_new.end());
	sem_train_cases_new.clear();
	// validation set
	fit_val.clear();
	fit_val.assign(fit_new_val.begin(),fit_new_val.end());
	fit_new_val.clear();
	sem_val_cases.clear();
	sem_val_cases.assign(sem_val_cases_new.begin(),sem_val_cases_new.end());
	sem_val_cases_new.clear();
	// test set
	fit_test.clear();
	fit_test.assign(fit_new_test.begin(),fit_new_test.end());
	fit_new_test.clear();
	sem_test_cases.clear();
	sem_test_cases.assign(sem_test_cases_new.begin(),sem_test_cases_new.end());
	sem_test_cases_new.clear();
}

int update_repulsors(int num_gen){
	// check if population is old enough, otherwise discard repulsors
	if (num_gen < config.repulsor_min_age){
		clog<<"\tDiscarding "<<sem_repulsors_new.size()<<" repulsors due to population not old enough (min age = "<<config.repulsor_min_age<<")"<<endl;
		sem_repulsors_new.clear();
		return 0;
	}

	// add to repulsor table, if overfitting worse than the best of the repulsers
	long N = 0;
	if (config.semantic_repulsor_max_number > -1 && sem_repulsors.size() > config.semantic_repulsor_max_number)
		N = sem_repulsors.size()-config.semantic_repulsor_max_number;
	
	int idx_best = 0;
	for (int b = 1; b < overfit_severity.size(); b++){
		if (overfit_severity[b] > overfit_severity[idx_best])
			idx_best=b;
	}
	bool changed = false;
	for (int r = 0; r < sem_repulsors_new.size(); r++){
		if (changed){ // search the least overfitting again
			for (int b = 1; b < overfit_severity.size(); b++){
				if (overfit_severity[b] > overfit_severity[idx_best])
					idx_best=b;
			}			
			changed = false;
		}	
		// replace if the repulsor is overfitting worse, add if list is not full yet, ignore if none
		if (sem_repulsors.size() < config.semantic_repulsor_max_number || config.semantic_repulsor_max_number < 0){
			sem_repulsors.push_back(sem_repulsors_new[r]);
			overfit_severity.push_back(overfit_severity_new[r]);
			repulsor_map.push_back(make_tuple(repulsor_map.size(), sem_repulsors.size()-1, num_gen));
			clog<<"\tAdded repulser"<<endl;
			changed = true;
		} else if (overfit_severity[idx_best] > overfit_severity_new[r]){
			// replace repulser with new one
			if (config.merge_repulsors == 0){
				sem_repulsors[idx_best] = sem_repulsors_new[r];
				overfit_severity[idx_best] = overfit_severity_new[r];
				repulsor_map.push_back(make_tuple(repulsor_map.size(), idx_best, num_gen));
				changed = true;
				clog<<"\tReplaced repulser "<<idx_best<<" with repulser "<<r<<"."<<endl;
			} else {
				// merge two most similar repulsers and add the other one
				int freed_index = merge_repulsors(2);
				sem_repulsors[freed_index] = sem_repulsors_new[r];
				overfit_severity[freed_index] = overfit_severity_new[r];
				repulsor_map.push_back(make_tuple(repulsor_map.size(), freed_index, num_gen));
				changed = true;
			}
		} else {
			clog<<"\tIgnored repulser"<<endl;
		}
	}
	
	// re-evaluate whole population only if there are new repulsors and update repulsor_distances
	if (sem_repulsors.size() > 0 &&  sem_repulsors_new.size() > 0){
		repulsor_distances.clear();
		for (int i = 0; i < config.population_size; i++){
			vector <double> rds;
			vector <double> inds_sems = sem_train_cases[i];
			inds_sems.insert(inds_sems.end(), sem_val_cases[i].begin(), sem_val_cases[i].end());
			double d = 0;
			for (int r=0; r < sem_repulsors.size(); r++){
				d = 0;
				for (int v=0; v < sem_repulsors[r].size(); v++){
					d += (sem_repulsors[r][v] - inds_sems[v])*(sem_repulsors[r][v] - inds_sems[v]);
				}
				rds.push_back(sqrt(d / inds_sems.size()));
			}
			repulsor_distances.push_back(rds);
		}
	}
	clog<<"\t"<<"Processsed "<<sem_repulsors_new.size()<<" new semantic repulsors because of being worse than the average validation elite, total: "<<sem_repulsors.size()<<endl;
	sem_repulsors_new.clear();
	overfit_severity_new.clear();

	return N;
}

int merge_repulsors(int count){
	if (count < 2 || sem_repulsors.size() < 2)
		return -1;

	int idx1 = -1;
	int idx2 = -1;
	int distance = -1;
	for (int r1 = 0; r1 < sem_repulsors.size()-1; r1++){
		for (int r2 = r1+1; r2 < sem_repulsors.size(); r2++){
			// get distance between r1 and r2
			double dist = 0;
			for (int s = 0; s < sem_repulsors[r1].size(); s++){
				dist += (sem_repulsors[r1][s] - sem_repulsors[r2][s])*(sem_repulsors[r1][s] - sem_repulsors[r2][s]);
			}
			dist = sqrt(dist / sem_repulsors[r1].size());
			if (distance < 0){ // first round
				distance = dist;
				idx1 = r1;
				idx2 = r2;
			} else if (distance > dist){ // current two are closer than so far closest found
				distance = dist;
				idx1 = r1;
				idx2 = r2;
			}
		}
	}
	// get means
	vector<double> avg_sems;
	for (int s = 0; s < sem_repulsors[idx1].size(); s++){
		avg_sems.push_back((sem_repulsors[idx1][s] + sem_repulsors[idx2][s]) / 2);
	}
	double avg_severity = (overfit_severity[idx1] + overfit_severity[idx2]) / 2;
	// replace data
	sem_repulsors[idx1] = avg_sems;
	overfit_severity[idx1] = avg_severity;
	clog<<"\tMerged repulser "<<idx1<<" with repulser "<<idx2<<"."<<endl;
	// recurse
	merge_repulsors(count-1);
	// return first freed index
	return idx2;
}


void read_input_data(char *train_file, char *test_file){
	fstream in(train_file,ios::in);
	if (!in.is_open()) {
		cout<<endl<<"ERROR: TRAINING FILE NOT FOUND." << endl;
		exit(-1);
	}
	fstream in_test(test_file,ios::in);
	if (!in_test.is_open()) {
		cout<<endl<<"ERROR: TEST FILE NOT FOUND." << endl;
		exit(-1);
	}
	char str[255];
	in >> str;
	nvar = atoi(str);
	in_test >> str;
	nvar_test = atoi(str);
	in >> str;
	nrow = atoi(str);
	in_test >> str;
	nrow_test = atoi(str);
	set = new Instance[nrow+nrow_test];
	for (int i=0;i<nrow;i++) {
		set[i].vars = new double[nvar];
		for (int j=0; j<nvar; j++) {
			in >> str;
			set[i].vars[j] = atof(str);
		}
		in >> str;
		set[i].y_value = atof(str);
	}
	in.close();
	for (int i=nrow;i<nrow+nrow_test;i++) {
		set[i].vars = new double[nvar];
		for (int j=0; j<nvar; j++) {
			in_test >> str;
			set[i].vars[j] = atof(str);
		}
		in_test >> str;
		set[i].y_value = atof(str);
	}
	in_test.close();
	
	clog<<"\tRead Data Files: "<<endl;
	clog<<"\t\t"<<train_file<<" (training data with "<<nrow<<" instances)"<<endl;
	clog<<"\t\t"<<test_file<<" (test data with "<<nrow_test<<" instances)"<<endl;

	if (config.shuffle_validation_split==1){
		clog<<"\tShuffling training data"<<endl;
		int curIdx = nrow;
		Instance tempV;
		int rndIdx = 0;
		while (curIdx != 0){
			rndIdx = rand() % curIdx;
			curIdx = curIdx - 1; 
			tempV = set[curIdx];
			set[curIdx] = set[rndIdx];
			set[rndIdx] = tempV;
		}
	}

	clog<<"\t"<<"Splitting training data - validation set proportion = "<<config.validation_set_size<<endl;
	nrow_val = floor(config.validation_set_size*nrow);
	nrow = nrow-nrow_val;
	clog<<"\t\tTraining Instances: "<<nrow<<endl;
	clog<<"\t\tValidation Instances: "<<nrow_val<<endl;

}

bool better (double f1, double f2){
	if(config.minimization_problem==1){
		if(f1<f2)
			return true;
		else
			return false;
	}
	else{
		if(f1>f2)
			return true;
		else
			return false;
	}
}

bool nsga_II_better (fitness_data i1, fitness_data i2){
	if (sem_repulsors.size()>0){
		if (get<1>(i1) < get<1>(i2)){
			return true;
		} else if (get<1>(i1) == get<1>(i2)){
			return better(get<0>(i1), get<0>(i2));
		}
		return false;
	} else {
		return better(get<0>(i1), get<0>(i2));
	}
}

void calculateMaxDistance (){
	// calculate population max distance
	double maxD = 0;
	double cMaxD = 0;
	for (int a = 0; a < sem_train_cases_new.size()-1; a++){
		for (int b = a+1; b < sem_train_cases_new.size(); b++){
			double d = 0;
			for (int s = 0; s < nrow; s++){
				d += (sem_train_cases_new[a][s]-sem_train_cases_new[b][s])*(sem_train_cases_new[a][s]-sem_train_cases_new[b][s]);
			}
			double d_c = d;
			d = sqrt(d / nrow);
			if (d > maxD){
				maxD = d;
			}
			// add validation distance to d_c
			for (int s = 0; s < nrow_val; s++){
				d_c += (sem_val_cases_new[a][s]-sem_val_cases_new[b][s])*(sem_val_cases_new[a][s]-sem_val_cases_new[b][s]);	
			}
			d_c = sqrt(d_c / (nrow+nrow_val));
			if (d_c > cMaxD)
				cMaxD = d_c;
		}
	}
	population_max_distance = maxD;
	population_combined_max_distance = cMaxD;
}

bool isEqualToAnyRepulsor (int i){
	// return true of distance of individual i to any repulsor is less than maximum distance * equality_delta
	for (int r = 0; r < sem_repulsors.size(); r++){
		double d = 0;
		for (int s=0; s < nrow; s++){
			d += (sem_repulsors[r][s] - sem_train_cases_new[i][s])*(sem_repulsors[r][s] - sem_train_cases_new[i][s]);
		}
		for (int s=0; s < nrow_val; s++){
			d += (sem_repulsors[r][s+nrow] - sem_val_cases_new[i][s])*(sem_repulsors[r][s+nrow] - sem_val_cases_new[i][s]);
		}
		d = sqrt(d / (nrow+nrow_val));
		if (d < population_combined_max_distance*config.equality_delta){
			return true;
		}
	}
	return false;
}

double calculateSmallestRepulsorDistance (int i){
	// return the smallest found distance to any repulsor
	double minD = population_combined_max_distance;
	for (int r = 0; r < sem_repulsors.size(); r++){
		double d = 0;
		for (int s=0; s < nrow; s++){
			d += (sem_repulsors[r][s] - sem_train_cases_new[i][s])*(sem_repulsors[r][s] - sem_train_cases_new[i][s]);
		}
		for (int s=0; s < nrow_val; s++){
			d += (sem_repulsors[r][s+nrow] - sem_val_cases_new[i][s])*(sem_repulsors[r][s+nrow] - sem_val_cases_new[i][s]);
		}
		d = sqrt(d / (nrow+nrow_val));
		if (d < minD){
			minD = d;
		}
	}
	return minD;
}

void log_semantics (ofstream *csem, int num_gen){
	// csem<<"gen\tidx\tisRep\tsemantics on training data"<<endl;
	vector < vector<double> > semantics = sem_train_cases_new;
	// no variation has been executed yet
	if (num_gen == 0){
		semantics = sem_train_cases;
	}
	// log semantics of individuals on training data
	for (int i = 0; i < semantics.size(); i++){
		(*csem)<<num_gen<<";"<<i<<";0";
		for (int s = 0; s < semantics[i].size(); s++){
			(*csem)<<";"<<semantics[i][s];
		}
		(*csem)<<endl;
	}
	// log repulsor semantics
	for (int r = 0; r < sem_repulsors.size(); r++){
		(*csem)<<num_gen<<";"<<r<<";1";
		for (int s = 0; s < sem_repulsors[r].size(); s++){
			(*csem)<<";"<<sem_repulsors[r][s];
		}
		(*csem)<<endl;
	}	
}