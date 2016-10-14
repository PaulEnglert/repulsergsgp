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
	/// variable that indicates whether each individual will be checked for overfitting or only the best of the population
	int use_only_best_as_rep_candidate;
	/// variable that indicates whether overfitting is determined by the median or the average of the validation elites fitness
	int overfit_by_median;
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

/// tuple containing a fitness measure, a non-dominated sorting rank and a crowded distance measure
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

/*!
 * \fn                 void read_config_file(cfg *config)
 * \brief             function that reads the configuration file
 * \param          cfg *config: pointer to the struct containing the variables needed to run the program
 * \return           void
 * \date             01/09/2012
 * \author          Mauro Castelli
 * \file               GP.h
 */
void read_config_file(cfg *config);


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
 * \fn               void update_repulsors()
 * \brief            function that updates the tables used to store the semantics of the repulsors, as well as recalculates the distances of the individuals to the repulsors
 * \return           void
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
void update_repulsors();


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
 * \fn               void create_fake_repulors(int num_repulsors)
 * \brief            function that randomly creates semantic repulsors for testing purposes
 * \param            int num_repulsors: number of repulsors to create
 * \return           void
 * \date             TODO add date
 * \author           Paul Englert
 * \file             GP.h
 */
void create_fake_repulsors(int num_repulsors);

void read_config_file(cfg *config){
	clog<<"\t"<<"Reading configuration data"<<endl;
	fstream f("configuration.ini", ios::in);
	if (!f.is_open()) {
		cerr<<"CONFIGURATION FILE NOT FOUND." << endl;
		exit(-1);
	}
	int k=0;
	while(!f.eof()){
		char str[100]="";
		char str2[100]="";
		int j=0;
		f.getline(str,100);
		if(str[0]!='\0'){
			while(str[j]!='='){
				j++;
			}
			j++;
			int i=0;
			while(str[j]==' '){
				j++;
			}
			while(str[j]!='\0'){
				str2[i] = str[j];
				j++;
				i++;
			}
		}
		clog<<"\t\t"<<str<<endl;
		if(k==0)
			config->population_size = atoi(str2);
		if(k==1)
			config->max_number_generations=atoi(str2);
		if(k==2)
			config->init_type=atoi(str2);
		if(k==3)
			config->p_crossover=atof(str2);
		if(k==4)
			config->p_mutation=atof(str2);
		if(k==5)
			config->max_depth_creation=atoi(str2);
		if(k==6)
			config->tournament_size=atoi(str2);
		if(k==7)
			config->zero_depth=atoi(str2);
		if(k==8)
			config->mutation_step=atof(str2);
		if(k==9){
			config->num_random_constants=atoi(str2);
			NUM_CONSTANT_SYMBOLS=config->num_random_constants;
		}
		if(k==10)
			config->min_random_constant=atof(str2);
		if(k==11)
			config->max_random_constant=atof(str2);
		if(k==12)
			config->minimization_problem=atoi(str2);
		if(k==13)
			config->validation_set_size=atof(str2);
		if(k==14)
			config->repulsor_min_age=atoi(str2);
		if(k==15)
			config->semantic_repulsor_max_number=atoi(str2);
		if(k==16)
			config->validation_elite_size=atoi(str2);
		if(k==17)
			config->use_only_best_as_rep_candidate=atoi(str2);
		if(k==18)
			config->overfit_by_median=atoi(str2);
		k++;
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
	while (domination_front.size()!=0){
		clog<<"\t"<<"Number of Individuals in front "<<front<<": "<<domination_front.size()<<endl;
		
		vector <int> next_front;
		extract_next_front(front, &next_front, p, &domination_front, (int**)&domination_counts, &dominated_individuals);
		
		// update iteration data
		front++;
		domination_front = next_front;
	}
	delete[] domination_counts;
}

void perform_fast_non_domination_sort(population **p, vector<int> *d_front, int **d_counts, vector< vector<int> > *d_individuals){
	clog<<"\t"<<"Calculating domination count and set of dominated individuals for each individual in p"<<endl;
	for(int i=0; i<config.population_size; i++){
		vector<int> d_inds;
		for(int j=0; j<config.population_size; j++){
			if (i == j) continue; // no need to compare to oneself
			// determine domination of i over j, or vice versa based on fitness and all repulsor distances
			bool iDominatesJ = better(get<0>(fit_new[i]), get<0>(fit_new[j]));
			bool jDominatesI = !iDominatesJ;
			for (int r = 0; r < sem_repulsors.size(); r++){
				iDominatesJ = (iDominatesJ && repulsor_distances_new[i][r] > repulsor_distances_new[j][r]);
				jDominatesI = (jDominatesI && repulsor_distances_new[i][r] < repulsor_distances_new[j][r]);
			}
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
	for(int i=0;i<config.tournament_size;i++){
		index[i]=int(frand()*(config.population_size-1));
	}
	fitness_data best=fit_[index[0]];
	int best_index=index[0];
	for(int j=1;j<config.tournament_size;j++){
		fitness_data opponent=fit_[index[j]];
		if (nsga_II_better(opponent, best)){
			best=opponent;
			best_index=index[j];
		}
	}
	delete[] index;
	return best_index;
}


void reproduction(int i){
	int p1 = i;
	if(i!=index_best){
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
	if(i!=index_best){
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
		update_repulsor_distances(val,1);
		
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
	if(i!=index_best){
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
		update_repulsor_distances(sem_train_cases_new[i],0);
		
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
	if (config.use_only_best_as_rep_candidate == 0 && is_overfitting(get<0>(fit_new_val[idx]))){
		sem_repulsors_new.push_back(sem_train_cases_new[idx]);
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
		for (int v=0; v < semantic_values.size(); v++){
			d += (sem_repulsors[r][v] - semantic_values[v])*(sem_repulsors[r][v] - semantic_values[v]);
		}
		rds.push_back(sqrt(d));
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
	int best_index=0;
	for(int i=0;i<fit_.size();i++){
		if(nsga_II_better(fit_[i],best)){
			best=fit_[i];
			best_index=i;
		}
	}
	
	// check if best individual is overfitting or not
	if (config.use_only_best_as_rep_candidate == 1 && is_overfitting(get<0>(fit_val[best_index]))){
		clog<<"\tBest individual has been found to be overfitting."<<endl;
		sem_repulsors_new.push_back(sem_train_cases[best_index]);
	}
	return best_index;
}

bool is_overfitting(double fit){
	// compare to average fitness of validation elite set
	if (sem_val_elite.size() < config.validation_elite_size)
		return false;
	return (better(val_elite_avg_fit, fit));
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

void update_repulsors(int num_gen){
	// check if population is old enough, otherwise discard repulsors
	if (num_gen < config.repulsor_min_age){
		clog<<"\tDiscarding "<<sem_repulsors_new.size()<<" repulsors due to population not old enough (min age = "<<config.repulsor_min_age<<")"<<endl;
		sem_repulsors_new.clear();
		return;
	}
	// add to repulsor table
	for (int r = 0; r < sem_repulsors_new.size(); r++){
		sem_repulsors.push_back(sem_repulsors_new[r]);
	}
	
	// enforce size constraint of configuration (-1 means no maximum limit is set)
	if (config.semantic_repulsor_max_number > -1 && sem_repulsors.size() > config.semantic_repulsor_max_number){
		long N = sem_repulsors.size()-config.semantic_repulsor_max_number;
		vector<decltype(sem_repulsors)::vector<double>>(sem_repulsors.begin()+N, sem_repulsors.end()).swap(sem_repulsors);
	}
	
	// re-evaluate whole population and update repulsor_distances
	repulsor_distances.clear();
	if (sem_repulsors.size() > 0){
		for (int i = 0; i < config.population_size; i++){
			vector <double> rds;
			double d = 0;
			for (int r=0; r < sem_repulsors.size(); r++){
				d = 0;
				for (int v=0; v < nrow; v++){
					d += (sem_repulsors[r][v] - sem_train_cases[i][v])*(sem_repulsors[r][v] - sem_train_cases[i][v]);
				}
				rds.push_back(sqrt(d));
			}
			repulsor_distances.push_back(rds);
		}
	}
	clog<<"\t"<<"Added "<<sem_repulsors_new.size()<<" semantic repulsors because of being worse than the average validation elite, total: "<<sem_repulsors.size()<<endl;
	sem_repulsors_new.clear();
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

void create_fake_repulsors(int num_repulsors){
	for (int i = 0; i < num_repulsors; i++){
		vector<double> sems;
		clog<<"\t"<<"Adding semantic repulsor with:";
		for (int v = 0; v < nvar; v++){
			sems.push_back(frand()*100*frand());
			clog<<" "<<sems[sems.size()-1];
		}
		clog<<endl;
		sem_repulsors.push_back(sems);
		// enforce size constraint
		if (sem_repulsors.size() > config.semantic_repulsor_max_number){
			long N = sem_repulsors.size()-config.semantic_repulsor_max_number;
			vector<decltype(sem_repulsors)::vector<double>>(sem_repulsors.begin()+N, sem_repulsors.end()).swap(sem_repulsors);
		}
	}
	
	// re-evaluate whole population and update repulsor_distances
	repulsor_distances.clear();
	for (int i = 0; i < config.population_size; i++){
		vector <double> rds;
		double d = 0;
		for (int r=0; r < sem_repulsors.size(); r++){
			d = 0;
			for (int v=0; v < nvar; v++){
				d += (sem_repulsors[r][v] - sem_train_cases[i][v])*(sem_repulsors[r][v] - sem_train_cases[i][v]);
			}
			rds.push_back(sqrt(d));
		}
		repulsor_distances.push_back(rds);
	}
}