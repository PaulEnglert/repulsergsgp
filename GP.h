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
typedef tuple<double, int, double> fitness_data;
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
* \fn               void add_semantic_repulsor(vector <double> semantic_values)
* \brief            function that adds the given semantics to the repulsor structure
* \param            vector <double> semantic_values: semantics of the repulsor to add
* \return           void
* \date             TODO add date
* \author           Paul Englert
* \file             GP.h
*/
void add_semantic_repulsor(vector <double> semantic_values);

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
		return	(num/den);
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
		(*p)->fitness[0]=make_tuple(Myevaluate((*p)->individuals[0]), 0, 0);
		(*p)->fitness_val[0]=make_tuple(Myevaluate_val((*p)->individuals[0]), 0, 0);
		(*p)->fitness_test[0]=make_tuple(Myevaluate_test((*p)->individuals[0]), 0, 0);
		(*p)->index_best=0;
		fit_.push_back((*p)->fitness[0]);
		fit_val.push_back((*p)->fitness_val[0]);
		fit_test.push_back((*p)->fitness_test[0]);
    	for(int i=1; i<config.population_size; i++){
    		(*p)->fitness[i]=make_tuple(Myevaluate((*p)->individuals[i]), 0, 0);
    		(*p)->fitness_val[i]=make_tuple(Myevaluate_val((*p)->individuals[i]), 0, 0);
    		(*p)->fitness_test[i]=make_tuple(Myevaluate_test((*p)->individuals[i]), 0, 0);
    		fit_.push_back((*p)->fitness[i]);
	       	fit_val.push_back((*p)->fitness_val[i]);
	       	fit_test.push_back((*p)->fitness_test[i]);

            if(better(get<0>((*p)->fitness[i]),get<0>((*p)->fitness[(*p)->index_best]))){
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
		cout<<"No semantic repulsors collected yet - skipping nsga_II_sort()"<<endl;
		return;
	}

	cout<<"Updating nondominated rank and crowded distance measure"<<endl;

	vector <int> domination_front;

	// for all individuals in p, calculate domination count and set fo dominated solutions
	int *domination_counts = new int [config.population_size] {};
	vector< vector<int> > dominated_individuals;
	cout<<"Calculating domination count and set of dominated individuals for each individual in p"<<endl;
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
			else if (jDominatesI)	 // increment count of times that i has been dominated
				domination_counts[i]++;
		}	
		dominated_individuals.push_back(d_inds);
		if (domination_counts[i] == 0){
			domination_front.push_back(i);
			// add rank 1 to data
			fit_new[i]=make_tuple(get<0>(fit_new[i]),1,0);
			(*p)->fitness[i]=make_tuple(get<0>(fit_new[i]),1,0);
		}
	}
	int front = 1;
	while (domination_front.size()!=0){
		cout<<"Number of Individuals in front "<<front<<": "<<domination_front.size()<<endl;
		vector <int> next_front;
		for (int i = 0; i < domination_front.size(); i++){
			// extract next front and update ranks
			for (int d = 0; d < dominated_individuals[domination_front[i]].size(); d++){
				int q = dominated_individuals[domination_front[i]][d];
				if (find(next_front.begin(), next_front.end(), q) == next_front.end()) { // check if the index has been worked already
					fit_new[q]=make_tuple(get<0>(fit_new[q]),front+1,0);
					(*p)->fitness[q]=make_tuple(get<0>(fit_new[q]),front+1,0);
					next_front.push_back(q);
				}
			}
		}
		cout<<"Calculating crowded distance for each individual in front "<<front<<endl;
		// estimate the average cuboid around an individual formed by the nearest neihbours
		// QUESTION: should this also include the fitness, or just the distances to a repulsor? I'd say yes...
		for (int o = 0; o <= sem_repulsors.size(); o++){ 
			// gather values from the objective and add into ascendingly sorted vector front_values
			vector< tuple<int, double> > front_values;
			for (int i = 0; i < domination_front.size(); i++){
				tuple<int, double> value;
				if (o == sem_repulsors.size()){
					// use fitness as objective
					value = make_tuple(domination_front[i], get<0>(fit_new[i]));
				} else {
					// use distance to repulsor objective o
					value = make_tuple(domination_front[i], repulsor_distances_new[i][o]);
				}
				auto it = lower_bound(front_values.begin(), front_values.end(), value,
                          [](tuple<int, double> const &t1, tuple<int, double> const &t2)
                          { return get<1>(t1) < get<1>(t2); });
				front_values.insert(it, value);
			}

			// print out vector to check sorting
			// cout<<"sorted vector of values of front "<<front<<" corresponding to objective "<<o<<" (max "<<sem_repulsors.size()<<"):"<<endl;
			// for (int i = 0; i < front_values.size(); i++){
			// 	cout<<"("<<get<0>(front_values[i])<<","<<get<1>(front_values[i])<<"),";
			// }
			// cout<<endl;

			// set distance of first and last element to infinite
			int q = get<0>(front_values[0]);
			fit_new[q]=make_tuple(get<0>(fit_new[q]),get<1>(fit_new[q]),numeric_limits<double>::infinity());
			(*p)->fitness[q]=make_tuple(get<0>(fit_new[q]),get<1>(fit_new[q]),numeric_limits<double>::infinity());
			q = get<0>(front_values[front_values.size()-1]);
			fit_new[q]=make_tuple(get<0>(fit_new[q]),get<1>(fit_new[q]),numeric_limits<double>::infinity());
			(*p)->fitness[q]=make_tuple(get<0>(fit_new[q]),get<1>(fit_new[q]),numeric_limits<double>::infinity());
			for (int q = 1; q < front_values.size()-1; q++){
				// cout<<"calculating cd: "<<get<2>(fit_new[q])<<"+"<<get<1>(front_values[q+1])<<"-"<<get<1>(front_values[q-1])<<"/"<<get<1>(front_values[front_values.size()-1])<<"-"<<get<1>(front_values[0])<<endl;
				double dist = get<2>(fit_new[q]) + (get<1>(front_values[q+1])-get<1>(front_values[q-1]))/(get<1>(front_values[front_values.size()-1])-get<1>(front_values[0]));
				fit_new[q]=make_tuple(get<0>(fit_new[q]),get<1>(fit_new[q]),dist);
				(*p)->fitness[q]=make_tuple(get<0>(fit_new[q]),get<1>(fit_new[q]),dist);
			}

		}

		// update iteration data
		front++;
		domination_front = next_front;
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
		fitness_data data=fit_[index[j]];
		if (get<1>(data) < get<1>(best)){
			// data has a better pareto rank
			best=data;
			best_index=index[j];
		} else if (get<1>(data) == get<1>(best)){
			if (get<2>(data) > get<2>(best)){
				// data has the same pareto rank, but lies in a less crowded region
				best=data;
				best_index=index[j];				
			}
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
        update_validation_fitness(sem_val_cases_new[i],0);

        // test
        for(int j=0;j<nrow_test;j++){
    	    double sigmoid_test_1=1.0/(1+exp(-(sem_RT1_test[j])));
    	    double sigmoid_test_2=1.0/(1+exp(-(sem_RT2_test[j])));
            sem_test_cases_new[i][j]=sem_test_cases_new[i][j]+config.mutation_step*(sigmoid_test_1-sigmoid_test_2);
        }
        update_test_fitness(sem_test_cases_new[i],0);
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



void update_training_fitness(vector <double> semantic_values, bool crossover){
    double d=0;
    for(int j=0;j<nrow;j++){
        d+=(semantic_values[j]-set[j].y_value)*(semantic_values[j]-set[j].y_value);
    }
    if(crossover==1)
        fit_new.push_back(make_tuple(sqrt(d/nrow), 0, 0));
    else    
        fit_new[fit_new.size()-1]=make_tuple(sqrt(d/nrow), 0, 0);
}

void update_validation_fitness(vector <double> semantic_values, bool crossover){
    double d=0;
    for(int j=nrow;j<nrow+nrow_val;j++){
        d+=(semantic_values[j-nrow]-set[j].y_value)*(semantic_values[j-nrow]-set[j].y_value);
    }
    if(crossover==1)
        fit_new_val.push_back(make_tuple(sqrt(d/nrow_val), 0, 0));
    else    
        fit_new_val[fit_new_val.size()-1]=make_tuple(sqrt(d/nrow_val), 0, 0);
}


void update_test_fitness(vector <double> semantic_values, bool crossover){
    double d=0;
    for(int j=nrow+nrow_val;j<nrow+nrow_val+nrow_test;j++){
        d+=(semantic_values[j-nrow-nrow_val]-set[j].y_value)*(semantic_values[j-nrow-nrow_val]-set[j].y_value);
    }
    if(crossover == 1)
        fit_new_test.push_back(make_tuple(sqrt(d/nrow_test), 0, 0));
    else    
        fit_new_test[fit_new_test.size()-1]=make_tuple(sqrt(d/nrow_test), 0, 0);
}

void update_repulsor_distances(vector <double> semantic_values, bool crossover){
    /// calculate distances to all repulsors and update repulsor_distances_new
    if (sem_repulsors.size()==0) 
    	return;
    // cout<<"updating distances to semantic repulsors"<<endl;
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

int best_individual(){
    double best_fitness=get<0>(fit_[0]);
    int best_index=0;
    for(unsigned int i=0;i<fit_.size();i++){
        if(better(get<0>(fit_[i]),best_fitness)){
            best_fitness=get<0>(fit_[i]);
            best_index=i;
        }
   }
   return best_index;
}


void update_tables(){
	// training set
    fit_.clear();
   	fit_.assign(fit_new.begin(),fit_new.end());
   	fit_new.clear();
   	sem_train_cases.clear();
    sem_train_cases.assign(sem_train_cases_new.begin(),sem_train_cases_new.end());
    sem_train_cases_new.clear();
    if (sem_repulsors.size()>0) {
    	cout<<"Updating Tables: "<<endl;
    	cout<<"repulsor_distances size: "<<repulsor_distances.size()<<endl;
    	cout<<"repulsor_distances_new size: "<<repulsor_distances_new.size()<<endl;
    	repulsor_distances.clear();
	    repulsor_distances.assign(repulsor_distances_new.begin(), repulsor_distances_new.end());
	    repulsor_distances_new.clear();
	}
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

	cout<<"Splitting train set; validation set proportion = "<<config.validation_set_size<<endl;
	nrow_val = floor(config.validation_set_size*nrow);
	nrow = nrow-nrow_val;
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

void add_semantic_repulsor(vector <double> semantic_values){
	sem_repulsors.push_back(semantic_values);
	// re-evaluate whole population and update repulsor_distances
	repulsor_distances.clear();
	for (int i = 0; i < config.population_size; i++){
		vector <double> rds;
	    double d = 0;
	    for (int r=0; r < sem_repulsors.size(); r++){
	    	d = 0;
	    	for (int v=0; v < semantic_values.size(); v++){
	    		d += (sem_repulsors[r][v] - sem_train_cases[i][v])*(sem_repulsors[r][v] - sem_train_cases[i][v]);
	    	}
	    	rds.push_back(sqrt(d));
	    }
	    repulsor_distances.push_back(rds);
	}
	cout<<"Updated repulsor_distances to: "<<repulsor_distances.size()<<endl;
}

void create_fake_repulsors(int num_repulsors){
	for (int i = 0; i < num_repulsors; i++){
		vector<double> sems;
		cout<<"Adding semantic repulsor with:";
		for (int v = 0; v < nvar; v++){
			sems.push_back(frand()*100*frand());
			cout<<" "<<sems[sems.size()-1];
		}
		cout<<endl;
		add_semantic_repulsor(sems);
	}
}