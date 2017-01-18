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
/// variable that indicates wether the semantics should be written to a log file
    int log_semantics;

/// mutation variation configurations
    char m_func_t[4];	// type of function 1 of the mutation
    char m_func_mb[4];	// type of function 2 of the mutation
    int m_deg;			// degree of the linear combination
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
		/// array of training fitness values
		double *fitness;
		/// array of test fitness values
		double *fitness_test;
		/// class constructor
		population(){individuals=new node* [config.population_size]; num_ind=0;
		fitness=new double [config.population_size];
		fitness_test=new double [config.population_size];
		};
		/// class destructor
		~population() { delete[] individuals;}
};

/// array of training fitness values at generation g
vector <double> fit_;
/// array of test fitness values at generation g
vector <double> fit_test;
/// array of training fitness values at the current generation g+1
vector <double> fit_new;
/// array of test fitness values at the current generation g+1
vector <double> fit_new_test;

/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the training set at generation g
vector < vector<double> > sem_train_cases;
/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the training set at generation g+1
vector < vector<double> > sem_train_cases_new;
/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the test set at generation g
vector < vector<double> > sem_test_cases;
/// array where each element (that is also an array) contains the semantics of an individual of the population, computed on the test set at generation g+1
vector < vector<double> > sem_test_cases_new;

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
* \fn                void mutation_linear_combination(double t, double mb)
* \brief             function that calculates the linear combination for the mutation
* \param            double t: value of an individuals semantic
* \param            double mb: mutation base which will be added to the individuals semantic value
* \return           double
* \date             01/18/2017
* \author          Paul Englert
* \file             GP.h
*/
double mutation_linear_combination(double t, double mb);


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
		if(strcmp(str1, "m_func_1") == 0){
			for (int i = 0; i < 3; i++)
				config->m_func_t[i] = str2[i];
			config->m_func_t[3] = '\0';
		}
		if(strcmp(str1, "m_func_2") == 0){
			for (int i = 0; i < 3; i++)
				config->m_func_mb[i] = str2[i];
			config->m_func_mb[3] = '\0';
		}
		if(strcmp(str1, "m_deg") == 0)
			config->m_deg=atoi(str2);
	}	
    f.close();
    if(config->p_crossover<0 || config->p_mutation<0 || config->p_crossover+config->p_mutation>1){
        cout<<"ERROR: CROSSOVER RATE AND MUTATION RATE MUST BE GREATER THAN (OR EQUAL TO) 0 AND THEIR SUM SMALLER THAN (OR EQUAL TO) 1.";
        exit(-1);
    }
    static const vector<string> validValues {"pow", "log", "exp", "non"};
    if (std::find(validValues.begin(), validValues.end(), string(config->m_func_t)) == validValues.end()  ){
	    cout<<"ERROR: VALUE OF m_func_1 ("<<config->m_func_t<<") NOT KNOWN - PLEASE CHOOSE FROM {pow, exp, log, none} (none = no transformation).";
        exit(-1);
    }
    if (std::find(validValues.begin(), validValues.end(), string(config->m_func_mb)) == validValues.end()  ){
	    cout<<"ERROR: VALUE OF m_func_2 ("<<config->m_func_mb<<") NOT KNOWN - PLEASE CHOOSE FROM {pow, exp, log, none} (none = no transformation).";
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
		(*p)->fitness[0]=Myevaluate((*p)->individuals[0]);
		(*p)->index_best=0;
		fit_.push_back((*p)->fitness[0]);
		fit_test.push_back(Myevaluate_test((*p)->individuals[0]));
    	for(int i=1; i<config.population_size; i++){
    		(*p)->fitness[i]=Myevaluate((*p)->individuals[i]);
    		fit_.push_back((*p)->fitness[i]);
	       	fit_test.push_back(Myevaluate_test((*p)->individuals[i]));
            if(better((*p)->fitness[i],(*p)->fitness[(*p)->index_best])){
                (*p)->index_best=i;
            }
        }
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


double Myevaluate_test (node *el) {
	double d=0;
    vector <double> val;
    for(int i=nrow;i<nrow+nrow_test;i++){
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

void Myevaluate_random_test(node *el, vector <double> & sem) {
    for(int i=nrow;i<nrow+nrow_test;i++){
        update_terminal_symbols(i);
        set[i].res=eval(el);
        sem.push_back(set[i].res);
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
	double best_fitness=fit_[index[0]];
	int best_index=index[0];
	for(int j=1;j<config.tournament_size;j++){
		double fit=fit_[index[j]];
		if(better(fit,best_fitness)){
			best_fitness=fit;
			best_index=index[j];
		}
	}
	delete[] index;
	return best_index;
}


void reproduction(int i){
    if(i!=index_best){
        int p1=tournament_selection();
        sem_train_cases_new.push_back(sem_train_cases[p1]);
        fit_new.push_back(fit_[p1]);
        sem_test_cases_new.push_back(sem_test_cases[p1]);
        fit_new_test.push_back(fit_test[p1]);
    }
    else{
        sem_train_cases_new.push_back(sem_train_cases[i]);
        fit_new.push_back(fit_[i]);
        sem_test_cases_new.push_back(sem_test_cases[i]);
        fit_new_test.push_back(fit_test[i]);
    }
}


void geometric_semantic_crossover(int i){
    if(i!=index_best){
        int p1=tournament_selection();
        int p2=tournament_selection();

        node* RT=new node;
        create_grow_tree((node**)&(RT),0, NULL, config.max_depth_creation);
        
        vector <double> sem_RT;
        vector <double> sem_RT_test;

        Myevaluate_random(RT, sem_RT);
        Myevaluate_random_test(RT, sem_RT_test);
        
        delete_individual(RT);

        vector <double> val;
        vector <double> val_test;
        for(int j=0;j<nrow;j++){
            double sigmoid=1.0/(1+exp(-(sem_RT[j])));
	        val.push_back(sem_train_cases[p1][j]*(sigmoid)+sem_train_cases[p2][j]*(1-sigmoid));
        }
        sem_train_cases_new.push_back(val);
        update_training_fitness(val,1);

        for(int j=0;j<nrow_test;j++){
            double sigmoid_test=1.0/(1+exp(-(sem_RT_test[j])));
	        val_test.push_back(sem_test_cases[p1][j]*(sigmoid_test)+sem_test_cases[p2][j]*(1-sigmoid_test));
        }
            sem_test_cases_new.push_back(val_test);
            update_test_fitness(val_test,1);
        }

    else{
        sem_train_cases_new.push_back(sem_train_cases[i]);
        fit_new.push_back(fit_[i]);
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
        vector <double> sem_RT1_test;
        vector <double> sem_RT2;
        vector <double> sem_RT2_test;
        
        Myevaluate_random(RT,sem_RT1);
        Myevaluate_random_test(RT,sem_RT1_test);
        Myevaluate_random(RT_2,sem_RT2);
        Myevaluate_random_test(RT_2,sem_RT2_test);
        delete_individual(RT);
        delete_individual(RT_2);

        for(int j=0;j<nrow;j++){
            double sigmoid_1=1.0/(1+exp(-(sem_RT1[j])));
            double sigmoid_2=1.0/(1+exp(-(sem_RT2[j])));
            sem_train_cases_new[i][j]=mutation_linear_combination(sem_train_cases_new[i][j], config.mutation_step*(sigmoid_1-sigmoid_2));
        }

        update_training_fitness(sem_train_cases_new[i],0);

        for(int j=0;j<nrow_test;j++){
    	    double sigmoid_test_1=1.0/(1+exp(-(sem_RT1_test[j])));
    	    double sigmoid_test_2=1.0/(1+exp(-(sem_RT2_test[j])));
            sem_test_cases_new[i][j]=mutation_linear_combination(sem_test_cases_new[i][j], config.mutation_step*(sigmoid_test_1-sigmoid_test_2));
        }
        update_test_fitness(sem_test_cases_new[i],0);
    }
    else{
        sem_train_cases_new.push_back(sem_train_cases[i]);
        fit_new.push_back(fit_[i]);
        sem_test_cases_new.push_back(sem_test_cases[i]);
        fit_new_test.push_back(fit_test[i]);
    }
}

double mutation_linear_combination(double t, double mb){
	if (config.m_deg == 0){ // don't do anything
		return t+mb;
	}
	double result = ( frand() * 2 ) - 1; // rand -1 < x < 1
	for (int c = 0; c < config.m_deg; c++){
		double a = ( frand() * 2 ) - 1; 
		double b = ( frand() * 2 ) - 1;
		double t_n = t;
		// cout<<"Applying ";
		if (strcmp(config.m_func_t, "non") == 0){ // only multiply with random constant
			// cout<<"non";
			result += a*t_n;
		}
		if (strcmp(config.m_func_t, "pow") == 0){ // polynomal
			// cout<<"pow"<<c+1;
			result += a*pow(t_n, c+1);
		}
		if (strcmp(config.m_func_t, "exp") == 0){ // exponential
			// cout<<"exp";
			result += a*exp(t_n);
		}
		if (strcmp(config.m_func_t, "log") == 0){ // logarithmic
			// cout<<"log";
			result += a*log(t_n);
		}
		// cout<<" to "<<a<<"*"<<t_n<<" ";
		double mb_n = mb;
		// cout<<"Applying ";
		if (strcmp(config.m_func_mb, "non") == 0){
			// cout<<"non";
			result += b*mb_n;
		}
		if (strcmp(config.m_func_mb, "pow") == 0){
			// cout<<"pow"<<c+1;
			result += b*pow(mb_n, c+1);
		}
		if (strcmp(config.m_func_mb, "exp") == 0){
			// cout<<"exp";
			result += b*exp(mb_n);
		}
		if (strcmp(config.m_func_mb, "log") == 0){
			// cout<<"log";
			result += b*log(mb_n);
		}
		// cout<<" to "<<b<<"*"<<mb_n<<"\n";
	}
	return result;
}



void update_training_fitness(vector <double> semantic_values, bool crossover){
    double d=0;
    for(int j=0;j<nrow;j++){
        d+=(semantic_values[j]-set[j].y_value)*(semantic_values[j]-set[j].y_value);
    }
    if(crossover==1)
        fit_new.push_back(sqrt(d/nrow));
    else    
        fit_new[fit_new.size()-1]=sqrt(d/nrow);
}


void update_test_fitness(vector <double> semantic_values, bool crossover){
    double d=0;
    for(int j=nrow;j<nrow+nrow_test;j++){
        d+=(semantic_values[j-nrow]-set[j].y_value)*(semantic_values[j-nrow]-set[j].y_value);
    }
    if(crossover == 1)
        fit_new_test.push_back(sqrt(d/nrow_test));
    else    
        fit_new_test[fit_new_test.size()-1]=sqrt(d/nrow_test);
}


int best_individual(){
    double best_fitness=fit_[0];
    int best_index=0;
    for(unsigned int i=0;i<fit_.size();i++){
        if(better(fit_[i],best_fitness)){
            best_fitness=fit_[i];
            best_index=i;
        }
   }
   return best_index;
}


void update_tables(){
    fit_.clear();
   	fit_.assign(fit_new.begin(),fit_new.end());
   	fit_new.clear();
   	sem_train_cases.clear();
    sem_train_cases.assign(sem_train_cases_new.begin(),sem_train_cases_new.end());
    sem_train_cases_new.clear();
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

	clog<<"\tRead Data Files: "<<endl;
	clog<<"\t\t"<<train_file<<" (training data with "<<nrow<<" instances)"<<endl;
	clog<<"\t\t"<<test_file<<" (test data with "<<nrow_test<<" instances)"<<endl;
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


void log_semantics (ofstream *csem, int num_gen){
	vector < vector<double> > semantics = sem_train_cases_new;
	// no variation has been executed yet
	if (num_gen == 0){
		semantics = sem_train_cases;
	}
	// log semantics of individuals on training data
	for (int i = 0; i < semantics.size(); i++){
		(*csem)<<num_gen<<";"<<i;
		for (int s = 0; s < semantics[i].size(); s++){
			(*csem)<<";"<<semantics[i][s];
		}
		(*csem)<<endl;
	}
}