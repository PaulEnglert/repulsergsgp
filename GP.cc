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

//!  \file            GP.cc
//! \brief            file containing the main with the genetic programming algorithm
//! \date            created on 01/09/2012

#include "GP.h"

using namespace std;


/*!
* \fn              int main(int argc, const char **argv)
* \brief           main method that runs the GP algorithm
* \param           int argc: number of parameters of the program
* \param           const char **argv: array of strings that contains the parameters of the program
* \return          int: 0 if the program ends without errors
* \date            01/09/2012
* \author          Mauro Castelli
* \file            GP.cc
*/
int main(int argc, const char **argv){
    
    // name of the file with training instances 
    char path_in[50]="";
    // name of the file with test instances
    char path_test[50]="";
   	for (int i=1; i<argc-1; i++) {
        if(strncmp(argv[i],"-train_file",11) == 0) {
            strcat(path_in,argv[++i]);
        }
        else{
            if(strncmp(argv[i],"-test_file",10) == 0) {
                strcat(path_test,argv[++i]);
	      	}        
       	}
   	}

    /*
    pointer to the file fitnesstrain.txt containing the training fitness of the best individual at each generation
    */
    ofstream fitness_train("fitnesstrain.txt",ios::out);
    /*
    pointer to the file fitnesstest.txt containing the training fitness of the best individual at each generation
    */
    ofstream fitness_test("fitnesstest.txt",ios::out);
    
    // initialization of the seed for the generation of random numbers
	srand(time (NULL));
	// reading the parameters of the GP algorithm
	read_config_file(&config);
	// reading training and test files
	read_input_data(path_in,path_test);
	// creation of terminal and functional symbols
    create_T_F();
    // creation of an empty population
	population *p=new population();
	// initialization of the population
    create_population((population **)&p, config.init_type);	
    // evaluation of the individuals in the initial population
    evaluate((population**)&p);
    // writing the  training fitness of the best individual on the file fitnesstrain.txt
    fitness_train<<Myevaluate(p->individuals[p->index_best])<<endl;
	// writing the  test fitness of the best individual on the file fitnesstest.txt
    fitness_test<<Myevaluate_test(p->individuals[p->index_best])<<endl;
	// index of the best individual stored in the variable best_index
    index_best=best_individual();
	// main GP cycle
	for(int num_gen=0; num_gen<config.max_number_generations; num_gen++){	
        
        cout<<"Generation "<<num_gen+1<<endl;
        // creation of a new population (without building trees!!)
		for(int k=0;k<config.population_size;k++){
            double rand_num=frand();
            // geometric semantic crossover
            if(rand_num<config.p_crossover)
                geometric_semantic_crossover(k);
            // geometric semantic mutation    
            if(rand_num>=config.p_crossover && rand_num<config.p_crossover+config.p_mutation){
                reproduction(k);
	       		geometric_semantic_mutation(k);
            }
            // reproduction
            if(rand_num>=config.p_crossover+config.p_mutation){
                reproduction(k);
            }
        }
        
        // updating the tables used to store semantics and fitness values
		update_tables();
		// index of the best individual stored in the variable best_index
       	index_best=best_individual(); 
        // writing the  training fitness of the best individual on the file fitnesstrain.txt       
        fitness_train<<fit_[index_best]<<endl;
        // writing the  test fitness of the best individual on the file fitnesstest.txt
        fitness_test<<fit_test[index_best]<<endl;
    }    
    
    // at the end of the execution all the data structures are deleted in order to deallocate memory
	for(int k=0; k<config.population_size; k++){
        delete_individual(p->individuals[k]);
	}
	delete[] p->fitness;
	delete[] p->fitness_test;
	delete p;	
	for(int i=0; i<nrow+nrow_test; i++){
        delete[] set[i].vars;
	}
	delete[] set;
	for(int i=symbols.size()-1;i>=0;i--){
        delete symbols[i];
		symbols.erase(symbols.begin()+i);
	}
	symbols.clear();
	return 0;
}
