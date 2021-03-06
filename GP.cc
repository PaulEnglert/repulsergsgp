/*  Copyright � 2012 Mauro Castelli
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
	
	time_t start_time = time(nullptr);
	stringstream strm;
	strm << start_time;
	string stamp = strm.str();
	
	// redirect cout to GSGP.log
	ofstream clog("results/"+stamp+"-gsgp.log");
	auto old_rdbuf = std::clog.rdbuf();
	std::clog.rdbuf(clog.rdbuf());
	
    // name of the file with configuration data
    char path_config[150]="";
    // name of the file with training instances 
    char path_in[150]="";
    // name of the file with test instances
    char path_test[150]="";
    for (int i=1; i<argc-1; i++) {
        if(strncmp(argv[i],"-train_file",11) == 0) {
            strcat(path_in,argv[++i]);
        }
        else if(strncmp(argv[i],"-test_file",10) == 0) {
                strcat(path_test,argv[++i]);
        }        
        else if(strncmp(argv[i],"-config",7) == 0) {
                strcat(path_config,argv[++i]);
        } 
    }

	clog<<"Algorithm: Repulsor GSGP"<<endl;
	clog<<"Timestamp: "<<stamp<<endl<<endl;	
	/*
	 pointer to the file fitnesstrain.txt containing the training fitness of the best individual at each generation
	 */
	ofstream fitness_train("results/"+stamp+"-fitnesstrain.txt",ios::out);
	fitness_train<<"gen;idx;fitness;pr;#rep"<<endl;
	/*
	 pointer to the file fitnesstest.txt containing the validation fitness of the best individual at each generation
	 */
	ofstream fitness_val("results/"+stamp+"-fitnessvalidation.txt",ios::out);
    /*
     pointer to the file fitnesstest.txt containing the training fitness of the best individual at each generation
     */
    ofstream fitness_test("results/"+stamp+"-fitnesstest.txt",ios::out);
	/*
	 pointer to the file repulserdistances.txt containing the data on repulsor distances at each generation
	 */
	ofstream repulsor_log("results/"+stamp+"-repulserdistances.txt",ios::out);
    repulsor_log<<"gen;id;distance to best individual"<<endl;
	
	clog<<"Starting Setup Phase"<<endl;
	// initialization of the seed for the generation of random numbers
	srand((int)time (NULL));
	// reading the parameters of the GP algorithm
	read_config_file(&config, path_config);
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
	fitness_train<<"0;"<<index_best<<";"<<Myevaluate(p->individuals[p->index_best])<<endl;
	// writing the validation fitness of the best individual on the file fitnesstest.txt
	fitness_val<<"0;"<<Myevaluate_val(p->individuals[p->index_best])<<endl;
	// writing the test fitness of the best individual on the file fitnesstest.txt
	fitness_test<<"0;"<<Myevaluate_test(p->individuals[p->index_best])<<endl;
	// index of the best individual stored in the variable best_index
	index_best=best_individual();
	
	int reps_lost = 0;

	ofstream csem;
	if (config.log_semantics==1){
		csem.open("results/"+stamp+"-Semantics.txt");
		csem<<"gen;idx;isRep;semantics on training data"<<endl;
		log_semantics(&csem, 0);
	}

	clog<<"Finished Setup Phase"<<endl<<endl;
	
	// main GP cycle
	for(int num_gen=0; num_gen<config.max_number_generations; num_gen++){
		
		cout<<"Generation "<<num_gen+1<<endl;
		clog<<endl<<endl<<"\t GENERATION \t "<<num_gen+1<<endl<<endl;
		// clog<<"\tRerunning NSGA II Sort"<<endl;
		// nsga_II_sort((population**)&p);
		// creation of a new population (without building trees!!)
		clog<<"Starting Variation Phase"<<endl;
		int recreateCount = 0;
		bool best_updated = false;
		for(int k=0;k<config.population_size;k++){
			bool individual_accepted=true;
			bool first_repetition=true;
			do{
				individual_accepted = true;
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
				if (!first_repetition){
					recreateCount++;
				}
				if (config.force_avoid_repulsors==1 && isEqualToAnyRepulsor(sem_train_cases_new.size()-1)){
					individual_accepted=false;
					first_repetition=false;
					// remove k from data tables
					fit_new.erase(fit_new.begin()+k);
					fit_new_val.erase(fit_new_val.begin()+k);
					fit_new_test.erase(fit_new_test.begin()+k);
					sem_train_cases_new.erase(sem_train_cases_new.begin()+k);
					sem_val_cases_new.erase(sem_val_cases_new.begin()+k);
					sem_test_cases_new.erase(sem_test_cases_new.begin()+k);
					if (sem_repulsors.size() > 0)
						repulsor_distances_new.erase(repulsor_distances_new.begin()+k);
				}
			} while (recreateCount<50000 && !individual_accepted);
		}
		clog<<"\tRecreated "<<recreateCount<<" individuals."<<endl;
		if (config.force_avoid_repulsors == 1){
			calculateMaxDistance();
		}
		clog<<"Finished Variation Phase"<<endl<<endl;
		clog<<"Starting Non-Dominated Sorting Phase"<<endl;
		// update non-domination rank and crowded distance measure
		nsga_II_sort((population**)&p);
		clog<<"Finished Non-Dominated Sorting Phase"<<endl<<endl;
		clog<<"Starting Structure Update Phase"<<endl;

		// log semantics before anything gets updated for the next generation
		if (config.log_semantics==1 && (num_gen+1)%50==0){
			log_semantics(&csem, num_gen+1);
		}
		// updating the tables used to store semantics and fitness values
		update_tables();
		// index of the best individual stored in the variable best_index and overfitting check
		index_best=best_individual();
		// update the repulsors table and reevaluate the distances
		reps_lost = reps_lost + update_repulsors(num_gen);
		clog<<"Finished Updating of tables and updating repulsors"<<endl<<endl;
		
		
		clog<<"Outputting Generation Results"<<endl<<endl;
		// writing the  training fitness of the best individual on the file fitnesstrain.txt
		fitness_train<<num_gen+1<<";"<<index_best<<";"<<get<0>(fit_[index_best])<<";"<<get<1>(fit_[index_best])<<";"<<sem_repulsors.size()<<endl;
		// writing the repulsor distances to repulserdistances.txt
        for (int r = 0; r < sem_repulsors.size(); r++){
            // get last  known id
            int id = -1;
            for (int i = repulsor_map.size()-1; i >= 0; i--){
                if (get<1>(repulsor_map[i]) == r){
                    id = get<0>(repulsor_map[i]);
                    break;
                }
            }
            repulsor_log<<""<<num_gen+1<<";"<<id<<";"<<repulsor_distances[index_best][r]<<endl;
        }
		// writing the  validation fitness of the best individual on the file fitnesstest.txt
		fitness_val<<num_gen+1<<";"<<get<0>(fit_val[index_best])<<endl;
		// writing the  test fitness of the best individual on the file fitnesstest.txt
		fitness_test<<num_gen+1<<";"<<get<0>(fit_test[index_best])<<endl;

	}
	
	clog<<endl<<"Starting Cleanup"<<endl;
	// at the end of the execution all the data structures are deleted in order to deallocate memory
	for(int k=0; k<config.population_size; k++){
		delete_individual(p->individuals[k]);
	}
	delete[] p->fitness;
	delete[] p->fitness_val;
	delete[] p->fitness_test;
	delete p;
	for(int i=0; i<nrow+nrow_val+nrow_test; i++){
		delete[] set[i].vars;
	}
	delete[] set;
	for(long i=symbols.size()-1;i>=0;i--){
		delete symbols[i];
		symbols.erase(symbols.begin()+i);
	}
	symbols.clear();
	clog<<endl<<"Finished Cleanup"<<endl;
	
	// print out runtime
	clog<<endl<<endl<<"================================"<<endl;
	time_t end_time = time(nullptr);
	clog<<endl<<"Finished in "<<(end_time-start_time)<<"s"<<endl;
	
	if (config.log_semantics==1){
		csem.close();
	}

	// reset log buffer
	std::clog.rdbuf(old_rdbuf);
	return 0;
}