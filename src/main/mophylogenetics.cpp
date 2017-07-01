/** 
 * MO-Phylogenetics (version 1.0.0) a software tool for multi-objective 
 * phylogenetic inference. This software integrates features 
 * of the jMetalCpp, Bio++ and PLL frameworks.
 * 
 * Copyright (C) 2017 Cristian Zambrano-Vega, Antonio J. Nebro.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For any other enquiries send an Email to Cristian Zambrano
 * czambrano@uteq.edu.ec
 *
 * When publishing work that uses this software please cite us.
 * 
 * @file main.cpp
 */

#include <cstdlib>

#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

#include <Problem.h>
#include <Phylogeny.h>
#include <Solution.h>
#include <TreeCrossover.h>
#include <PhylogeneticMutation.h>
#include <BinaryTournament2.h>
#include <RandomSelection.h>


#include <NSGAII.h>
#include <ssNSGAII.h>
#include <SMSEMOA.h>
#include <MOEAD.h>
#include <paes.h>
#include <PhyloMOCHC.h>


using namespace std;
using namespace bpp;


void Message()
{
  (*ApplicationTools::message <<"* Multi-Objective Phylogenetic Inference Software, jMetalC++ & Bio++ & PLL*").endLine();
  (*ApplicationTools::message  <<"*                                                                        *").endLine();
  (*ApplicationTools::message <<"* Authors:                                             Last Modif. 13/04/16*").endLine();
  (*ApplicationTools::message <<"* Cristian Zambrano-Vega                                                  *").endLine();
  (*ApplicationTools::message <<"* Antonio J. Nebro                                                        *").endLine();
  (*ApplicationTools::message <<"***************************************************************************").endLine();
  cout << endl;
}



int main(int argc, char** argv) {

  clock_t t_ini, t_fin;
  cout.precision(5);
  cout << fixed;

  Problem   * problem   ; // The problem to solve
  Algorithm * algorithm ; // The algorithm to use
  Operator  * crossover ; // Crossover operator
  Operator  * mutation  ; // Mutation operator
  Operator  * selection ; // Selection operator
  
  srand (time(NULL));
   
  Message();
  
  BppApplication * objApp = new BppApplication(argc, argv, "Params");

  string NumExp =  ApplicationTools::getStringParameter("experimentid", objApp->getParams(), "1", "", false, false);

  //Problem
  problem =  new Phylogeny(objApp);
  
 //MOEA  
  string AlgorithmName =   ApplicationTools::getStringParameter("algorithm", objApp->getParams(), "NSGAII", "", false, false);
  int populationSize,maxEvaluations, biSections,archiveSize,IntervalOptSubsModel;
  double offset;
  string dataDirectory;
  double initialConvergenceCount;
  double preservedPopulation;
  int convergenceValue;
  
   populationSize= ApplicationTools::getIntParameter("populationsize", objApp->getParams(), 100, "", false, false);
   maxEvaluations = ApplicationTools::getIntParameter("maxevaluations", objApp->getParams(), 2000, "", false, false);
   IntervalOptSubsModel= ApplicationTools::getIntParameter("intervalupdateparameters", objApp->getParams(), 500, "", false, false);
  
  if(AlgorithmName=="NSGAII"){
        algorithm = new NSGAII(problem);
  
  }else if(AlgorithmName=="ssNSGAII"){
        algorithm = new ssNSGAII(problem);
  
  }else if(AlgorithmName=="SMSEMOA"){
        offset = ApplicationTools::getDoubleParameter("offset", objApp->getParams(), 100, "", false, false);
        algorithm = new SMSEMOA(problem);
        algorithm->setInputParameter("offset",&offset);
        
        
  }else if(AlgorithmName=="MOEAD"){
        dataDirectory = ApplicationTools::getStringParameter("datadirectory", objApp->getParams(), "", "", false, false);
        algorithm = new MOEAD(problem);
        algorithm->setInputParameter("dataDirectory",&dataDirectory);
      
  }else if(AlgorithmName=="PAES"){
        biSections = ApplicationTools::getIntParameter("bisections", objApp->getParams(), 5, "", false, false);
        archiveSize = ApplicationTools::getIntParameter("archivesize", objApp->getParams(), 100, "", false, false);
        algorithm = new paes(problem);
        algorithm->setInputParameter("biSections",&biSections);
        algorithm->setInputParameter("archiveSize",&archiveSize);
      
  }else if(AlgorithmName=="MOCHC"){

      initialConvergenceCount = ApplicationTools::getDoubleParameter("initialconvergencecount", objApp->getParams(), 0.25, "", false, false);
      preservedPopulation = ApplicationTools::getDoubleParameter("preservedPopulation", objApp->getParams(), 0.05, "", false, false);
      convergenceValue = ApplicationTools::getIntParameter("convergencevalue", objApp->getParams(), 3, "", false, false);
      
      algorithm = new PhyloMOCHC(problem);
        
      algorithm->setInputParameter("initialConvergenceCount",&initialConvergenceCount);
      algorithm->setInputParameter("preservedPopulation",&preservedPopulation);
      algorithm->setInputParameter("convergenceValue",&convergenceValue);


  }else {
        algorithm = new NSGAII(problem);
  }
 
  algorithm->setInputParameter("populationSize",&populationSize);
  algorithm->setInputParameter("maxEvaluations",&maxEvaluations);
  algorithm->setInputParameter("intervalupdateparameters",&IntervalOptSubsModel);
  
        
  //Operator CrossOver
  map<string, void *> parameters;
  double crossoverProbability =  ApplicationTools::getDoubleParameter("crossover.probability", objApp->getParams(), 0.8, "", false, false); //0.8;
  int NumDescendientes =   ApplicationTools::getIntParameter("crossover.offspringsize", objApp->getParams(), 2, "", false, false);;
  parameters["probability"] =  &crossoverProbability;
  parameters["numDescendientes"] =  &NumDescendientes;
  crossover = new TreeCrossover(parameters);


  //Operator Mutation
  parameters.clear();
  double mutationProbability = ApplicationTools::getDoubleParameter("mutation.probability", objApp->getParams(), 0.2, "", false, false);
  double mutationDistributionIndex = 20;
  string OperadorMutacion= ApplicationTools::getStringParameter("mutation.method", objApp->getParams(), "NNI", "", false, false);
  
  parameters["probability"] = &mutationProbability;
  parameters["distributionIndex"] = &mutationDistributionIndex;
  parameters["metodo"] = &OperadorMutacion;
  mutation = new PhylogeneticMutation(parameters);
  
  
  //Selection Operator
  string OperadorSeleccion = ApplicationTools::getStringParameter("selection.method", objApp->getParams(), "binarytournament", "", false, false);
  
  parameters.clear();
  if(OperadorSeleccion=="binarytournament"){
      selection = new BinaryTournament2(parameters);
  }else if(OperadorSeleccion=="randomselection"){
        selection = new RandomSelection(parameters);     
  }else
      selection = new BinaryTournament2(parameters);
      
  algorithm->addOperator("crossover",crossover);
  algorithm->addOperator("mutation",mutation);
  algorithm->addOperator("selection",selection);
    
  
  cout <<"********************PARAMETERS**********************************" << endl;
  cout << "NumExp: " << NumExp << endl << endl;
  
  cout <<"**************************MOEA**********************************" << endl ;
  cout << "AlgorithmName: " << AlgorithmName << endl;
  cout << "populationSize: " << populationSize << endl;
  cout << "IntervalUpdateParameters: " << IntervalOptSubsModel << endl;
  
  if(AlgorithmName=="SMSEMOA"){
        cout << "offset: " << offset << endl;      
  }else if(AlgorithmName=="MOEAD"){
        cout << "dataDirectory: " << dataDirectory << endl;      
  }else if(AlgorithmName=="PAES"){
      cout << "biSections: " << biSections << endl;      
      cout << "archiveSize: " << archiveSize << endl;      
  } else if(AlgorithmName=="MOCHC"){
      cout << "initialConvergenceCount: " << initialConvergenceCount << endl;      
      cout << "preservedPopulation: " << preservedPopulation << endl;  
      cout << "convergenceValue: " << convergenceValue << endl;  
  }
  
  cout << "maxEvaluations: " << maxEvaluations << endl << endl;
  
  
  ((Phylogeny*)problem)->printParameters();
  ((TreeCrossover*)crossover)->printParameters();
  ((PhylogeneticMutation*)mutation)->printParameters();
  
  cout << "********************* Selection Operator ********************* " << endl;
  cout << "Method: " << OperadorSeleccion << endl;
  
  cout << endl ;
  
  
   cout << "****************** Start of Algorithm ***************" << endl;
    
   time_t timer = time(NULL);    printf("Start:  %s\n", ctime(&timer));
   t_ini = clock();

   SolutionSet * population = algorithm->execute();

   
   Phylogeny* p = (Phylogeny*) problem;
   if(p->FinalOptRamas){
        Solution *solution;
        double NewLik ;
        cout << "Optimizing Final Solutions - Population Size: " << population->size() << endl;
        
        for(int i=0; i<population->size(); i++){
                solution =  population->get(i);
                cout << "Optimizing Solution " << i + 1 << " Likelihood = " << solution->getObjective(1);
          
                Solution * Soltmp = new Solution(solution);
                
                NewLik=  p->BranchLengthOptimization(Soltmp,p->FinalMetodoOptRamas,p->FinalNumIterOptRamas,p->FinalTolerenciaOptRamas)* -1;
         
                if(NewLik < solution->getObjective(1)){
                        delete solution;
                        Soltmp->setObjective(1, NewLik);  
                        population->replace(i,Soltmp) ;
                        cout << " final: " << NewLik << endl;
                }else{
                        delete Soltmp;
                        cout << " final: -- NOT Improved -- "  << endl;
                }
        }
   
       timer = time(NULL);    printf("Final Optimization Ends %s\n", ctime(&timer));
   }
    
     /************************ SOLO OPBTENGO FRENTE NO DOMINADO *****************/
    cout << "Optimal Pareto front Approximation (only non dominated solutions)." << endl;
    Ranking *ranking = new Ranking(population);
    SolutionSet * FrenteOP = new SolutionSet(ranking->getSubfront(0)->size());
    for (int i=0;i<ranking->getSubfront(0)->size();i++) {
            FrenteOP->add(new Solution(ranking->getSubfront(0)->get(i)));
    }
    delete ranking;
    cout << "Pareto front Size: " <<  FrenteOP->size() << endl;
    
    cout << "Ordering Solutions" << endl;
    FrenteOP->sort( new DominanceComparator());
    cout << "Variables values have been written to file VAR" << endl;
    FrenteOP->printVariablesToFile("VAR" + NumExp);
    cout << "Objectives values have been written to file FUN" << endl;
    FrenteOP->printObjectivesToFile2("FUN" + NumExp);

    t_fin = clock(); 
    double secs = ((double) (t_fin - t_ini))/ CLOCKS_PER_SEC;
    cout << "Total execution time: " << secs << "s" << endl;
    
    timer = time(NULL);    printf("End of Algorithm %s\n", ctime(&timer));
    
    /*cout << "Plot Pareto Front" << endl;
    const char * PathFilePlot = "../resultados/frente"; 
    strcpy(PathFilePlot,NumExp.c_str());
    FrenteOP->plotObjectivesToFile(PathFilePlot, "FUN", 
                                 "Frente de soluciones NO Dominadas Problema 55 Taxas resuelto con NSGAII", 
                                 "/usr/bin/gnuplot");*/

    //((Phylogeny*)problem)->CloseScores();
    
  
  delete objApp;
  delete crossover;
  delete mutation;
  delete selection;
  //delete problem;
  delete algorithm;
  
    
}



