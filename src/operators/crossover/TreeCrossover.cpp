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
 * @file TreeCrossover.cpp
 */

#include <TreeCrossover.h>

/**
  * @class TreeCrossover
  * @brief This class is aimed at representing a TreeCrossover operator
**/


/**
 * Constructor
 * Create a new Tree crossover operator whit a default
 * index given by <code>DEFAULT_INDEX_CROSSOVER</code>
 */
TreeCrossover::TreeCrossover(map<string, void *> parameters)
: Crossover(parameters) {


  if (parameters["probability"] != NULL)
    crossoverProbability_ = *(double *)parameters["probability"];

  if (parameters["numDescendientes"] != NULL)
    numDescendientes_ = *(int *)parameters["numDescendientes"];
  else
    numDescendientes_ =2;

  

} // TreeCrossover


void TreeCrossover::printParameters() { 
  
  cout << "************** CrossOver Operator Parameters ************ " << endl;
  cout << "CrossOver: Prune-Delete-Graft" << endl;
  cout << "crossoverProbability: " << crossoverProbability_ << endl;
  cout << "Number of Offsprings:  " << numDescendientes_  << endl << endl;

}
/**
 * Destructor
 */
TreeCrossover::~TreeCrossover() { } // ~TreeCrossover


/**
* Perform the crossover operation.
* @param probability Crossover probability
* @param parent1 The first parent
* @param parent2 The second parent
* @return An array containing the two offsprings
**/
Solution * TreeCrossover::doCrossover(double probability, Solution *parent1, Solution *parent2) {

  Solution * offSpring;
  PhyloTree *copyPT1, *copyPT2, *offspringTree;
  Phylogeny * Problem = (Phylogeny*)parent1->getProblem();
  
  if(PseudoRandom::randDouble()<0.5)  {
      offSpring = new Solution(parent1);
      copyPT1 = (PhyloTree *)parent2->getDecisionVariables()[0];
      
  }else{
      offSpring = new Solution(parent2);
      copyPT1 = (PhyloTree *)parent1->getDecisionVariables()[0];
  }

  copyPT2 = (PhyloTree *)offSpring->getDecisionVariables()[0];
  copyPT2->setModificada(false);

  if (PseudoRandom::randDouble() < probability){

            TreeTemplate<Node> * tree;
            bool b=true;
            do{
                 offspringTree = new PhyloTree(copyPT2);  //Copy of PT2
               
                 //Cross PT1 y offspringTree, offspringTree is really affected
                 CrossTrees(copyPT1, offspringTree);
                 offspringTree->setModificada(true);
             
                 tree = offspringTree->getTree();
                 
                 if(!Problem->PLLisTreeValidate(tree)){
                     delete offspringTree;
                     //cout << "Invalited Crossed Tree " << endl ;
                 }else { 
                     b=false;
                 }
                 
             }while(b);
          
            delete copyPT2;
            offSpring->getDecisionVariables()[0]=offspringTree;
          
  }

  return offSpring;

} // doCrossover

/**
* Perform the crossover operation.
* @param PtMon  y PtDad deben ser COPIAS y PtDad resulta el ArbolCruzado
* @param parent2 The second parent
* @return One offspring
**/

/***************PTDad es afectado por el cruce entre PtMon y PtDad*****************/
void TreeCrossover::CrossTrees(PhyloTree * PtMon, PhyloTree * PtDad) {

	TreeTemplate<Node> * tree1 = PtMon->getTree();
	TreeTemplate<Node> * tree2 = PtDad->getTree();

        
	Node *Nodo1;
	Node *Nodo2;
        Node *SubTree;
        vector<string> hojas;
        vector<int> nodosIDs;

	int NumHojasT1 = tree1->getNumberOfLeaves();
        nodosIDs = tree1->getNodesId();

         bool b=true;
	 do{
		  Nodo1 = selectNodeToCross(tree1, nodosIDs );
		  if(TreeTemplateTools::getNumberOfLeaves(*Nodo1) < NumHojasT1-4){

                          SubTree=TreeTemplateTools::cloneSubtree<Node>(*Nodo1);
			  b=false;
   		  }
         }while(b);


       //Borra los nodos que se repiten en la Madre
	hojas = TreeTemplateTools::getLeavesNames(*SubTree);
        
//        cout << "SubTree to Cross" << endl;
//        for (unsigned int i = 0; i < hojas.size(); i++){
//            cout << hojas[i] << " " ;
//        }cout << endl;
        
	for (unsigned int i = 0; i < hojas.size(); i++)
            TreeTemplateTools::dropLeaf(*tree2,hojas[i]);


	Nodo2 = selectNodeToCross(tree2, tree2->getNodesId()); //Extrae PUNTERO refrencia a NODO, NO es copia


	Node * padre = Nodo2->getFather();
	double distancetofather = Nodo2->getDistanceToFather();
	int PosNodo2= padre->getSonPosition(Nodo2);


	Node * nodo = new Node();
	// agrega vector SON los punteros y establce hijos->father=nodo
	nodo->addSon(SubTree); nodo->addSon(Nodo2);
	Nodo2->setDistanceToFather(distancetofather/2);

	//Agrego Nuevo Nodo al Padre
	padre->setSon(PosNodo2, nodo);
	nodo->setDistanceToFather(distancetofather/2);



	tree2->resetNodesId(); //Uniion de otro subNodos se deben resetear ID

	 //ToDO :: RECALCULAR BRANCHS LENGTHS
        
        //string treenewick = TreeTemplateTools::treeToParenthesis(*tree2) ; cout << "Crossed Tree " << endl << treenewick << endl;

} // CrossTrees

/*Node * TreeCrossover::selectNodeToCross(TreeTemplate<Node> * tree_ ){
  	Node * nodoSel;
  	vector<int> nodosIDs = tree_->getNodesId();
  	do{
  		nodoSel =  tree_->getNode( nodosIDs[PseudoRandom::randInt(0, nodosIDs.size() - 1)]);
  	}while(!nodoSel ->hasFather() || nodoSel ->isLeaf() );

  	return  TreeTemplateTools::cloneSubtree<Node>(*tree_, nodoSel->getId());

  	//TreeTemplate<Node> *S = new TreeTemplate<Node>();
  	//S->setRootNode(nodo);
  	//cout << " Seleccion SubTree inside " << S->getNumberOfLeaves() << " Hojas " << endl;

  	//return nodo;
}*/

 Node * TreeCrossover::selectNodeToCross(TreeTemplate<Node> * tree_, vector<int> nodosIDs ){
   	Node * nodoSel;
   	do{
  		nodoSel =  tree_->getNode( nodosIDs[PseudoRandom::randInt(0, nodosIDs.size() - 1)]);
   	}while(!nodoSel->hasFather() || nodoSel->isLeaf() );

   	return  nodoSel;
 }

/**
* Executes the operation
* @param object An object containing an array of two parents
* @return An object containing the offSprings
*/
void * TreeCrossover::execute(void *object) {


  if(numDescendientes_==1)
  {
    void ** objects = (void **) object;
    Solution ** parents = (Solution **) objects[1];
    return  doCrossover(crossoverProbability_, parents[1], parents[2]);

  }else{

   Solution ** parents = (Solution **) object;
   Solution ** offSpring = new Solution*[2];
   
   offSpring[0] = doCrossover(crossoverProbability_, parents[0], parents[1]);
   offSpring[1] = doCrossover(crossoverProbability_, parents[1], parents[0]);
   
   return offSpring;
  }
   
} // execute


