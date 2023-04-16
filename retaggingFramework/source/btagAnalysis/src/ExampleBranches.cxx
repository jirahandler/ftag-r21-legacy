#include "../btagAnalysis/ExampleBranches.hh"
#include "../btagAnalysis/ExampleBranchBuffer.hh"

#include "xAODJet/Jet.h"
#include "AthContainers/exceptions.h"
#include "TTree.h"



//!-----------------------------------------------------------------------------------------------------------------------------!//
ExampleBranches::ExampleBranches():
  m_branches(new ExampleBranchBuffer)
{
  // instantiate all the vectors here ...
  m_branches->example_vector_vector_branch = new std::vector<std::vector<int> >();
  m_branches->example_vector_branch = new std::vector<float>();

}

//!-----------------------------------------------------------------------------------------------------------------------------!//
ExampleBranches::~ExampleBranches() {
  // delete all the vectors here ...
  delete m_branches->example_vector_vector_branch;
  delete m_branches->example_vector_branch;

  delete m_branches;
}

void ExampleBranches::set_tree(TTree& output_tree) const {

  output_tree.Branch( "example_vector_vector"       , &m_branches->example_vector_vector_branch );
  output_tree.Branch( "example_vector"       , &m_branches->example_vector_branch );


}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void ExampleBranches::fill(const xAOD::Jet& jet) {

  float Modified_jetpT  = ExampleFunction( jet.pt() );

  m_branches->example_vector_branch->push_back(Modified_jetpT);

  std::vector<int> example_vector;

  for (int i = 0; i < 5; ++i)
  {
    example_vector.push_back(i);
  }

  m_branches->example_vector_vector_branch->push_back(example_vector);

}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void ExampleBranches::clear() {
  // clear vectors

  m_branches->example_vector_vector_branch->clear();
  m_branches->example_vector_branch->clear();


}


float ExampleBranches::ExampleFunction(float number){

  return number*3.0;

}






