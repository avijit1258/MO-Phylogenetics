//
// File: JTT92.cpp
// Created by: Julien Dutheil
// Created on: Wed Jan 21 14:09:43 2004
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "JTT92.h"

//From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

JTT92::JTT92(const ProteicAlphabet* alpha) :
  AbstractParameterAliasable("JTT92."),
  AbstractSubstitutionModel(alpha, "JTT92."),
  AbstractReversibleSubstitutionModel(alpha, "JTT92."),
  freqSet_(0)
{
  #include "__JTT92ExchangeabilityCode"
  #include "__JTT92FrequenciesCode"
  freqSet_ = new FixedProteinFrequenciesSet(alpha, freq_);
  updateMatrices();  
}

JTT92::JTT92(const ProteicAlphabet* alpha, ProteinFrequenciesSet* freqSet, bool initFreqs) :
  AbstractParameterAliasable("JTT92+F."),
  AbstractSubstitutionModel(alpha, "JTT92+F."),
  AbstractReversibleSubstitutionModel(alpha, "JTT92+F."),
  freqSet_(freqSet)
{
  #include "__JTT92ExchangeabilityCode"
  #include "__JTT92FrequenciesCode"
  if (initFreqs) freqSet_->setFrequencies(freq_);
  else freq_ = freqSet_->getFrequencies();
  freqSet_->setNamespace("JTT92+F."+freqSet_->getNamespace());
  addParameters_(freqSet_->getParameters());
  updateMatrices();  
}

/******************************************************************************/

void JTT92::setFreqFromData(const SequenceContainer& data)
{
  std::map<int, double> freqs;
  SequenceContainerTools::getFrequencies(data, freqs);
  double t = 0;
  for (unsigned int i = 0; i < size_; i++) t += freqs[i];
  for (unsigned int i = 0; i < size_; i++) freq_[i] = freqs[i] / t;
  freqSet_->setFrequencies(freq_);
  //Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/

