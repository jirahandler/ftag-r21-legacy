/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "JetTagTools/JetTagUtils.h"

std::string JetTagUtils::getJetAuthor(xAOD::Jet& jetToTag) {

  xAOD::JetAlgorithmType::ID jetAlgID =  jetToTag.getAlgorithmType();
  //std::string name = xAOD::JetAlgorithmType::algName(jetAlgID);
  std::string name = xAOD::JetAlgorithmType::algName(jetAlgID);
  xAOD::JetInput::Type jetAlgType = jetToTag.getInputType();
  std::string type =  xAOD::JetInput::typeName(jetAlgType);
  std::string size = std::to_string(int(jetToTag.getSizeParameter()*10));

  // Special test for HI jet collections
  // (the attribute JetUnsubtractedScaleMomentum is specific to them)
  xAOD::JetFourMom_t v;
  if ( jetToTag.getAttribute<xAOD::JetFourMom_t>("JetUnsubtractedScaleMomentum",v) ) {
    type = "HI";
  }
    
  std::string author = name;
  author.append(size);
  author.append(type);

  return author;
}



