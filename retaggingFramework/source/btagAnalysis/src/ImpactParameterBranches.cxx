#include "../btagAnalysis/ImpactParameterBranches.hh"
#include "../btagAnalysis/ImpactParameterBranchBuffer.hh"

#include "xAODJet/Jet.h"
#include "AthContainers/exceptions.h"
#include "TTree.h"

//!-----------------------------------------------------------------------------------------------------------------------------!//
ImpactParameterBranches::ImpactParameterBranches():
  m_branches(new ImpactParameterBranchBuffer)
{
  // instantiate all the vectors here ...
  m_branches->v_jet_ip2d_pb  = new std::vector<float>();
  m_branches->v_jet_ip2d_pc  = new std::vector<float>();
  m_branches->v_jet_ip2d_pu  = new std::vector<float>();
  m_branches->v_jet_ip2d_llr = new std::vector<float>();
  m_branches->v_jet_ip3d_pb  = new std::vector<float>();
  m_branches->v_jet_ip3d_pc  = new std::vector<float>();
  m_branches->v_jet_ip3d_pu  = new std::vector<float>();
  m_branches->v_jet_ip3d_llr = new std::vector<float>();

  m_branches->v_jet_ip2 = new std::vector<float>();
  m_branches->v_jet_ip2_c = new std::vector<float>();
  m_branches->v_jet_ip2_cu = new std::vector<float>();

  m_branches->v_jet_ip2_nan    = new std::vector<float>();
  m_branches->v_jet_ip2_c_nan  = new std::vector<float>();
  m_branches->v_jet_ip2_cu_nan = new std::vector<float>();

  m_branches->v_jet_ip3 = new std::vector<float>();
  m_branches->v_jet_ip3_c = new std::vector<float>();
  m_branches->v_jet_ip3_cu = new std::vector<float>();

  m_branches->v_jet_ip3_nan    = new std::vector<float>();
  m_branches->v_jet_ip3_c_nan  = new std::vector<float>();
  m_branches->v_jet_ip3_cu_nan = new std::vector<float>();

  m_branches->v_jet_rnnip_pb = new std::vector<float>();
  m_branches->v_jet_rnnip_pc = new std::vector<float>();
  m_branches->v_jet_rnnip_pu = new std::vector<float>();
  m_branches->v_jet_rnnip_ptau = new std::vector<float>();
}

//!-----------------------------------------------------------------------------------------------------------------------------!//
ImpactParameterBranches::~ImpactParameterBranches() {
  // delete all the vectors here ...
  delete m_branches->v_jet_ip2d_pb;
  delete m_branches->v_jet_ip2d_pc;
  delete m_branches->v_jet_ip2d_pu;
  delete m_branches->v_jet_ip2d_llr;
  delete m_branches->v_jet_ip3d_pb;
  delete m_branches->v_jet_ip3d_pc;
  delete m_branches->v_jet_ip3d_pu;
  delete m_branches->v_jet_ip3d_llr;

  delete m_branches->v_jet_ip2;
  delete m_branches->v_jet_ip2_c;
  delete m_branches->v_jet_ip2_cu;

  delete m_branches->v_jet_ip2_nan;
  delete m_branches->v_jet_ip2_c_nan;
  delete m_branches->v_jet_ip2_cu_nan;

  delete m_branches->v_jet_ip3;
  delete m_branches->v_jet_ip3_c;
  delete m_branches->v_jet_ip3_cu;

  delete m_branches->v_jet_ip3_nan;
  delete m_branches->v_jet_ip3_c_nan;
  delete m_branches->v_jet_ip3_cu_nan;


  delete m_branches->v_jet_rnnip_pb ;
  delete m_branches->v_jet_rnnip_pc ;
  delete m_branches->v_jet_rnnip_pu ;
  delete m_branches->v_jet_rnnip_ptau;

  delete m_branches;
}

void ImpactParameterBranches::set_tree(TTree& output_tree, std::map<std::string, double > defaultDict, bool replaceDefaults){

  m_defaultDict = defaultDict;
  m_replaceDefaults = replaceDefaults;

  output_tree.Branch( "jet_ip2d_pb" , &m_branches->v_jet_ip2d_pb );
  output_tree.Branch( "jet_ip2d_pc" , &m_branches->v_jet_ip2d_pc );
  output_tree.Branch( "jet_ip2d_pu" , &m_branches->v_jet_ip2d_pu );
  output_tree.Branch( "jet_ip2d_llr" , &m_branches->v_jet_ip2d_llr);
  output_tree.Branch( "jet_ip3d_pb" , &m_branches->v_jet_ip3d_pb );
  output_tree.Branch( "jet_ip3d_pc" , &m_branches->v_jet_ip3d_pc );
  output_tree.Branch( "jet_ip3d_pu" , &m_branches->v_jet_ip3d_pu );
  output_tree.Branch( "jet_ip3d_llr" , &m_branches->v_jet_ip3d_llr);

  output_tree.Branch( "jet_ip2", &m_branches->v_jet_ip2    );
  output_tree.Branch( "jet_ip2_c", &m_branches->v_jet_ip2_c  );
  output_tree.Branch( "jet_ip2_cu", &m_branches->v_jet_ip2_cu );

  output_tree.Branch( "jet_ip2_nan", &m_branches->v_jet_ip2_nan    );
  output_tree.Branch( "jet_ip2_c_nan", &m_branches->v_jet_ip2_c_nan  );
  output_tree.Branch( "jet_ip2_cu_nan", &m_branches->v_jet_ip2_cu_nan );

  output_tree.Branch( "jet_ip3", &m_branches->v_jet_ip3    );
  output_tree.Branch( "jet_ip3_c", &m_branches->v_jet_ip3_c  );
  output_tree.Branch( "jet_ip3_cu", &m_branches->v_jet_ip3_cu );

  output_tree.Branch( "jet_ip3_nan", &m_branches->v_jet_ip3_nan    );
  output_tree.Branch( "jet_ip3_c_nan", &m_branches->v_jet_ip3_c_nan  );
  output_tree.Branch( "jet_ip3_cu_nan", &m_branches->v_jet_ip3_cu_nan );

  output_tree.Branch("jet_rnnip_pb",&m_branches->v_jet_rnnip_pb );
  output_tree.Branch("jet_rnnip_pc",&m_branches->v_jet_rnnip_pc );
  output_tree.Branch("jet_rnnip_pu",&m_branches->v_jet_rnnip_pu );
  output_tree.Branch("jet_rnnip_ptau",&m_branches->v_jet_rnnip_ptau);


}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void ImpactParameterBranches::fill(const xAOD::Jet& jet) {

    double ip2d_pb   = m_replaceDefaults ? m_defaultDict["ip2_pb"] : NAN;;
    double ip2d_pc   = m_replaceDefaults ? m_defaultDict["ip2_pc"] : NAN;;
    double ip2d_pu   = m_replaceDefaults ? m_defaultDict["ip2_pu"] : NAN;;

    float ip2        = m_replaceDefaults ? m_defaultDict["ip2"   ] : NAN;;
    float ip2_c      = m_replaceDefaults ? m_defaultDict["ip2_c" ] : NAN;;
    float ip2_cu     = m_replaceDefaults ? m_defaultDict["ip2_cu"] : NAN;;
    float ip2_nan    = NAN;
    float ip2_c_nan  = NAN;
    float ip2_cu_nan = NAN;


    double ip3d_pb   = m_replaceDefaults ? m_defaultDict["ip3_pb"] : NAN;;
    double ip3d_pc   = m_replaceDefaults ? m_defaultDict["ip3_pc"] : NAN;;
    double ip3d_pu   = m_replaceDefaults ? m_defaultDict["ip3_pu"] : NAN;;

    float ip3        = m_replaceDefaults ? m_defaultDict["ip3"   ] : NAN;;
    float ip3_c      = m_replaceDefaults ? m_defaultDict["ip3_c" ] : NAN;;
    float ip3_cu     = m_replaceDefaults ? m_defaultDict["ip3_cu"] : NAN;;
    float ip3_nan    = NAN;
    float ip3_c_nan  = NAN;
    float ip3_cu_nan = NAN;

    const xAOD::BTagging *bjet = jet.btagging();

    std::vector< ElementLink< xAOD::TrackParticleContainer > > IP2DTracks;
    IP2DTracks = bjet->auxdata<std::vector<ElementLink< xAOD::TrackParticleContainer> > >("IP2D_TrackParticleLinks");
    if (IP2DTracks.size() > 0) {

      ip2d_pb = bjet->IP2D_pb();
      ip2d_pc = bjet->IP2D_pc();
      ip2d_pu = bjet->IP2D_pu();

      ip2    = bjet->calcLLR(ip2d_pb,ip2d_pu);
      ip2_c  = bjet->calcLLR(ip2d_pb,ip2d_pc);
      ip2_cu = bjet->calcLLR(ip2d_pc,ip2d_pu);

      if(ip2d_pb<=0. or ip2d_pu<=0. or ip2d_pb==NAN or ip2d_pu==NAN) {
          ip2_nan = NAN;
      } else {
          ip2_nan = log(ip2d_pb/ip2d_pu);
      }

      if(ip2d_pb<=0. or ip2d_pc<=0. or ip2d_pb==NAN or ip2d_pc==NAN) {
          ip2_c_nan = NAN;
      } else {
          ip2_c_nan = log(ip2d_pb/ip2d_pc);
      }

      if(ip2d_pc<=0. or ip2d_pu<=0. or ip2d_pc==NAN or ip2d_pu==NAN) {
          ip2_cu_nan = NAN;
      } else {
          ip2_cu_nan = log(ip2d_pc/ip2d_pu);
      }

      m_branches->v_jet_ip2d_pb->push_back(ip2d_pb);
      m_branches->v_jet_ip2d_pc->push_back(ip2d_pc);
      m_branches->v_jet_ip2d_pu->push_back(ip2d_pu);
      m_branches->v_jet_ip2d_llr->push_back(bjet->IP2D_loglikelihoodratio());

      m_branches->v_jet_ip2->push_back( ip2        );
      m_branches->v_jet_ip2_c->push_back( ip2_c      );
      m_branches->v_jet_ip2_cu->push_back( ip2_cu     );
      m_branches->v_jet_ip2_nan->push_back( ip2_nan    );
      m_branches->v_jet_ip2_c_nan->push_back( ip2_c_nan  );
      m_branches->v_jet_ip2_cu_nan->push_back( ip2_cu_nan );


    }
    else {
      m_branches->v_jet_ip2d_pb->push_back(-99);
      m_branches->v_jet_ip2d_pc->push_back(-99);
      m_branches->v_jet_ip2d_pu->push_back(-99);
      m_branches->v_jet_ip2d_llr->push_back(-99);

      m_branches->v_jet_ip2->push_back( ip2        );
      m_branches->v_jet_ip2_c->push_back( ip2_c      );
      m_branches->v_jet_ip2_cu->push_back( ip2_cu     );
      m_branches->v_jet_ip2_nan->push_back( ip2_nan    );
      m_branches->v_jet_ip2_c_nan->push_back( ip2_c_nan  );
      m_branches->v_jet_ip2_cu_nan->push_back( ip2_cu_nan );
    }

    // IP3D
    std::vector< ElementLink< xAOD::TrackParticleContainer > > IP3DTracks;
    IP3DTracks = bjet->auxdata<std::vector<ElementLink< xAOD::TrackParticleContainer> > >("IP3D_TrackParticleLinks");
    if (IP3DTracks.size()) {
      ip3d_pb = bjet->IP3D_pb();
      ip3d_pc = bjet->IP3D_pc();
      ip3d_pu = bjet->IP3D_pu();

      ip3    = bjet->calcLLR(ip3d_pb,ip3d_pu);
      ip3_c  = bjet->calcLLR(ip3d_pb,ip3d_pc);
      ip3_cu = bjet->calcLLR(ip3d_pc,ip3d_pu);

      if(ip3d_pb<=0. or ip3d_pu<=0. or ip3d_pb==NAN or ip3d_pu==NAN) {
          ip3_nan = NAN;
      } else {
          ip3_nan = log(ip3d_pb/ip3d_pu);
      }

      if(ip3d_pb<=0. or ip3d_pc<=0. or ip3d_pb==NAN or ip3d_pc==NAN) {
          ip3_c_nan = NAN;
      } else {
          ip3_c_nan = log(ip3d_pb/ip3d_pc);
      }

      if(ip3d_pc<=0. or ip3d_pu<=0. or ip3d_pc==NAN or ip3d_pu==NAN) {
          ip3_cu_nan = NAN;
      } else {
          ip3_cu_nan = log(ip3d_pc/ip3d_pu);
      }

      m_branches->v_jet_ip3d_pb->push_back(ip3d_pb);
      m_branches->v_jet_ip3d_pc->push_back(ip3d_pc);
      m_branches->v_jet_ip3d_pu->push_back(ip3d_pu);
      m_branches->v_jet_ip3d_llr->push_back(bjet->IP3D_loglikelihoodratio());

      m_branches->v_jet_ip3->push_back( ip3        );
      m_branches->v_jet_ip3_c->push_back( ip3_c      );
      m_branches->v_jet_ip3_cu->push_back( ip3_cu     );
      m_branches->v_jet_ip3_nan->push_back( ip3_nan    );
      m_branches->v_jet_ip3_c_nan->push_back( ip3_c_nan  );
      m_branches->v_jet_ip3_cu_nan->push_back( ip3_cu_nan );


    }
    else {
      m_branches->v_jet_ip3d_pb->push_back(-99);
      m_branches->v_jet_ip3d_pc->push_back(-99);
      m_branches->v_jet_ip3d_pu->push_back(-99);
      m_branches->v_jet_ip3d_llr->push_back(-99);

      m_branches->v_jet_ip3->push_back( ip3        );
      m_branches->v_jet_ip3_c->push_back( ip3_c      );
      m_branches->v_jet_ip3_cu->push_back( ip3_cu     );
      m_branches->v_jet_ip3_nan->push_back( ip3_nan    );
      m_branches->v_jet_ip3_c_nan->push_back( ip3_c_nan  );
      m_branches->v_jet_ip3_cu_nan->push_back( ip3_cu_nan );
    }

    double rnnip_pb   = bjet->auxdata<double>("rnnip_pb")  ;
    double rnnip_pc   = bjet->auxdata<double>("rnnip_pc")  ;
    double rnnip_pu   = bjet->auxdata<double>("rnnip_pu")  ;

    //We are currently saving ptau for PFlow but not VR, so...
    double rnnip_ptau;
    try {
        bjet->auxdata<double>("rnnip_ptau");
    } catch (...) {
        rnnip_ptau = -99;
    }

    if(m_replaceDefaults){

      if(std::isnan(rnnip_pb  )){ rnnip_pb   = m_defaultDict["rnnip_pb"  ]; }
      if(std::isnan(rnnip_pc  )){ rnnip_pc   = m_defaultDict["rnnip_pc"  ]; }
      if(std::isnan(rnnip_pu  )){ rnnip_pu   = m_defaultDict["rnnip_pu"  ]; }
      if(std::isnan(rnnip_ptau)){ rnnip_ptau = m_defaultDict["rnnip_ptau"]; }

    }

    //RNNIP
    m_branches->v_jet_rnnip_pb->push_back(   rnnip_pb   );
    m_branches->v_jet_rnnip_pc->push_back(   rnnip_pc   );
    m_branches->v_jet_rnnip_pu->push_back(   rnnip_pu   );
    m_branches->v_jet_rnnip_ptau->push_back( rnnip_ptau     );
}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void ImpactParameterBranches::clear() {
  // clear vectors

  m_branches->v_jet_ip2d_pb->clear();
  m_branches->v_jet_ip2d_pc->clear();
  m_branches->v_jet_ip2d_pu->clear();
  m_branches->v_jet_ip2d_llr->clear();
  m_branches->v_jet_ip3d_pb->clear();
  m_branches->v_jet_ip3d_pc->clear();
  m_branches->v_jet_ip3d_pu->clear();
  m_branches->v_jet_ip3d_llr->clear();

  m_branches->v_jet_ip2->clear();
  m_branches->v_jet_ip2_c->clear();
  m_branches->v_jet_ip2_cu->clear();
  m_branches->v_jet_ip3->clear();
  m_branches->v_jet_ip3_c->clear();
  m_branches->v_jet_ip3_cu->clear();

  m_branches->v_jet_ip2_nan->clear();
  m_branches->v_jet_ip2_c_nan->clear();
  m_branches->v_jet_ip2_cu_nan->clear();
  m_branches->v_jet_ip3_nan->clear();
  m_branches->v_jet_ip3_c_nan->clear();
  m_branches->v_jet_ip3_cu_nan->clear();

  m_branches->v_jet_rnnip_pb->clear();
  m_branches->v_jet_rnnip_pc->clear();
  m_branches->v_jet_rnnip_pu->clear();
  m_branches->v_jet_rnnip_ptau->clear();
}






