#ifndef IMPACTPARAMETER_BRANCH_BUFFER_HH
#define IMPACTPARAMETER_BRANCH_BUFFER_HH

#include <vector>

struct ImpactParameterBranchBuffer {

	std::vector<float> *v_jet_ip2d_pb ;
	std::vector<float> *v_jet_ip2d_pc ;
	std::vector<float> *v_jet_ip2d_pu ;
	std::vector<float> *v_jet_ip2d_llr;
	std::vector<float> *v_jet_ip3d_pb ;
	std::vector<float> *v_jet_ip3d_pc ;
	std::vector<float> *v_jet_ip3d_pu ;
	std::vector<float> *v_jet_ip3d_llr;

	std::vector<float> *v_jet_ip2;
	std::vector<float> *v_jet_ip2_c;
	std::vector<float> *v_jet_ip2_cu;
	std::vector<float> *v_jet_ip3;
	std::vector<float> *v_jet_ip3_c;
	std::vector<float> *v_jet_ip3_cu;

	std::vector<float> *v_jet_ip2_nan;
	std::vector<float> *v_jet_ip2_c_nan;
	std::vector<float> *v_jet_ip2_cu_nan;
	std::vector<float> *v_jet_ip3_nan;
	std::vector<float> *v_jet_ip3_c_nan;
	std::vector<float> *v_jet_ip3_cu_nan;

	std::vector<float> *v_jet_rnnip_pb;
	std::vector<float> *v_jet_rnnip_pc;
	std::vector<float> *v_jet_rnnip_pu;
	std::vector<float> *v_jet_rnnip_ptau;
};

#endif // IMPACTPARAMETER_BRANCH_BUFFER_HH
