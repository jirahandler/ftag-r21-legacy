#ifndef TAGGERSCORE_BRANCH_BUFFER_HH
#define TAGGERSCORE_BRANCH_BUFFER_HH

#include <vector>

struct TaggerScoreBranchBuffer {

	std::vector<double> *v_jet_dl1_pb;
	std::vector<double> *v_jet_dl1_pc;
	std::vector<double> *v_jet_dl1_pu;
	std::vector<double> *v_jet_dl1rmu_pb;
	std::vector<double> *v_jet_dl1rmu_pc;
	std::vector<double> *v_jet_dl1rmu_pu;
	std::vector<double> *v_jet_dl1r_pb;
	std::vector<double> *v_jet_dl1r_pc;
	std::vector<double> *v_jet_dl1r_pu;
	std::vector<double> *v_jet_mv2c10;
	std::vector<double> *v_jet_mv2c10mu;
	std::vector<double> *v_jet_mv2c10rnn;
	std::vector<double> *v_jet_mv2c100;
	std::vector<double> *v_jet_mv2cl100;

};

#endif // TAGGERSCORE_BRANCH_BUFFER_HH
