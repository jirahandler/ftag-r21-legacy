#mv2c10 working points listed here:
#https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarksRelease21#MV2c10_tagger
#FixedCutBEff_60	0.94	61.14	22	150	1204
#FixedCutBEff_70	0.83	70.84	8	39	313
#FixedCutBEff_77	0.64	77.53	4	16	113
#FixedCutBEff_85	0.11	85.23	2	6	28


root -l -b -q "../../macros/effPlot_eos.cxx+gO(\"effPlot_FixedCutBEff_ttbar-PF-test\",0.11,\"root://hepgrid11.ph.liv.ac.uk//dpm/ph.liv.ac.uk/home/atlas/atlaslocalgroupdisk/rucio/user/awychan/adjMC/PFlow/FTAG1-ttbar/\",-1)"
root -b -q -l "../../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_ttbar-PF-test.root\",\"pt85_ttbar-PF-test\", \"pt_LF\")"
root -b -q -l "../../macros/extract_prj.cxx(\"pt85_ttbar-PF-test.root\",\"LF_WP85_ttbar-PF\")"