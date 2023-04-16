#!/bin/bash

root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_50_dijetPy8_eval.root\",\"pt_LF_comb50\", \"pt_LF\",\"\",\"SF_vs_pt_wp50\",\"effPlot_FixedCutBEff_noJZ0t_50_dijetHw_eval.root\")"
root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_60_dijetPy8_eval.root\",\"pt_LF_comb60\", \"pt_LF\",\"\",\"SF_vs_pt_wp60\",\"effPlot_FixedCutBEff_noJZ0t_60_dijetHw_eval.root\")"
root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_70_dijetPy8_eval.root\",\"pt_LF_comb70\", \"pt_LF\",\"\",\"SF_vs_pt_wp70\",\"effPlot_FixedCutBEff_noJZ0t_70_dijetHw_eval.root\")"
root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_77_dijetPy8_eval.root\",\"pt_LF_comb77\", \"pt_LF\",\"\",\"SF_vs_pt_wp77\",\"effPlot_FixedCutBEff_noJZ0t_77_dijetHw_eval.root\")"
root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_85_dijetPy8_eval.root\",\"pt_LF_comb85\", \"pt_LF\",\"\",\"SF_vs_pt_wp85\",\"effPlot_FixedCutBEff_noJZ0t_85_dijetHw_eval.root\")"

root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_50_dijetPy8_eval.root\",\"pteta50_LF_comb\", \"etapt_LF\",\"SF_wp50.txt\",\"\",\"effPlot_FixedCutBEff_noJZ0t_50_dijetHw_eval.root\")"
root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_60_dijetPy8_eval.root\",\"pteta60_LF_comb\", \"etapt_LF\",\"SF_wp60.txt\",\"\",\"effPlot_FixedCutBEff_noJZ0t_60_dijetHw_eval.root\")"
root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_70_dijetPy8_eval.root\",\"pteta70_LF_comb\", \"etapt_LF\",\"SF_wp70.txt\",\"\",\"effPlot_FixedCutBEff_noJZ0t_70_dijetHw_eval.root\")"
root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_77_dijetPy8_eval.root\",\"pteta77_LF_comb\", \"etapt_LF\",\"SF_wp77.txt\",\"\",\"effPlot_FixedCutBEff_noJZ0t_77_dijetHw_eval.root\")"
root -b -q -l "../macros/combine_effects.cxx+g(\"effPlot_FixedCutBEff_noJZ0t_85_dijetPy8_eval.root\",\"pteta85_LF_comb\", \"etapt_LF\",\"SF_wp85.txt\",\"\",\"effPlot_FixedCutBEff_noJZ0t_85_dijetHw_eval.root\")"

root -b -q -l "../macros/extract_prj.cxx(\"pteta50_LF_comb.root\",\"LF_WP_Eff50\",1)"
root -b -q -l "../macros/extract_prj.cxx(\"pteta60_LF_comb.root\",\"LF_WP_Eff60\",1)"
root -b -q -l "../macros/extract_prj.cxx(\"pteta70_LF_comb.root\",\"LF_WP_Eff70\",1)"
root -b -q -l "../macros/extract_prj.cxx(\"pteta77_LF_comb.root\",\"LF_WP_Eff77\",1)"
root -b -q -l "../macros/extract_prj.cxx(\"pteta85_LF_comb.root\",\"LF_WP_Eff85\",1)"


# python ../macros/makeSystTable.py SF_wp85.txt 1 > AdjMC_SF_wp85_eta1.tex
# python ../macros/makeSystTable.py SF_wp85.txt 2 > AdjMC_SF_wp85_eta2.tex 
# python ../macros/makeSystTable.py SF_wp77.txt 1 > AdjMC_SF_wp77_eta1.tex
# python ../macros/makeSystTable.py SF_wp77.txt 2 > AdjMC_SF_wp77_eta2.tex 
# python ../macros/makeSystTable.py SF_wp70.txt 1 > AdjMC_SF_wp70_eta1.tex
# python ../macros/makeSystTable.py SF_wp70.txt 2 > AdjMC_SF_wp70_eta2.tex 
# python ../macros/makeSystTable.py SF_wp60.txt 1 > AdjMC_SF_wp60_eta1.tex
# python ../macros/makeSystTable.py SF_wp60.txt 2 > AdjMC_SF_wp60_eta2.tex 
# python ../macros/makeSystTable.py SF_wp50.txt 1 > AdjMC_SF_wp50_eta1.tex
# python ../macros/makeSystTable.py SF_wp50.txt 2 > AdjMC_SF_wp50_eta2.tex 
