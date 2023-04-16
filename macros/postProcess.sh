
InitStr="effPlot_FixedCutBEff_noJZ0_"

hadd -f ${InitStr}ALL_dijetPy8_tot.root ${InitStr}ALL_dijetPy8_JobN*of6.root
hadd -f ${InitStr}85_dijetPy8_tot.root ${InitStr}85_dijetPy8_JobN*of6.root
hadd -f ${InitStr}77_dijetPy8_tot.root ${InitStr}77_dijetPy8_JobN*of6.root
hadd -f ${InitStr}70_dijetPy8_tot.root ${InitStr}70_dijetPy8_JobN*of6.root
hadd -f ${InitStr}60_dijetPy8_tot.root ${InitStr}60_dijetPy8_JobN*of6.root
hadd -f ${InitStr}50_dijetPy8_tot.root ${InitStr}50_dijetPy8_JobN*of6.root

hadd -f ${InitStr}ALL_dijetHw_tot.root ${InitStr}ALL_dijetHw_JobN*of3.root
hadd -f ${InitStr}85_dijetHw_tot.root ${InitStr}85_dijetHw_JobN*of3.root
hadd -f ${InitStr}77_dijetHw_tot.root ${InitStr}77_dijetHw_JobN*of3.root
hadd -f ${InitStr}70_dijetHw_tot.root ${InitStr}70_dijetHw_JobN*of3.root
hadd -f ${InitStr}60_dijetHw_tot.root ${InitStr}60_dijetHw_JobN*of3.root
hadd -f ${InitStr}50_dijetHw_tot.root ${InitStr}50_dijetHw_JobN*of3.root

root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}ALL_dijetPy8_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,1,-1,-1,\"${InitStr}ALL_dijetPy8_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}85_dijetPy8_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,1,-1,-1,\"${InitStr}85_dijetPy8_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}77_dijetPy8_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,1,-1,-1,\"${InitStr}77_dijetPy8_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}70_dijetPy8_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,1,-1,-1,\"${InitStr}70_dijetPy8_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}60_dijetPy8_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,1,-1,-1,\"${InitStr}60_dijetPy8_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}50_dijetPy8_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,1,-1,-1,\"${InitStr}50_dijetPy8_tot.root\")" 

root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}ALL_dijetHw_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,0,-1,-1,\"${InitStr}ALL_dijetHw_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}85_dijetHw_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,0,-1,-1,\"${InitStr}85_dijetHw_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}77_dijetHw_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,0,-1,-1,\"${InitStr}77_dijetHw_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}70_dijetHw_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,0,-1,-1,\"${InitStr}70_dijetHw_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}60_dijetHw_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,0,-1,-1,\"${InitStr}60_dijetHw_tot.root\")" 
root -l -b -q "../macros/effPlot_newBinning.cxx+gO(\"${InitStr}50_dijetHw_eval\",-1,\"/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-09/dijet_py8/\",40000000,-999.,-999.,0,-1,-1,\"${InitStr}50_dijetHw_tot.root\")" 

