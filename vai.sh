#!/bin/bash
./mcFit.py -a

cd CBParamsPlots/
./cbParamsPlot.py
cd -

./dataFit.py -a

# NB imprtant to have done dataMCAgreement with 1S, 2S and 3S

./efficiency.py --makePlots -f --reweighPt 1S
./efficiency.py --makePlots -f --reweighPt 2S
./efficiency.py --makePlots -f --reweighPt 3S
./efficiency.py -a
./efficiency.py --makeTable

cd sistematici
./sistematici.py --setup
./variaCBparams.py -a
./variaFunz.py -a
cd -
./efficiency.py --sistematici -f
cd -
./sistematici.py --makeTable
cd ../provaVariCuts/
./PUbins.py
cd ..
./rcs.py --makeTable --makePlot

cd tables/
pdflatex vediTabelle
cd -

./dataFit.py -a --pas

