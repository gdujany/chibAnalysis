Per fare l'analisi occorre:

./efficiency -f --makeHistos -c1
./efficiency -f --makeHistos -c2
./efficiency -f --makeHistos --reweighPt 1S -c1
./efficiency -f --makeHistos --reweighPt 1S -c2
./efficiency -f --makeHistos --reweighPt 2S -c1
./efficiency -f --makeHistos --reweighPt 2S -c2

Crea le rootuple con gli istogrammi necessari per il calcolo delle efficienze nella cartella effHistos
Bisogna che in ../store/ ci siano anche le rootuple con i generati (di solito lo faccio su pccmsto04)
Ognuno ci mette circa due ore ma possono essere lanciati tutti insieme in un unico screen

Poi

cd sistematici
./variaCBparams.sh

che fa i 100 fit facendo variare i parametri delle DSCB entro gli errori. Produce i files nella cartella trees/ che poi vanno ultriormente uniti
es. hadd numChib_refit_20_40.root numChib_refit_20_40_?.root

Ci mette una ventina di ore su pccmsto04 e contiene il memory leak (ma ora dovrebbe girare senza danni)


Poi, eventualmente anche su una macchina locale:
cd ../AgreementDataMC
./agreementDataMC.py 
per 3 volte cambiando il valore della variabile globale 
mc_ptSpectrum = 1S
mc_ptSpectrum = 2S
mc_ptSpectrum = 3S

Poi lanciare, eventualmente anche su una macchina locale dopo aver copiato i files ottenuti nei passi precedenti

./vai.sh

che fa tutto il resto dell'analisi, ci mette circa mezzora e salva le tabelle in tables e i plots un po' in giro
per i dettagli sui vari passi aprire vai.sh

Per vedere dove sono le tabelle e i plots da includere nell'AN guardare copy.sh 