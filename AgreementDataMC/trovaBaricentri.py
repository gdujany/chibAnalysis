#!/usr/bin/python
from ROOT import gSystem, gROOT,  gStyle, TFile, TCanvas, TH1, TH1D
import ROOT, sys
sys.path.append('..')
sys.path.append('../Polarization')
from pyUtils import *
from polarization import *

TH1.SetDefaultSumw2(True) #Cosi` in automatico per tutti gli isto dovrebbe abilitare il calcolo degli errori come si deve (come somma dei pesi)

ptBins = [7,11,16,20,40]

inFile_name = 'agreementDataMC_2S_plots.root'

inFile = TFile(inFile_name, 'read')
h_dimuon_pt = inFile.Get('dimuon_pt_signal')

for ptBin in range(1,len(ptBins)):
    
    binMin = h_dimuon_pt.FindBin(ptBins[ptBin-1])
    binMax = h_dimuon_pt.FindBin(ptBins[ptBin])-1
    
    num = 0
    den = 0
    for i in range(binMin, binMax+1):
        num += h_dimuon_pt.GetBinCenter(i) * h_dimuon_pt.GetBinContent(i)
        den += h_dimuon_pt.GetBinContent(i)
    print num/den
