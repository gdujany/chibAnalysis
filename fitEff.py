#!/usr/bin/python
from ROOT import gSystem, gROOT, gStyle, TFile, TTree, TH1, TH1D, TCanvas, THStack, TF1
import ROOT, sys, getopt, math
sys.path.append('Polarization')
from efficiency import *
from array import array

ptBins = [7,11,16,20,40]

def makeEffFit(inFile1_name, inFile2_name):
    
    ratios = []
    s_ratios = []

    for ptBin in range(1, len(ptBins)):
        eff1 = makeEfficiency(inFile_name=inFile1_name, ptMin=ptBins[ptBin-1], ptMax=ptBins[ptBin])
        eff2 = makeEfficiency(inFile_name=inFile2_name, ptMin=ptBins[ptBin-1], ptMax=ptBins[ptBin])
        ratios[ptBin-1] = eff1.eff/eff2.eff
        s_ratios[ptBin-1] = ratio * sqrt((eff1.s_eff/eff1.eff)**2 + (eff2.s_eff/eff2.eff)**2)
    
    x = array('d', [8.71, 12.94, 17.51, 26.02])
    exl = array('d',[(x[ptBin-1]-ptBins[ptBin-1]) for ptBin in range(1, len(ptBins))])
    exh = array('d',[(ptBins[ptBin]-x[ptBin-1]) for ptBin in range(1, len(ptBins))])
    
    y = array('d',ratios)
    ey = array('d',s_ratios)
    graph = TGraphAsymmErrors(len(x), x, y, exl, exh, ey, ey) 
    
    canvas = TCanvas('canvas','canvas')
    graph.Draw()
    canvas.SaveAs('eff.png')
    


if __name__ == '__main__':
    
    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style

    varLabel = 'qValue'
    ptSpectrum_rw = '2S'
    try:
        opts, args = getopt.getopt(sys.argv[1:],'f',['refit', 'pt='])
    except getopt.GetoptError:
        print './fitEff.py [-f] [--pt]'
        print '-f = --refit'  
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-f', '--refit'):
            varLabel = 'refit'
        if opt == '--pt':
            ptSpectrum_rw = arg

    if ptSpectrum_rw in ('1S', '1s'):
        ptSpectrum_rw = '1S'
    elif ptSpectrum_rw in ('2S', '2s'):
        ptSpectrum_rw = '2S'
    elif ptSpectrum_rw in ('3S', '3S'):
        ptSpectrum_rw = '3S'

    inFile1_name = 'effHistos/chib1_hx_2S_'+ptSpectrum_rw+'_'+varLabel+'.root' 
    inFile2_name = 'effHistos/chib2_hx_2S_'+ptSpectrum_rw+'_'+varLabel+'.root'

    makeEffFit(inFile1_name, inFile2_name)
