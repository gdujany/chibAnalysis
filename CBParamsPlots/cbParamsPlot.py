#!/usr/bin/python
from __future__ import division
from ROOT import gSystem, gROOT,  gStyle, TH1, TH1D, TCanvas, THStack, TLine, TGaxis
import ROOT, sys
sys.path.append('..')
from pyUtils import *
from array import array

ROOT.gROOT.ProcessLine('.L ../tdrstyle.C')
ROOT.gROOT.Reset()
ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style
TGaxis.SetMaxDigits(3)

ptBins = [7, 11, 16, 20, 40]

ptLabel = 'UpsilonPt'
outputFile_name = 'cbParamsPlot.ps'

Chib1_params = []
Chib2_params = []
for ptBin in range(len(ptBins)-1):
    ptBin_label = '_'+str(ptBins[ptBin])+'_'+str(ptBins[ptBin+1])
    params1 = CB_parameters()
    params2 = CB_parameters()
    params1 = params1.loadFromFile('../txt/CB_ChiB1{0}_refit.txt'.format(ptBin_label))
    params2 = params2.loadFromFile('../txt/CB_ChiB2{0}_refit.txt'.format(ptBin_label))
    Chib1_params.append(params1)
    Chib2_params.append(params2)

Chib1_params_all = CB_parameters()
Chib1_params_all = Chib1_params_all.loadFromFile('../txt/CB_ChiB1_refit.txt')
Chib2_params_all = CB_parameters()
Chib2_params_all = Chib2_params_all.loadFromFile('../txt/CB_ChiB2_refit.txt')

keys = ['mean', 'sigma', 'a1', 'a2']#, 'n1', 'n2']
axisLabels = dict(
    mean = ';p_{T}(#Upsilon);#mu [GeV]',
    sigma = ';p_{T}(#Upsilon);#sigma [GeV]',
    a1 = ';p_{T}(#Upsilon);#alpha_{1}',
    a2 = ';p_{T}(#Upsilon);#alpha_{2}',
    #n1 = ';p_{T}(#Upsilon);n_{1}',
    #n2 = ';p_{T}(#Upsilon);n_{2}',
    )
estremiGrafico = dict(
    mean = [9.88, 9.92],
    sigma = [0., 0.01],
    a1 = [0., 1.],
    a2 = [1., 3.],
    #n1 = 0,
    #n2 = 0,
    )
h_chib1 = dict()
h_chib2 = dict()

for key in keys:
    h_chib1[key] = TH1D(key+'chib1', '#chi_{b1}'+axisLabels[key], len(ptBins)-1, array('d', ptBins))
    h_chib2[key] = TH1D(key+'chib2', '#chi_{b2}'+axisLabels[key], len(ptBins)-1, array('d', ptBins))
    for ptBin in range(len(ptBins)-1):
        h_chib1[key].SetBinContent(ptBin+1, getattr(Chib1_params[ptBin], key))
        h_chib1[key].SetBinError(ptBin+1, getattr(Chib1_params[ptBin], 's_'+key))
        h_chib2[key].SetBinContent(ptBin+1, getattr(Chib2_params[ptBin], key))
        h_chib2[key].SetBinError(ptBin+1, getattr(Chib2_params[ptBin], 's_'+key))
    

canvas=TCanvas("canvas","canvas")
canvas.Print(outputFile_name+"[") #apre file .ps

for key in keys:
    h_chib1[key].SetLineColor(ROOT.kBlue)
    h_chib1[key].SetMarkerColor(ROOT.kBlue)
    h_chib1[key].SetMarkerStyle(22)
    h_chib2[key].SetLineColor(ROOT.kRed)
    h_chib2[key].SetMarkerColor(ROOT.kRed)
    h_chib2[key].SetMarkerStyle(23)
    hs = THStack('hs', key+axisLabels[key])
    hs.Add(h_chib1[key])
    hs.Add(h_chib2[key])
    hs.SetMinimum(estremiGrafico[key][0])
    hs.SetMaximum(estremiGrafico[key][-1])
    hs.Draw('nostack')
    leg=canvas.BuildLegend(.8,.95,.92,.8)
    line1 = TLine(ptBins[0],getattr(Chib1_params_all, key),ptBins[-1],getattr(Chib1_params_all, key))
    line2 = TLine(ptBins[0],getattr(Chib2_params_all, key),ptBins[-1],getattr(Chib2_params_all, key))
    line1.SetLineColor(ROOT.kBlue)
    line2.SetLineColor(ROOT.kRed)
    line1.SetLineStyle(2)
    line2.SetLineStyle(2)
    print line1.GetTitle()
    line1.Draw('same')
    line2.Draw('same')
    leg.SetFillColor(0)
    leg.Draw("SAME")
    canvas.Update()
    canvas.Print(outputFile_name)
    canvas.Print('plots/'+key+'.pdf')
    canvas.Print('plots/'+key+'.root')

canvas.Print(outputFile_name+"]") #chiude file .ps
