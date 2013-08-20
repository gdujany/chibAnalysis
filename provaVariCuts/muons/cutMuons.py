#!/usr/bin/python
from __future__ import division
from ROOT import gSystem, gROOT, gStyle, TFile, TTree, TH1, TH1D,TH2D, TCanvas, THStack
import ROOT, sys, getopt, math
sys.path.append('../..')
sys.path.append('../../Polarization')
from pyUtils import *
from polarization import *

rootupla_gen_1 = '../../../store/GenParticles_ChiB_1P_1_Upsilon2SPt.root'
rootupla_gen_2 = '../../../store/GenParticles_ChiB_1P_2_Upsilon2SPt.root'


outputFile_name = 'cutMuons'

cuts = Cuts()
cuts = cuts.loadFromFile("../../cuts.txt")

def printAcceptance():
    def calcAcceptance(chib_state):
        rootupla = rootupla_gen_1 if chib_state == 1 else rootupla_gen_2
        in_file = TFile.Open(rootupla,'READ')
        tree = TTree()
        in_file.GetObject('rootuple/GenParticlesTree', tree)

        N0 = 0
        N1 = 0

        cont = 0
        for event in tree:
            if cont == 10000:
                break
            if event.Upsilon_pt > cuts.upsilon_pt_lcut and event.Upsilon_pt < cuts.upsilon_pt_hcut and abs(event.photon_eta) < cuts.photon_eta_cut and abs(event.Upsilon_rapidity) < cuts.upsilon_rapidity_cut and abs(event.muP_p4.Eta()) < 2.4 and abs(event.muM_p4.Eta()) < 2.4:
                muP_pt = event.muP_p4.Pt()
                muP_eta = event.muP_p4.Eta()
                ups_dir, mu_dir = upsilonMuDirections(event.chib_p4, event.Upsilon_p4, event.muP_p4,'hx')
                weight = angDist(ups_dir, mu_dir, chib_state, None)
                N0 += weight
                cont += 1
                if  abs(muP_eta) < 2.4 and abs(event.muM_p4.Eta()) < 2.4 and muP_pt > 2.5 and event.muM_p4.Pt() > 2.5:
                    N1 += weight

        return N0, N1


    N0_1, N1_1 = calcAcceptance(1)
    N0_2, N1_2 = calcAcceptance(2)

    print 'Acceptance due to muon selection for chi_b1 = ' + str(N1_1/N0_1)
    print 'Acceptance due to muon selection for chi_b2 = ' + str(N1_2/N0_2)



def makePlots():
    
    
    histos_chib1 = dict()
    histos_chib2 = dict()

    histos_chib1['mu_pt'] = TH1D('mu_pt_chib1', '#chi_{b1};p_{T}(#mu);',60, 0, 30)
    histos_chib1['mu_eta'] = TH1D('mu_eta_chib1', '#chi_{b1};#eta(#mu);',60, -3, 3)

    histos_chib2['mu_pt'] = TH1D('mu_pt_chib2', '#chi_{b2};p_{T}(#mu);',60, 0, 30)
    histos_chib2['mu_eta'] = TH1D('mu_eta_chib2', '#chi_{b2};#eta(#mu);',60, -3, 3)

    h2D_chib1 = TH2D('pt_vs_eta_chib1', '#chi_{b1};#eta(#mu);p_{T}(#mu)',60, -3, 3, 60, 0, 30)
    h2D_chib2 = TH2D('pt_vs_eta_chib2', '#chi_{b2};#eta(#mu);p_{T}(#mu)',60, -3, 3, 60, 0, 30)

    titles = dict( mu_pt = ';p_{T}(#mu);',
                   mu_eta = ';#eta(#mu);')


    def fillHistos(chib_state, histos, h2D):
        rootupla = rootupla_gen_1 if chib_state == 1 else rootupla_gen_2
        in_file = TFile.Open(rootupla,'READ')
        tree = TTree()
        in_file.GetObject('rootuple/GenParticlesTree', tree)

        cont = 0
        for event in tree:
            if cont == 10000:
                break
            if event.Upsilon_pt > cuts.upsilon_pt_lcut  and event.Upsilon_pt < cuts.upsilon_pt_hcut and abs(event.photon_eta) < cuts.photon_eta_cut and abs(event.Upsilon_rapidity) < cuts.upsilon_rapidity_cut:# and abs(event.muP_p4.Eta()) < 2.4 and abs(event.muM_p4.Eta()) < 2.4:
                muP_pt = event.muP_p4.Pt()
                muP_eta = event.muP_p4.Eta()
                ups_dir, mu_dir = upsilonMuDirections(event.chib_p4, event.Upsilon_p4, event.muP_p4,'hx')
                weight = angDist(ups_dir, mu_dir, chib_state, None)
                histos['mu_pt'].Fill(muP_pt, weight)
                histos['mu_pt'].Fill(event.muM_p4.Pt(), weight)
                histos['mu_eta'].Fill(muP_eta, weight)
                histos['mu_eta'].Fill(event.muM_p4.Eta(), weight)
                h2D.Fill(muP_eta, muP_pt, weight)
                h2D.Fill(event.muM_p4.Eta(), event.muM_p4.Pt(), weight)
                cont += 1

    fillHistos(1, histos_chib1, h2D_chib1)
    fillHistos(2, histos_chib2, h2D_chib2)

    #for histo in histos_chib1.values()+histos_chib2.values():
    #    histo.Scale(1./histo.Integral())   

    for histo in histos_chib1.values():
        histo.SetLineColor(ROOT.kBlue)
        histo.SetLineWidth(2)
        histo.SetMarkerColor(ROOT.kBlue)
        histo.SetMarkerStyle(20)

    for histo in histos_chib2.values():
        histo.SetLineColor(ROOT.kRed)
        histo.SetLineWidth(2)
        histo.SetMarkerColor(ROOT.kRed)
        histo.SetMarkerStyle(21)



    canvas=TCanvas("canvas","canvas")
    canvas.Print(outputFile_name+".ps[") #apre file .ps
    for key in histos_chib1.keys():
        hs = THStack('hs', titles[key])
        hs.Add(histos_chib1[key])
        hs.Add(histos_chib2[key])
        hs.Draw('nostack')
        leg = canvas.BuildLegend(.65,.93,.95,.7)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw("SAME")
        canvas.Update()
        canvas.Print(outputFile_name+'.ps')
        canvas.Print(key+'.png')
        h2D_chib1.Draw('colz')
        canvas.Print(outputFile_name+'.ps')
        h2D_chib2.Draw('colz')
        canvas.Print(outputFile_name+'.ps')
        canvas.Print(outputFile_name+".ps]") #chiude file .ps


if __name__ == '__main__':
    #printAcceptance()
    makePlots()

