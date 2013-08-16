#!/usr/bin/python
from __future__ import division
from ROOT import gSystem
gSystem.Load('../My_double_CB/My_double_CB_cxx')
from ROOT import gROOT,  gStyle, TCanvas, TFile, TTree, TH1, TH1D, RooRandom, TLegend
from array import array
import sys
sys.path.append('..')
from pyUtils import *
from mcFit import doMCFit
from dataFit import doDataFit
import ROOT, getopt, time

ROOT.SetMemoryPolicy( ROOT.kMemoryStrict )

def makeTree(outFile_name, cuts, fittedVariable):

    rootupla_Chib1_Sel = '../../store/ChiB_1P_1_UpsilonPt.root'
    rootupla_Chib2_Sel = '../../store/ChiB_1P_2_UpsilonPt.root'
    rootupla_2012Data = "../../store/2012_AllData.root"

    mass_chib1_1P =  9.89278
    mass_chib2_1P =  9.91221
   

    rootFile = TFile(outFile_name, 'recreate')
    tree = TTree('numChib', 'tree with number of chibs')

    numChib1 = array('i', [0])
    numChib2 = array('i', [0])
    ratio = array('d', [0])
    
    mean1 = array('d', [0])
    sigma1 = array('d', [0])
    a1_1 = array('d', [0])
    a2_1 = array('d', [0])
        
    mean2 = array('d', [0])
    sigma2 = array('d', [0])
    a1_2 = array('d', [0])
    a2_2 = array('d', [0])
        
    reducedChi2 = array('d', [0])
    
    

    tree.Branch('numChib1', numChib1, 'numChib1/I')
    tree.Branch('numChib2', numChib2, 'numChib2/I')
    tree.Branch('ratio', ratio, 'ratio/D')

    tree.Branch('mean1', mean1, 'mean1/D')
    tree.Branch('sigma1', sigma1, 'sigma1/D')
    tree.Branch('a1_1', a1_1, 'a1_1/D')
    tree.Branch('a2_1', a2_1, 'a2_1/D')
    

    tree.Branch('mean2', mean2, 'mean2/D')
    tree.Branch('sigma2', sigma2, 'sigma2/D')
    tree.Branch('a1_2', a1_2, 'a1_2/D')
    tree.Branch('a2_2', a2_2, 'a2_2/D')
    
    tree.Branch('reducedChi2', reducedChi2, 'reducedChi2/D')

    

    Chib1_fitResult = doMCFit(inputfile_name=rootupla_Chib1_Sel, mass_chib=mass_chib1_1P, cuts=cuts, fittedVariable=fittedVariable, returnOnlyFitResult = True)
    Chib2_fitResult = doMCFit(inputfile_name=rootupla_Chib2_Sel, mass_chib=mass_chib2_1P, cuts=cuts, fittedVariable=fittedVariable, returnOnlyFitResult = True)

    # make dataset outside the loop only once to speed up the process
    print "Creating DataSet from file "+str(rootupla_2012Data)
    dataSet = makeRooDataset(rootupla_2012Data)
    RooRandom.randomGenerator().SetSeed(int(time.time()*256))

    for i in range(100):
        print '\nstarting data fit number ' + str(i) + '\n'
        Chib1_listParams = Chib1_fitResult.randomizePars()
        Chib2_listParams = Chib2_fitResult.randomizePars()

        Chib1_parameters = CB_parameters(a1 = Chib1_listParams[0].getVal(),
                                         a2 = Chib1_listParams[1].getVal(),
                                         mean = Chib1_listParams[2].getVal(),
                                         n1 = 2.5,
                                         n2 = 3,
                                         sigma = Chib1_listParams[3].getVal(),
                                         s_sigma = Chib1_listParams[3].getError())

        Chib2_parameters = CB_parameters(a1 = Chib2_listParams[0].getVal(),
                                         a2 = Chib2_listParams[1].getVal(),
                                         mean = Chib2_listParams[2].getVal(),
                                         n1 = 2.5,
                                         n2 = 3,
                                         sigma = Chib2_listParams[3].getVal(),
                                         s_sigma = Chib2_listParams[3].getError())

        numChib, chi2r = doDataFit(RooDataSet = dataSet,Chib1_parameters=Chib1_parameters,Chib2_parameters=Chib2_parameters, cuts=cuts, fittedVariable=fittedVariable, noPlots = True)

        numChib1[0] = int(numChib.numChib1)
        numChib2[0] = int(numChib.numChib2)
        ratio[0] = numChib.calcRatio()
       
        a1_1[0] = Chib1_listParams[0].getVal()
        a2_1[0] = Chib1_listParams[1].getVal()
        mean1[0] = Chib1_listParams[2].getVal()
        sigma1[0] = Chib1_listParams[3].getVal()
        
        a1_2[0] = Chib2_listParams[0].getVal()
        a2_2[0] = Chib2_listParams[1].getVal()
        mean2[0] = Chib2_listParams[2].getVal()
        sigma2[0] = Chib2_listParams[3].getVal()

        reducedChi2[0] = chi2r
        
        tree.Fill()

    rootFile.Write()
    rootFile.Close()
   


def calcSistematic(varLabel, ptBin):

    ptBin_label = '' if ptBin == None else '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])        
    inFile_name = 'trees/numChib_'+varLabel+ptBin_label+'.root'
    outFile_name = 'plotsVariaCB/ratios_'+varLabel+ptBin_label

       
    h_ratios = TH1D('ratios','Distribution of Num #chi_{b2}/Num #chi_{b1} changing DSCB parameters with covariance matrix;N_{2}/N_{1};frequency',30,1,1)# 0.48, 0.67)# 0.56, 0.68)
    # for refit allData use limits histo 0.52, 0.62, for bin2 0.5, 0.7
    # for qValue 

    inFile = TFile(inFile_name,'read')
    tree = TTree()
    inFile.GetObject('numChib', tree)
    
    def fillHisto(h_ratios, tree):
        for event in tree:
            if event.reducedChi2 < 1.5:
                h_ratios.Fill(event.ratio)
        fitResult = h_ratios.Fit('gaus','S')
        return fitResult

    fitResult = fillHisto(h_ratios, tree)
    mean, sigma = fitResult.Parameter(1), fitResult.Parameter(2)
    
    # h_ratios = TH1D('ratios','Distribution of Num #chi_{b2}/Num #chi_{b1} changing DSCB parameters with covariance matrix;Num #chi_{b2}/Num #chi_{b1};',30,mean-5.0*sigma,mean+5.0*sigma)

    # fitResult = fillHisto(h_ratios, tree)
    # mean, sigma = fitResult.Parameter(1), fitResult.Parameter(2)

    h_ratios = TH1D('ratios','Distribution of Num #chi_{b2}/Num #chi_{b1} changing DSCB parameters with covariance matrix;N_{2}/N_{1};frequency',30,mean-3.5*sigma,mean+3.5*sigma)
    fitResult = fillHisto(h_ratios, tree)
    inFile.Close

    print 'chi2/dof = ' + str(fitResult.Chi2()) + ' / ' + str(fitResult.Ndf())
    print 'chi2 probability (p-value) = ' + str(fitResult.Prob())
    print 'mean = ' + str(fitResult.Parameter(1))
    print 'sigma = ' + str(fitResult.Parameter(2))
    
    canvas = TCanvas("canvas","canvas")
    h_ratios.Draw()

    legend = TLegend(.05,.78,.55,.99)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    ptMin = str(ptBins[0]) if ptBin == None else str(ptBins[ptBin-1])
    ptMax = str(ptBins[-1]) if ptBin == None else str(ptBins[ptBin])
    legend.AddEntry(0,ptMin+' GeV < p_{T}(#Upsilon) < '+ptMax+' GeV','')
    legend.Draw('same')
    canvas.Update()
       
    canvas.SaveAs(outFile_name+'.png')
    canvas.SaveAs(outFile_name+'.root')
    return fitResult.Parameter(2)
    

if __name__ == '__main__':

    #gROOT.SetStyle('Plain') # white background
    #gStyle.SetOptStat(222211)
    # gStyle.SetOptStat(0)
    gStyle.SetOptFit(101)

    ROOT.gROOT.ProcessLine('.L ../tdrstyle.C')
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style
    gStyle.SetOptFit(11)
    gStyle.SetOptStat(10)
    gStyle.SetStatY(0.99)
    gStyle.SetStatX(0.99)
    gStyle.SetStatW(0.2)
    gStyle.SetStatH(0.4)


    ptBins = [7,11,16,20,40]
    
    # GRAB OPTIONS
    # Defaults
    fittedVariable='qValue'
    ptBin = None
    makeAll = False
    isMakeTree = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],'fab:',['refit', 'all', 'ptBin=', 'makeTree='])
    except getopt.GetoptError:
        print './variaCBparams.py [-f] [-a] [-b <1, 2, 3, 4>]'
        print '-f = --refit'  
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-f', '--refit'):
            fittedVariable='refittedMass'
        if opt in ('-a', '--all'):
            makeAll = True
        if opt in ('-b', '--ptBin'):
            ptBin = int(arg) 
        if opt == '--makeTree':
            isMakeTree = True
            if str(arg) == '':
                addTreeName = str(arg) 
            else:
                addTreeName = '_'+str(arg)
            
    
    varLabel = 'refit' 
    if fittedVariable != 'refittedMass': varLabel = 'qValue'

  
    def writeSistematic(varLabel, ptBin):
        
        error = calcSistematic(varLabel, ptBin)

        ptBin_label = '' 
        if ptBin != None:
            ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin]) 
        sistFile_name = '../txt/sistematici_'+varLabel+ptBin_label+'.txt'
        sistematici = Sistematici()
        sistematici = sistematici.loadFromFile(sistFile_name)
        sistematici.dscbParam = error/sistematici.value
        sistematici.saveToFile(sistFile_name)

###################

    if isMakeTree:

        cuts = Cuts()
        cuts = cuts.loadFromFile("../cuts.txt")

        if ptBin == None:
            ptBin_label = ''
            cuts.upsilon_pt_lcut = ptBins[0]
            cuts.upsilon_pt_hcut = ptBins[-1]
        else:
            ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
            cuts.upsilon_pt_lcut = ptBins[ptBin-1]
            cuts.upsilon_pt_hcut = ptBins[ptBin]
        
        if fittedVariable == 'refittedMass': cuts.varLabel = 'rf1S_'
        print str(cuts)
    
        outFile_name = 'trees/numChib_'+varLabel+ptBin_label+addTreeName+'.root'
        makeTree(outFile_name, cuts, fittedVariable)

    elif makeAll:
        for ptBin in [None] + range(1,len(ptBins)):
            varLabel = 'refit'
            writeSistematic(varLabel, ptBin)
    else:
        writeSistematic(varLabel, ptBin)
        
        
