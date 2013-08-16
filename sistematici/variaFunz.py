#!/usr/bin/python
from __future__ import division
from ROOT import gSystem
gSystem.Load('../My_double_CB/My_double_CB_cxx')
from ROOT import gROOT,  gStyle, TCanvas, TFile, TTree, TH1, TH1D
from array import array
import sys
sys.path.append('..')
from pyUtils import *
from mcFit import doMCFit
from dataFit import doDataFit
import ROOT, getopt

ptBins = [7,11,16,20,40]

rootupla_Chib1_Sel = '../../store/ChiB_1P_1_UpsilonPt.root'
rootupla_Chib2_Sel = '../../store/ChiB_1P_2_UpsilonPt.root'
rootupla_2012Data = "../../store/2012_AllData.root"

# make dataset outside the loop only once to speed up the process
print "Creating DataSet from file "+str(rootupla_2012Data)
dataSet = makeRooDataset(rootupla_2012Data)

def vai(varLabel = 'qValue', ptBin = None): 

    fittedVariable = 'refittedMass' if varLabel == 'refit' else 'qValue' 

    mass_chib1_1P =  9.89278
    mass_chib2_1P =  9.91221

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

    outputFile_name = 'plotsVariaFunz/variaFunz_'+varLabel+ptBin_label+'.ps'
    canvas=TCanvas("canvas","canvas")
    canvas.Print(outputFile_name+"[")

    Chib1_parameters, canvas =  doMCFit(inputfile_name=rootupla_Chib1_Sel, mass_chib=mass_chib1_1P, cuts=cuts, output_name='ChiB1', plotTitle='#chi_{b1}', fittedVariable=fittedVariable,useOtherSignalParametrization = False)
    canvas.Print(outputFile_name)
    #canvas.Print('plotsVariaFunz/chib1_'+varLabel+ptBin_label+'.pdf')

    Chib2_parameters, canvas =  doMCFit(inputfile_name=rootupla_Chib2_Sel, mass_chib=mass_chib2_1P, cuts=cuts, output_name='ChiB2', plotTitle='#chi_{b2}', fittedVariable=fittedVariable,useOtherSignalParametrization = False)
    canvas.Print(outputFile_name)
    #canvas.Print('plotsVariaFunz/chib2_'+varLabel+ptBin_label+'.pdf')

# Dataset made outside    

    numChib_DSCB_D0Background, canvas = doDataFit(RooDataSet=dataSet,Chib1_parameters=Chib1_parameters,Chib2_parameters=Chib2_parameters, cuts=cuts, fittedVariable=fittedVariable, plotTitle = 'Default', titleOnPlot=True, cmsOnPlot=False, useOtherSignalParametrization = False)
    canvas.Print(outputFile_name)
    canvas.Print('plotsVariaFunz/default_'+varLabel+ptBin_label+'.pdf')
    canvas.Print('plotsVariaFunz/default_'+varLabel+ptBin_label+'.root')

    numChib_CBGauss_D0Background, canvas = doDataFit(RooDataSet=dataSet,Chib1_parameters=Chib1_parameters,Chib2_parameters=Chib2_parameters, cuts=cuts, fittedVariable=fittedVariable, plotTitle = 'Alternative signal', titleOnPlot=True, cmsOnPlot=False, massFreeToChange = True, sigmaFreeToChange = True)
    canvas.Print(outputFile_name)
    canvas.Print('plotsVariaFunz/changeSgn_'+varLabel+ptBin_label+'.pdf')
    canvas.Print('plotsVariaFunz/changeSgn_'+varLabel+ptBin_label+'.root')

    numChib_DSCB_ChebBackground, canvas = doDataFit(RooDataSet=dataSet,Chib1_parameters=Chib1_parameters,Chib2_parameters=Chib2_parameters, cuts=cuts, fittedVariable=fittedVariable, plotTitle = 'Alternative background', titleOnPlot=True, cmsOnPlot=False, useOtherSignalParametrization = False, useOtherBackgroundParametrization = True)
    canvas.Print(outputFile_name)
    canvas.Print('plotsVariaFunz/changeBkg_'+varLabel+ptBin_label+'.pdf')
    canvas.Print('plotsVariaFunz/changeBkg_'+varLabel+ptBin_label+'.root')

    numChib_CBGauss_ChebBackground, canvas = doDataFit(RooDataSet=dataSet,Chib1_parameters=Chib1_parameters,Chib2_parameters=Chib2_parameters, cuts=cuts, fittedVariable=fittedVariable, plotTitle = 'Alternative signal and background', titleOnPlot=True, cmsOnPlot=False,  useOtherBackgroundParametrization = True, massFreeToChange = True, sigmaFreeToChange = True)
    canvas.Print(outputFile_name)
    canvas.Print('plotsVariaFunz/changeSgnBkg_'+varLabel+ptBin_label+'.pdf')
    canvas.Print('plotsVariaFunz/changeSgnBkg_'+varLabel+ptBin_label+'.root')


    canvas.Print(outputFile_name+"]")

    numChibs = [numChib_DSCB_D0Background, numChib_CBGauss_D0Background, numChib_DSCB_ChebBackground, numChib_CBGauss_ChebBackground]

      
    ratioDefault = numChib_DSCB_D0Background.ratio_21

    ratios = [numChib_CBGauss_D0Background.ratio_21, numChib_DSCB_ChebBackground.ratio_21, numChib_CBGauss_ChebBackground.ratio_21]
    
    maxDiff = 0
    for ratio in ratios:
        diff = abs(ratioDefault-ratio)
        if diff > maxDiff:
            maxDiff = diff

    return maxDiff, maxDiff/(ratioDefault), numChibs #error and relative error


if __name__ == '__main__':
    
    #gROOT.SetStyle('Plain') # white background
#gStyle.SetOptStat(222211)
    gStyle.SetOptFit(0101)

    ROOT.gROOT.ProcessLine('.L ../tdrstyle.C')
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style

    # GRAB OPTIONS
    # Defaults
    printAllData = False
    fittedVariable='qValue'
    ptBin = None
    makeAll = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],'fab:',['refit', 'all', 'ptBin='])
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
    
    varLabel = 'refit' if fittedVariable == 'refittedMass' else 'qValue'
    
    

    if not makeAll:
        err, rel_err, numChibs = vai(varLabel = varLabel, ptBin = ptBin)
        print 'error =          '+str(err)
        print 'relative error = '+str(rel_err)
    else:
        numChibs_dict = dict()
        with open('SistematiciVariaFunz.txt','w') as outFile:
            for varLabel in ('refit',):# 'qValue'):
                for ptBin in [None] + range(1,len(ptBins)):
                    err, rel_err, numChibs = vai(varLabel = varLabel, ptBin = ptBin)
                    outFile.write('\n'+varLabel + 'ptBin = '+str(ptBin)+'\n')
                    outFile.write('error =          '+str(err)+'\n'+'relative error = '+str(rel_err))

                    # write it to sistematici.txt
                    ptBin_label = '' if ptBin == None else '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
                    file_name = '../txt/sistematici_'+varLabel+ptBin_label+'.txt'
                    sistematici = Sistematici()
                    sistematici = sistematici.loadFromFile(file_name)
                    sistematici.funzParam = rel_err
                    sistematici.saveToFile(file_name)
                    
                    if varLabel == 'refit':
                        numChibs_dict[ptBin] = numChibs

        pickle.dump(numChibs_dict, open('numChibs_dict.pk','wb'))

        # write table
        varLabel = 'refit'
        rowsTitles = ['Default', 'Alternative signal', 'Alternative background', 'Alternative signal and background']
        with open('../tables/variaFunz_'+varLabel+'.tex','w') as outFile:
            numColumns = len(ptBins)
            if printAllData:
                numColumns += 1
            outFile.write('\\begin{center} \n \\begin{tabular}{p{0.25\linewidth}*{'+str(numColumns-1)+'}{c}}\n\\toprule\n')       
            outFile.write('$p_T(\\Upsilon)$ [GeV] ')
            for ptBin in range(1, len(ptBins)):
                outFile.write(' & '+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin])+' ')
            if printAllData:
                outFile.write(' & '+str(ptBins[0])+' - '+str(ptBins[-1])+' ')
            outFile.write('\\\\ \n \\midrule \n')
            outFile.write(' Fitting strategy &  \\multicolumn{'+str(numColumns-1)+'}{c}{$N_{2}/N_{1}$} \\\\ \n \\midrule \n')
            for riga in range(len(rowsTitles)):
                outFile.write(rowsTitles[riga])
                for ptBin in range(1, len(ptBins))+([None] if printAllData else []):
                    val, err = roundPair(numChibs_dict[ptBin][riga].ratio_21, numChibs_dict[ptBin][riga].s_ratio_21, 2)
                    outFile.write(' & '+val+' $\\pm$ '+err+' ')
                outFile.write('\\\\ \n')
            outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')
                
        
