#!/usr/bin/python
import sys, ROOT
sys.path.append('..')
from ROOT import gSystem
gSystem.Load('../My_double_CB/My_double_CB_cxx')
from ROOT import My_double_CB

from ROOT import  gROOT,  gStyle, TCanvas, TH1D, TF1, TGraphErrors
from array import array

from pyUtils import *
from mcFit import doMCFit
from dataFit import doDataFit

rootupla_Chib1_Sel = '../../store/ChiB_1P_1_UpsilonPt.root'
rootupla_Chib2_Sel = '../../store/ChiB_1P_2_UpsilonPt.root'
rootupla_2012Data = "../../store/2012_AllData.root"


if(len(sys.argv)>1 and sys.argv[1]=='2'):
    fittedVariable='qValue'
    label = 'qValue'
else:
    fittedVariable='refittedMass'
    label = 'refit'
    

mass_chib1_1P =  9.89278
mass_chib2_1P =  9.91221

# gROOT.SetStyle('Plain') # white background
# gStyle.SetOptStat(222211)
ROOT.gROOT.ProcessLine('.L ../tdrstyle.C')
ROOT.gROOT.Reset()
ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style




# for ptBin in (None, 1, 2):

#     if ptBin == None:
#         ptBin_label = ''
#         ptBinCut = ''
#     elif ptBin == 1:
#         ptBin_label = '_bin1'
#         ptBinCut = ' && dimuon_pt < 19'
#     elif ptBin == 2:
#         ptBin_label = '_bin2'
#         ptBinCut = ' && dimuon_pt > 19'


def makePlots():
    outputFile_name = 'PUbins_{0}.ps'.format(label)
    canvas=TCanvas("canvas","canvas")
    canvas.Print(outputFile_name+"[")

    lowEdgesOfBins = [0,11,14,18,40]

    lowEdgesOfBins.pop() # toglie 40. dalla lista cosi` ho solo gli inizi


    PU_cuts_labels = []
    PU_cuts_labels.append('# PV #leq {}'.format(lowEdgesOfBins[1]))
    for i in range(1,len(lowEdgesOfBins)-1):
        PU_cuts_labels.append('{} < # PV #leq {}'.format(lowEdgesOfBins[i],lowEdgesOfBins[i+1]))
    PU_cuts_labels.append('# PV > {}'.format(lowEdgesOfBins[len(lowEdgesOfBins)-1]))

    PU_cuts = []
    PU_cuts.append('numPrimaryVertices <= {}'.format(lowEdgesOfBins[1]))
    for i in range(1,len(lowEdgesOfBins)-1):
        PU_cuts.append('numPrimaryVertices > {} && numPrimaryVertices <= {}'.format(lowEdgesOfBins[i],lowEdgesOfBins[i+1]))
    PU_cuts.append('numPrimaryVertices > {}'.format(lowEdgesOfBins[len(lowEdgesOfBins)-1]))

    cuts = Cuts()
    cuts = cuts.loadFromFile("../cuts.txt")

    Chib1_parameters, canvas = doMCFit(inputfile_name=rootupla_Chib1_Sel, mass_chib=mass_chib1_1P, output_name='ChiB1', cuts=cuts, plotTitle='#chi_{b1}', fittedVariable=fittedVariable)
    canvas.Print(outputFile_name)
    Chib2_parameters, canvas = doMCFit(inputfile_name=rootupla_Chib2_Sel, mass_chib=mass_chib2_1P, output_name='ChiB2', cuts=cuts, plotTitle='#chi_{b2}', fittedVariable=fittedVariable)
    canvas.Print(outputFile_name)

    print "Creating DataSet from file "+str(rootupla_2012Data)
    dataSet = makeRooDataset(rootupla_2012Data)
    
    y = []
    ey = []
    
    for i in range(len(lowEdgesOfBins)):
        cuts.str_cut = '&& ' + PU_cuts[i] 
        plotTitle_add = PU_cuts_labels[i]         
        numChib, canvas = doDataFit(RooDataSet=dataSet,Chib1_parameters=Chib1_parameters,Chib2_parameters=Chib2_parameters, cuts=cuts, plotTitle ='#chi_{b}(1P)'+plotTitle_add, fittedVariable=fittedVariable)
        canvas.Print(outputFile_name)
        #canvas.Print('imgPU/{}_PU{}.png'.format(label,i+1))
        y.append(numChib.calcRatio())
        ey.append(numChib.calcRatioError())
        

    lowEdgesOfBins = lowEdgesOfBins + [40]
    x = array('d', [(lowEdgesOfBins[i-1]+lowEdgesOfBins[i])/2. for i in range(1, len(lowEdgesOfBins))])
    ex = array('d', [(x[i-1] - lowEdgesOfBins[i-1]) for i in range(1, len(lowEdgesOfBins))])

    y = array('d',y)
    ey = array('d',ey)
    
    
    h_ratios = TGraphErrors(len(x),x,y,ex,ey)
    h_ratios.SetTitle('Ratio #chi_{b2}/#chi_{b1};Number of primary vertices;N_{2}/N_{1}')
    

    gStyle.SetOptStat(0)
    h_ratios.SetMinimum(0)
    h_ratios.SetMaximum(1)
    h_ratios.Draw()
    canvas.Print(outputFile_name)
    # #canvas.Print('imgPU/{0}_PU.pdf'.format(label))
    # canvas.Print('imgPU/{0}_PU.root'.format(label))
    h_ratios.SaveAs('imgPU/{0}_PU.root'.format(label))

    canvas.Print(outputFile_name+"]")


def adjustPlot():

    ROOT.gROOT.ProcessLine('.L ../tdrstyle.C')
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine('setTDRStyle()') #Set CMS TDR style

    canvas=TCanvas("canvas","canvas")

    inputFile = TFile('imgPU/{0}_PU.root'.format(label),'READ')
    h_ratios = inputFile.Get('Graph')

    fa = TF1("fa","[0]",6,40)
    fa.SetParameter(0, 1)
    fa.SetParName(0,"Constant")

    h_ratios.Fit('fa')

    gStyle.SetOptStat(0)
    h_ratios.Draw('ap')  
   
    canvas.Print('imgPU/{0}_PU.pdf'.format(label))
    # canvas.Print('imgPU/{0}_PU.root'.format(label))
    
    

if __name__ == '__main__':
    makePlots()
    adjustPlot()
    
