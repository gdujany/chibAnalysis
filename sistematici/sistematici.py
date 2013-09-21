#!/usr/bin/python
from __future__ import division
import sys, getopt
sys.path.append('..')
from pyUtils import *


ptBins = [7,11,16,20,40]

def setupFileSistematici(varLabel, ptMin, ptMax, ptFileLabel, ptSpectrum_file):

    Chib1_eff_file = '../txt/eff_chib1'+'_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_file+ptBin_label+'.txt'
    Chib2_eff_file = '../txt/eff_chib2'+'_'+varLabel+'_'+ptFileLabel+'_'+ptSpectrum_file+ptBin_label+'.txt'
    numChib_file = '../txt/numChib_'+varLabel+ptBin_label+'.txt'
    outFile_name = '../txt/sistematici_'+varLabel+ptBin_label+'.txt'

    chib1_eff = Efficiency()
    chib1_eff = chib1_eff.loadFromFile(Chib1_eff_file)
        
    chib2_eff = Efficiency()
    chib2_eff = chib2_eff.loadFromFile(Chib2_eff_file)
        
    numChib = NumChib()
    numChib = numChib.loadFromFile(numChib_file)

    value = numChib.ratio_21
    statistic = numChib.s_ratio_21/numChib.ratio_21

    ratioEff = sqrt( (chib1_eff.s_eff/chib1_eff.eff)**2 +  (chib2_eff.s_eff/chib2_eff.eff)**2 )

    sistematici = Sistematici(value = value, statistic = statistic, ratioEff = ratioEff)
    sistematici.saveToFile(outFile_name)


def makeTable(varLabel,printAllData=False):

    tipiSistematici = dict(
        ratioEff = 'Statistical error on $\\epsilon_1/\\epsilon_2$',
        dscbParam = 'Signal parameterization',
        funzParam = 'Data fitting systematic',
        ptDistro = 'Choice of $\\chi_b$ $p_T$ spectrum',
        )
    
    sistematici = [0]*(len(ptBins)+1) # I leave place 0 empty
    for ptBin in range(1, len(ptBins)):
        ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
        sistematici[ptBin] = Sistematici()
        sistematici[ptBin] = sistematici[ptBin].loadFromFile('../txt/sistematici_'+varLabel+ptBin_label+'.txt')
    
    sistematici_allData = Sistematici()
    sistematici_allData = sistematici_allData.loadFromFile('../txt/sistematici_'+varLabel+'.txt')
    
    with open('../tables/sistematici_'+varLabel+'.tex','w') as outFile:
        numColumns = len(ptBins)
        if printAllData:
            numColumns += 1
        outFile.write('\\begin{center} \n \\begin{tabular}{p{0.43\linewidth}*{'+str(numColumns-1)+'}{c}}\n\\toprule\n')
        outFile.write('$p_T(\\Upsilon)$ [GeV] ')
        for ptBin in range(1, len(ptBins)):
            outFile.write(' & '+str(ptBins[ptBin-1])+' - '+str(ptBins[ptBin])+' ')
        if printAllData:
            outFile.write(' & '+str(ptBins[0])+' - '+str(ptBins[-1])+' ')
        outFile.write('\\\\ \n \\midrule \n')
        outFile.write('Source of uncertainty  &  \\multicolumn{'+str(numColumns-1)+'}{c}{Systematic uncertainty (\\%)} \\\\ \n \\midrule \n')
        for key in tipiSistematici.keys():
            outFile.write(tipiSistematici[key])
            for ptBin in range(1, len(ptBins)):
                outFile.write(' & '+'{:.1f}'.format(getattr(sistematici[ptBin], key)*100))
            if printAllData:
                outFile.write(' & '+'{:.1f}'.format(getattr(sistematici_allData, key)*100))
            outFile.write('\\\\ \n')
            
        outFile.write('\\midrule \n')
        outFile.write('Total uncertainty')
        for ptBin in range(1, len(ptBins)):
            outFile.write(' & '+'{:.1f}'.format(sistematici[ptBin].s_rel_total()* 100))
        if printAllData:
            outFile.write(' & '+'{:.1f}'.format(sistematici_allData.s_rel_total()*100))
        outFile.write('\\\\ \n \\midrule \n')
        outFile.write('Statistical error')
        for ptBin in range(1, len(ptBins)):
            outFile.write(' & '+'{:.1f}'.format(sistematici[ptBin].statistic* 100))
        if printAllData:
            outFile.write(' & '+'{:.1f}'.format(sistematici_allData.statistic*100))
        outFile.write('\\\\ \n')
        outFile.write('\\bottomrule \n \\end{tabular} \n \\end{center}')



if __name__ == '__main__':

    varLabel = 'qValue'
    ptBin = None
    isMakeTable = False
    isSetup = False
    ptSpectrum_file = '2S'
    try:
        opts, args = getopt.getopt(sys.argv[1:],'fb:',['refit', 'ptBin=', 'makeTable', 'setup', 'filePt='])
    except getopt.GetoptError:
        print './sistematici.py [-f] [--makeTable] [--setup] [--filePt <1S, 2S, 3S>]'
        print '-f = --refit'   
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-f', '--refit'):
            varLabel = 'refit'
        if opt in ('-b', '--ptBin'):
            ptBin = int(arg)
        if opt == '--makeTable':
            isMakeTable = True
        if opt == '--setup':
            isSetup = True  
        if opt == '--filePt':
            ptSpectrum_file = arg

    if ptSpectrum_file in ('1S', '1s'):
        ptFileLabel = 'UpsilonPt'
        ptSpectrum_file = '1S'
    elif ptSpectrum_file in ('2S', '2s'):
        ptFileLabel = 'Upsilon2SPt'
        ptSpectrum_file = '2S'
    elif ptSpectrum_file in ('3S', '3s'):
        ptFileLabel = 'Upsilon3SPt'
        ptSpectrum_file = '3S'
    
            
 
    if isMakeTable:
        for varLabel in ('refit',):# 'qValue'):
            makeTable(varLabel)

    if isSetup:
        for varLabel in ('refit', ):#'qValue'):
            for ptBin in [None] + range(1,len(ptBins)):
                
                if ptBin == None:
                    ptBin_label = ''
                    ptMin = ptBins[0]
                    ptMax = ptBins[-1]
                else:
                    ptBin_label = '_'+str(ptBins[ptBin-1])+'_'+str(ptBins[ptBin])
                    ptMin = ptBins[ptBin-1]
                    ptMax = ptBins[ptBin]
                    
                setupFileSistematici(varLabel, ptMin, ptMax, ptFileLabel, ptSpectrum_file)
