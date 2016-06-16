#!/usr/bin/env python
import argparse,ROOT,sys,pprint
from pointDict import pointDict
from massDict import massDict
class interpolator:
    def __init__(self,fileName,histName=''):
        self.nCuts = 10
        self.rootFile = ROOT.TFile.Open(fileName)
        self.goalPrec = 50
        self.lumi = 6E3
        if not self.rootFile:
            print 'File not found',fileName
            sys.exit(1)

        self.outFile = ROOT.TFile.Open('output.root','RECREATE')
        self.effHists = []
        if not histName:
            self.effHists = [ROOT.TH2F('h_'+str(i),'h_'+str(i),57,587.5,2012.5,9,-50,1750) for i in range(self.nCuts)]
            for key in self.rootFile.GetListOfKeys():
                if 'h_cutflow' in key.GetName():
                    dsid = int(key.GetName().split('_')[2])
                    h = self.rootFile.Get(key.GetName())
                    if h.GetNbinsX() != self.nCuts:
                        print 'Warning: wrong number of cuts!'
                        self.nCuts = h.GetNbinsX()
                    mG = pointDict[dsid][0]
                    mX = pointDict[dsid][1]
                    for i in range(self.nCuts):
                        self.effHists[i].Fill(mG,mX,h.GetBinContent(i+1)/h.GetBinContent(1))
        else:
            self.nCuts=1
            self.effHists=[self.rootFile.Get(histName)]
        for hist in self.effHists:
            hist.GetXaxis().SetTitle('m_{#tilde{g}} [GeV]')
            hist.GetYaxis().SetTitle('m_{#tilde{#chi}} [GeV]')
            hist.SetMarkerSize(1.5)
            hist.GetYaxis().SetTitleOffset(1.5)
            
    def interpolate(self,m1,m2,m3,e1,e3):
        return e1+(e3-e1)*(m2-m1)/(m3-m1)
    def loopAndFill(self,pList,goalPrec):
        while(self.getPrec(pList) > goalPrec):
            prec = self.getPrec(pList)
            for mi in range(len(pList)-1):
                m1 = pList[mi][0]
                m3 = pList[mi+1][0]
                m2 = (m1+m3)/2.0
                e1 = pList[mi][1]
                e3 = pList[mi+1][1]
                e2=self.interpolate(m1,m2,m3,e1,e3)

                if m3-m1 > goalPrec:
                    pList.insert(mi+1,(m2,e2))
    def fillSpaces(self):
        h= self.effHists[self.nCuts-1]
        for j in range(1,h.GetNbinsY()+1):
            pList = []
            for i in range(1,h.GetNbinsX()+1):
                if h.GetBinContent(i,j) != 0:
                    pList.append((h.GetXaxis().GetBinCenter(i),h.GetBinContent(i,j)))
            self.loopAndFill(pList,100)
            self.loopAndFill(pList,50)
#            self.loopAndFill(pList,25)
            for i in range(1,h.GetNbinsX()+1):
                intDict = dict(pList)# = [item[1] for item in pList if item[0] == h.GetBinCenter(i,j)][0]
                if h.GetXaxis().GetBinCenter(i) in intDict and h.GetBinContent(i,j) == 0:
                    h.SetBinContent(i,j,intDict[h.GetXaxis().GetBinCenter(i)])
    def write(self):
        self.outFile.Write()
    def getPrec(self,someList):
        prec = 0
        for i in range(len(someList)-1):
            if someList[i+1][0]-someList[i][0] > prec:
                prec = someList[i+1][0]-someList[i][0]
        return prec
    def makeYieldHist(self):

        self.yieldHist = self.effHists[self.nCuts-1].Clone('yieldhist')
        h= self.yieldHist

        for i in range(1,h.GetNbinsX()+1):
            for j in range(1,h.GetNbinsY()+1):
                mass = h.GetXaxis().GetBinCenter(i)
                if mass in massDict:
                    h.SetBinContent(i,j,self.effHists[self.nCuts-1].GetBinContent(i,j)*self.lumi*massDict[mass])

parser = argparse.ArgumentParser(add_help=False, description='Plot Yields')
parser.add_argument('input')
parser.add_argument('--model',dest='model',type=str,default='RPV10')
args = parser.parse_args()

foo = interpolator(args.input)    

foo.fillSpaces()
foo.makeYieldHist()
foo.write()
ROOT.gROOT.LoadMacro('~/atlasstyle/AtlasStyle.C')
ROOT.gROOT.LoadMacro('~/atlasstyle/AtlasLabels.C')
ROOT.SetAtlasStyle()
ROOT.gStyle.SetPaintTextFormat('2.2f')
c1=ROOT.TCanvas('c1','c1',800,600)
foo.yieldHist.Draw('text45')
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex=ROOT.TLatex()
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 6.0 fb^{-1}')
c2=ROOT.TCanvas('c2','c2',800,600)
foo.effHists[foo.nCuts-1].Draw('text45')
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 6.0 fb^{-1}')
