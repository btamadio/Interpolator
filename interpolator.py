#!/usr/bin/env python
import argparse,ROOT,sys,pprint,array
from pointDict import pointDict
from massDict import massDict
class interpolator:
    def __init__(self,fileName):
        #Number of events needed for 3sigma observation
        self.contours = [43.1,30.3,23.9,18.2,13.8,68.0,48.9,36.0,27.5,20.5,10.0,7.5,6.1,4.4,3.7,14.9,11.0,9.1,7.1,5.7]
        self.rootFile = ROOT.TFile.Open(fileName)
        if not self.rootFile:
            print 'File not found',fileName
            sys.exit(1)
        self.lumi = 4800
        self.outFile = ROOT.TFile.Open('output.root','RECREATE')
        self.effHists = [ROOT.TH2F('h_eff_SR'+str(i),'h_eff_SR'+str(i),57,587.5,2012.5,9,-50,1750) for i in range(1,21)]
        self.yieldHists = []
        self.h0 = [ROOT.TH2F('h0_'+str(i),'h0_'+str(i),57,587.5,2012.5,9,-50,1750) for i in range(1,21)]
        self.e0 = [ROOT.TH2F('e0_'+str(i),'e0_'+str(i),57,587.5,2012.5,9,-50,1750) for i in range(1,21)]
        for key in self.rootFile.GetListOfKeys():
            if 'h_sigyield_dy_' in key.GetName():
                sList = key.GetName().split('_')
                dsid = int(sList[-1])
                h = self.rootFile.Get(key.GetName())
                cutflow = self.rootFile.Get('h_cutflow_'+str(dsid))
                mG = pointDict[dsid][0]
                mX = pointDict[dsid][1]
                for i in range(20):
                    self.effHists[i].Fill(mG,mX,h.GetBinContent(i+1)/cutflow.GetBinContent(1))
                    self.e0[i].Fill(mG,mX,h.GetBinContent(i+1)/cutflow.GetBinContent(1))
                    self.h0[i].Fill(mG,mX,h.GetBinContent(i+1))
        for effHist in self.effHists:
            effHist.GetXaxis().SetTitle('m_{#tilde{g}} [GeV]')
            effHist.GetYaxis().SetTitle('m_{#tilde{#chi}} [GeV]')
            effHist.SetMarkerSize(1.0)
            effHist.GetYaxis().SetTitleOffset(1.5)
            self.fillSpaces(effHist)
        self.makeYieldHists()

    def interpolate(self,m1,m2,m3,e1,e3):
        return e1+(e3-e1)*(m2-m1)/(m3-m1)

    def loopAndFill(self,h,j,pList,goalPrec):
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
        for i in range(1,h.GetNbinsX()+1):
            intDict = dict(pList)
            if h.GetXaxis().GetBinCenter(i) in intDict and h.GetBinContent(i,j) == 0:
                h.SetBinContent(i,j,intDict[h.GetXaxis().GetBinCenter(i)])

    def fillSpaces(self,h):
        for j in range(1,h.GetNbinsY()+1):
            pList = []
            for i in range(1,h.GetNbinsX()+1):
                if h.GetBinContent(i,j) != 0:
                    pList.append((h.GetXaxis().GetBinCenter(i),h.GetBinContent(i,j)))
            self.loopAndFill(h,j,pList,100)
            self.loopAndFill(h,j,pList,50)
            self.loopAndFill(h,j,pList,25)

            if j == 6:
                for binx in range(27,33):
#                    print h.GetXaxis().GetBinCenter(binx)
#                    print h.GetYaxis().GetBinCenter(j)
                    m1 = 850.0
                    m2 = 1050.0
                    m3 = 1250.0
                    e0 = h.GetBinContent(binx,4)
                    e1 = h.GetBinContent(binx,5)
                    e2 = h.GetBinContent(binx,6)
                    e3 = h.GetBinContent(binx,7)
                    e3 = (m3-m2)*(e2-e1)/(m2-m1)+e2
#                    print 'mG,e0,e1,e2,e3=',h.GetXaxis().GetBinCenter(binx),e0,e1,e2,e3            
#                    if h.GetXaxis().GetBinCenter(binx) == 1300:
                        #h.SetBinContent(binx,7,e3)

                
    def write(self):
        self.outFile.Write()

    def getPrec(self,someList):
        prec = 0
        for i in range(len(someList)-1):
            if someList[i+1][0]-someList[i][0] > prec:
                prec = someList[i+1][0]-someList[i][0]
        return prec

    def makeYieldHists(self):
        for i in range(len(self.effHists)):
            self.yieldHists.append(self.effHists[i].Clone('yieldhist_'+str(i)))
        for k in range(len(self.yieldHists)):
            h = self.yieldHists[k]
            for i in range(1,h.GetNbinsX()+1):
                for j in range(1,h.GetNbinsY()+1):
                    mass = h.GetXaxis().GetBinCenter(i)
                    if mass in massDict:
                        h.SetBinContent(i,j,self.effHists[k].GetBinContent(i,j)*self.lumi*massDict[mass])

            h.SetContour(1,array.array('d',[self.contours[k]]))
parser = argparse.ArgumentParser(add_help=False, description='Plot Yields')
parser.add_argument('input')
parser.add_argument('--model',dest='model',type=str,default='RPV10')
args = parser.parse_args()
foo = interpolator(args.input)    
foo.write()

ROOT.gROOT.LoadMacro('~/atlasstyle/AtlasStyle.C')
ROOT.gROOT.LoadMacro('~/atlasstyle/AtlasLabels.C')
ROOT.SetAtlasStyle()
ROOT.gStyle.SetPaintTextFormat('2.3f')

hists_m4_b1 = foo.yieldHists[0:5]
hists_m4_b9 = foo.yieldHists[5:10]
hists_m5_b1 = foo.yieldHists[10:15]
hists_m5_b9 = foo.yieldHists[15:20]
cols = [ROOT.kRed,ROOT.kGreen,ROOT.kBlue,ROOT.kViolet,ROOT.kCyan]
labs = ['MJ > 600 GeV','MJ > 650 GeV','MJ > 700 GeV','MJ > 750 GeV','MJ > 800 GeV']
leg = ROOT.TLegend(0.75,0.65,0.8,0.9)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.03)
lumiLatex=ROOT.TLatex()

c1=ROOT.TCanvas('c1','c1',800,600)
for i in range(len(hists_m4_b1)):
    leg.AddEntry(hists_m4_b1[i],labs[i],'l')
    hists_m4_b1[i].SetLineColor(cols[i])
    hists_m4_b1[i].SetLineWidth(2)
    if i == 0:
        hists_m4_b1[i].Draw('cont3')
    else:
        hists_m4_b1[i].Draw('cont3 same')
    hists_m4_b1[i].GetXaxis().SetRangeUser(1200,1800)
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 4.8 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 4}{b-tag}')
leg.Draw()

c2=ROOT.TCanvas('c2','c2',800,600)
for i in range(len(hists_m4_b9)):
    hists_m4_b9[i].SetLineColor(cols[i])
    hists_m4_b9[i].SetLineWidth(2)
    if i == 0:
        hists_m4_b9[i].Draw('cont3')
    else:
        hists_m4_b9[i].Draw('cont3 same')
    hists_m4_b9[i].GetXaxis().SetRangeUser(1200,1800)
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 4.8 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 4}{b-inclusive}')
leg.Draw()

c3=ROOT.TCanvas('c3','c3',800,600)
for i in range(len(hists_m5_b1)):
    hists_m5_b1[i].SetLineColor(cols[i])
    hists_m5_b1[i].SetLineWidth(2)
    if i == 0:
        hists_m5_b1[i].Draw('cont3')
    else:
        hists_m5_b1[i].Draw('cont3 same')
    hists_m5_b1[i].GetXaxis().SetRangeUser(1200,1800)
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 4.8 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 5}{b-tag}')
leg.Draw()

c4=ROOT.TCanvas('c4','c4',800,600)
for i in range(len(hists_m5_b9)):
    hists_m5_b9[i].SetLineColor(cols[i])
    hists_m5_b9[i].SetLineWidth(2)
    if i == 0:
        hists_m5_b9[i].Draw('cont3')
    else:
        hists_m5_b9[i].Draw('cont3 same')
    hists_m5_b9[i].GetXaxis().SetRangeUser(1200,1800)
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 4.8 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 5}{b-inclusive}')
leg.Draw()
