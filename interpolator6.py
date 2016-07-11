#!/usr/bin/env python
import argparse,ROOT,sys,pprint,array
from pointDict import pointDict
from massDict import massDict
import subprocess
class interpolator:
    def __init__(self,fileName):
        #Number of events needed for 3sigma observation
        #self.contours =[ 50.9, 40.9, 30.5, 24.2, 19.6, 85.6, 66.0, 48.0, 37.6, 30.4, 11.5, 9.4, 7.5, 6.1, 4.4, 17.3, 13.9, 11.5, 9.5, 7.7]
        #number of events needed for 95%CL exclusion
        self.contours =[48.8,39.0,30.1,24.2,19.7,
                        98.9,69.9,51.5,38.1,29.1,
                        12.7,11.2,9.7,8.8,7.5,
                        17.5,14.8,12.8,11.2,9.9]
        self.lines = []
        self.intersections = []
        self.rootFile = ROOT.TFile.Open(fileName)
        if not self.rootFile:
            print 'File not found',fileName
            sys.exit(1)
        self.lumi = 5800
        self.outFile = ROOT.TFile.Open('output_rpv6.root','RECREATE')
        self.effHists = [ROOT.TH1F('h_eff_SR'+str(i),'h_eff_SR'+str(i),57,587.5,2012.5) for i in range(1,21)]
        self.yieldHists = []
        self.h0 = [ROOT.TH1F('h0_'+str(i),'h0_'+str(i),57,587.5,2012.5) for i in range(1,21)]
        self.e0 = [ROOT.TH1F('e0_'+str(i),'e0_'+str(i),57,587.5,2012.5) for i in range(1,21)]
        for key in self.rootFile.GetListOfKeys():
            if 'h_SRyield_' in key.GetName() and 'unweighted' not in key.GetName():
                sList = key.GetName().split('_')
                dsid = int(sList[-1])
                print dsid
                h = self.rootFile.Get(key.GetName())
                cutflow = self.rootFile.Get('h_cutflow_'+str(dsid))
                mG = pointDict[dsid][0]
                for i in range(20):
                    self.effHists[i].Fill(mG,h.GetBinContent(i+1)/cutflow.GetBinContent(2))
                    self.e0[i].Fill(mG,h.GetBinContent(i+1)/cutflow.GetBinContent(2))
                    self.h0[i].Fill(mG,h.GetBinContent(i+1))
        for effHist in self.effHists:
            effHist.GetXaxis().SetTitle('m_{#tilde{g}} [GeV]')
            effHist.GetYaxis().SetTitle('yield in signal region')
            effHist.SetMarkerSize(1.0)
            effHist.GetYaxis().SetTitleOffset(1.5)
            self.fillSpaces(effHist)
        self.makeYieldHists()
        self.makeLines()

    def interpolate(self,m1,m2,m3,e1,e3):
        return e1+(e3-e1)*(m2-m1)/(m3-m1)

    def loopAndFill(self,h,pList,goalPrec):
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
            if h.GetXaxis().GetBinCenter(i) in intDict and h.GetBinContent(i) == 0:
                h.SetBinContent(i,intDict[h.GetXaxis().GetBinCenter(i)])

    def fillSpaces(self,h):
        pList = []
        for i in range(1,h.GetNbinsX()+1):
            if h.GetBinContent(i) != 0:
                pList.append((h.GetXaxis().GetBinCenter(i),h.GetBinContent(i)))
            self.loopAndFill(h,pList,100)
            self.loopAndFill(h,pList,50)
            self.loopAndFill(h,pList,25)

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
                mass = h.GetXaxis().GetBinCenter(i)
                if mass in massDict:
                    h.SetBinContent(i,self.effHists[k].GetBinContent(i)*self.lumi*massDict[mass])
            self.yieldHists[k].SetMaximum(self.yieldHists[k].GetMaximum()*1.5)
            self.yieldHists[k].SetMinimum(0)
    def makeLines(self):
        for i in range(len(self.yieldHists)):
            h = self.yieldHists[i]
            c = self.contours[i]
            foundLine = False
            for j in range(1,h.GetNbinsX()):
                if h.GetBinContent(j) > c and h.GetBinContent(j+1) < c:
                    foundLine = True
                    m1 = h.GetBinCenter(j)
                    m3 = h.GetBinCenter(j+1)
                    y1 = h.GetBinContent(j)
                    y2 = c
                    y3 = h.GetBinContent(j+1)
                    m2 = (y2-y1)*(m3-m1)/(y3-y1)+m1
                    self.intersections.append((m2,h.Interpolate(m2)))
                    self.lines.append(ROOT.TLine(m2,0,m2,10))
            if not foundLine:
                self.lines.append(ROOT.TLine(2000,0,2000,10))
                self.intersections.append((0,0))
                
parser = argparse.ArgumentParser(add_help=False, description='Plot Yields')
parser.add_argument('input')
args = parser.parse_args()
foo = interpolator(args.input)    
foo.write()

ROOT.gROOT.LoadMacro('~/atlasstyle/AtlasStyle.C')
ROOT.gROOT.LoadMacro('~/atlasstyle/AtlasLabels.C')
ROOT.SetAtlasStyle()
ROOT.gStyle.SetPaintTextFormat('2.3f')

cLines = [ROOT.TLine(900,c,1400,c) for c in foo.contours]
cols = [ROOT.kRed,ROOT.kGreen,ROOT.kBlue,ROOT.kViolet,ROOT.kCyan]*4
labs = ['MJ > 600 GeV','MJ > 650 GeV','MJ > 700 GeV','MJ > 750 GeV','MJ > 800 GeV']
catlab = ['#splitline{n_{jet}#geq 4}{b-tag}',
          '#splitline{n_{jet}#geq 4}{b-inclusive}',
          '#splitline{n_{jet}#geq 5}{b-tag}',
          '#splitline{n_{jet}#geq 5}{b-inclusive}']

leg = ROOT.TLegend(0.65,0.65,0.8,0.9)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.03)
lumiLatex=ROOT.TLatex()

c = [ROOT.TCanvas('c'+str(i+1),'c'+str(i+1)) for i in range(4)]
j=0
for i in range(len(foo.yieldHists)):
    if i%5==0:
        c[j].cd()
        j+=1
    if i < 5:
        leg.AddEntry(foo.yieldHists[i],labs[i],'l')
    foo.yieldHists[i].SetLineColor(cols[i])
    foo.yieldHists[i].SetLineWidth(2)
    foo.yieldHists[i].SetLineStyle(ROOT.kDashed)
    if i%5 == 0:
        foo.yieldHists[i].Draw('c hist')
    else:
        foo.yieldHists[i].Draw('c hist same')
    foo.yieldHists[i].GetXaxis().SetRangeUser(900,1400)
    cLines[i].SetLineColor(cols[i])
    cLines[i].SetLineWidth(2)
    if foo.intersections[i][0] != 0:
        cLines[i].DrawLine(900,foo.intersections[i][1],foo.intersections[i][0],foo.intersections[i][1])
        cLines[i].DrawLine(foo.intersections[i][0],foo.yieldHists[i].GetMinimum(),foo.intersections[i][0],foo.intersections[i][1])
    else:
        cLines[i].DrawLine(900,foo.contours[i],950,foo.contours[i])

    if i%5 ==0:
        ROOT.ATLASLabel(0.2,0.85,'Internal')
        lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 5.8 fb^{-1}')
        lumiLatex.DrawLatexNDC(0.75,0.5,catlab[j-1])
        leg.Draw()

for suffix in ['.pdf','.png','.C']:
    c[0].Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/07_07_5p8fb/reach_RPV6_m4_b1_95CL'+suffix)
    c[1].Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/07_07_5p8fb/reach_RPV6_m4_b9_95CL'+suffix)
    c[2].Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/07_07_5p8fb/reach_RPV6_m5_b1_95CL'+suffix)
    c[3].Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/07_07_5p8fb/reach_RPV6_m5_b9_95CL'+suffix)
subprocess.call('chmod a+r /global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/07_07_5p8fb/*',shell=True)
