#!/usr/bin/env python
import argparse,ROOT,sys,pprint,array
from pointDict import pointDict
from massDict import massDict
class interpolator:
    def __init__(self,fileName):
        #Number of events needed for 3sigma observation
        #self.contours = [43.1,30.3,23.9,18.2,13.8,68.0,48.9,36.0,27.5,20.5,10.0,7.5,6.1,4.4,3.7,14.9,11.0,9.1,7.1,5.7]
        #RPV6: optimized SR is nJet >= 5, b-tag, MJ > 600
        #this corresponds to self.contours[10] = 9.6
        self.contours = [74.4,48.5,33.4,23.8,17.8,157.2,99.2,68.0,47.1,32.6,9.6,8.0,6.0,4.2,2.8,16.1,13.0,10.6,8.1,6.6]
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
            if 'h_sigyield_dy_' in key.GetName():
                sList = key.GetName().split('_')
                dsid = int(sList[-1])
                h = self.rootFile.Get(key.GetName())
                cutflow = self.rootFile.Get('h_cutflow_'+str(dsid))
                mG = pointDict[dsid][0]
                for i in range(20):
                    self.effHists[i].Fill(mG,h.GetBinContent(i+1)/cutflow.GetBinContent(2))
                    self.e0[i].Fill(mG,h.GetBinContent(i+1)/cutflow.GetBinContent(1))
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
parser.add_argument('--model',dest='model',type=str,default='RPV6')
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


#c[0].Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/reach_RPV6_m4_b1_dy14_5p8fb.pdf')
#c[1].Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/reach_RPV6_m4_b9_dy14_5p8fb.pdf')
#c[2].Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/reach_RPV6_m5_b1_dy14_5p8fb.pdf')
#c[3].Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/reach_RPV6_m5_b9_dy14_5p8fb.pdf')
