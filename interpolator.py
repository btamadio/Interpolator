#!/usr/bin/env python
import argparse,ROOT,sys,pprint,array
from pointDict import pointDict
from massDict import massDict
class interpolator:
    def __init__(self,fileName):
        #Number of events needed for 3sigma observation
        #self.contours = [43.1,30.3,23.9,18.2,13.8,68.0,48.9,36.0,27.5,20.5,10.0,7.5,6.1,4.4,3.7,14.9,11.0,9.1,7.1,5.7]
        #for RPV10, optimal SR is njet >=5, b-tag, MJ > 800
        #this correspond to self.contours[14] = 2.8

        #95% CL
        self.contours = [63.1,42.2,30.1,22.4,17.6,129.9,83.2,57.9,41.1,29.4,11.2,10.0,8.6,7.5,6.7,16.3,13.9,12.0,10.1,9.0]
        #3sigma
        #self.contours = [74.4,48.5,33.4,23.8,17.8,157.2,99.2,68.0,47.1,32.6,9.6,8.0,6.0,4.2,2.8,16.1,13.0,10.6,8.1,6.6]

        self.rootFile = ROOT.TFile.Open(fileName)
        if not self.rootFile:
            print 'File not found',fileName
            sys.exit(1)
        self.lumi = 5800
        self.outFile = ROOT.TFile.Open('output.root','RECREATE')
        self.effHists = [ROOT.TH2F('h_eff_SR'+str(i),'h_eff_SR'+str(i),57,587.5,2012.5,9,-50,1750) for i in range(1,21)]
        self.yieldHists = []
        self.yieldGraphs = []
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
                    self.effHists[i].Fill(mG,mX,h.GetBinContent(i+1)/cutflow.GetBinContent(2))
                    self.e0[i].Fill(mG,mX,h.GetBinContent(i+1)/cutflow.GetBinContent(2))
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
        yieldCan = ROOT.TCanvas('yieldCan','yieldCan',800,600)
        for k in range(len(self.yieldHists)):
            h = self.yieldHists[k]
            for i in range(1,h.GetNbinsX()+1):
                for j in range(1,h.GetNbinsY()+1):
                    mass = h.GetXaxis().GetBinCenter(i)
                    if mass in massDict:
                        h.SetBinContent(i,j,self.effHists[k].GetBinContent(i,j)*self.lumi*massDict[mass])

            h.SetContour(1,array.array('d',[self.contours[k]]))
            h.Draw('cont z list')
            yieldCan.Update()

            contourTList = ROOT.gROOT.GetListOfSpecials().FindObject('contours')
            for cont in contourTList:
                for c in cont:
                    self.yieldGraphs.append(c.Clone('yieldGraph_'+str(k)))
            self.yieldGraphs[k] = self.trimGraph(self.yieldGraphs[k])
            print self.yieldGraphs[k]
    def trimGraph(self,g):
        xarr = []
        yarr = []
        for iPoint in range(g.GetN()):
            x,y=ROOT.Double(0),ROOT.Double(0)
            g.GetPoint(iPoint,x,y)
            if x < 1200 or (x < 1400 and y > 1000) or (x < 1600 and y >1200) or (x < 1800 and y > 1400):
                pass
            else:
                xarr.append(x)
                yarr.append(y)
        if len(xarr) > 0:
            return ROOT.TGraph(len(xarr),array.array('d',xarr),array.array('d',yarr))
        else:
            return 0
    def getMaxXY(self,g):
        result = (0,0)
        yMax = 0
        xMax = 0
        yAtMaxX = 0
        xMin = 999999
        for iPoint in range(g.GetN()):
            x,y=ROOT.Double(0),ROOT.Double(0)
            g.GetPoint(iPoint,x,y)
            if x > xMax:
                xMax = x
                yAtMaxX = y
        for iPoint in range(g.GetN()):
            x,y=ROOT.Double(0),ROOT.Double(0)
            g.GetPoint(iPoint,x,y)
            if y > yAtMaxX and x < xMin:
                xMin = x
                result = x,y
        return result

    def makeLimitHist(self,index,uncert=0.1):
        self.limitHist_1sig = self.yieldHists[index].Clone('limit_1sig_'+str(index))
        self.limitHist_2sig = self.yieldHists[index].Clone('limit_2sig_'+str(index))
        centVal = self.contours[index]
        uncert*=centVal
        cList1 = [centVal-uncert,centVal,centVal+uncert]
        cList2 = [centVal-2*uncert,centVal,centVal+2*uncert]
        self.limitHist_1sig.SetContour(len(cList1),array.array('d',cList1))
        self.limitHist_2sig.SetContour(len(cList2),array.array('d',cList2))
        self.canvas=ROOT.TCanvas('canvas','canvas',800,600)

        self.limitHist_1sig.Draw('cont z list')
        self.canvas.Update()
        self.limitConts_1sig = []
        self.limitConts_2sig = []

        contourTList = ROOT.gROOT.GetListOfSpecials().FindObject('contours')
        ci = 0
        for cont in contourTList:
            for c in cont:
                self.limitConts_1sig.append(c.Clone('cont_1sig_'+str(ci)))
                ci+=1
        
        self.canvas.cd()
        self.limitHist_2sig.Draw('cont z list')
        self.canvas.Update()
        contourTList = ROOT.gROOT.GetListOfSpecials().FindObject('contours')
        ci = 0
        for cont in contourTList:
            for c in cont:
                self.limitConts_2sig.append(c.Clone('cont_2sig_'+str(ci)))
                ci+=1
        self.limit_nom = self.trimGraph(self.limitConts_1sig[1])
        self.limit_1up = self.trimGraph(self.limitConts_1sig[0])
        self.limit_1down = self.trimGraph(self.limitConts_1sig[2])
        self.limit_2up = self.trimGraph(self.limitConts_2sig[0])
        self.limit_2down = self.trimGraph(self.limitConts_2sig[2])
        self.canvas.Clear()
        self.limit_2up.SetFillColor(ROOT.kGreen)
        self.limit_1up.SetFillColor(ROOT.kYellow)
        self.limit_nom.SetFillColor(ROOT.kYellow)
        self.limit_nom.SetLineColor(ROOT.kBlack)
        self.limit_1down.SetFillColor(ROOT.kGreen)
        self.limit_1down.SetLineColor(ROOT.kBlack)
        self.limit_2down.SetFillColor(ROOT.kBlack)
        
        maxXY_2up = self.getMaxXY(self.limit_2up)
        maxXY_1down = self.getMaxXY(self.limit_1down)

#        self.limit_2down.SetPoint(self.limit_2down.GetN(),maxXY_1down[0],maxXY_1down[1])
#        self.limit_2down.SetPoint(self.limit_2down.GetN(),maxXY_2up[0],maxXY_2up[1])
#        self.limit_2down.SetPoint(self.limit_2down.GetN(),850,maxXY_2up[1])
        
parser = argparse.ArgumentParser(add_help=False, description='Plot Yields')
parser.add_argument('input')
parser.add_argument('--model',dest='model',type=str,default='RPV10')
args = parser.parse_args()
foo = interpolator(args.input)    
foo.makeLimitHist(14,0.05)
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
    hists_m4_b1[i].GetXaxis().SetRangeUser(1200,2000)
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 5.8 fb^{-1}')
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
    hists_m4_b9[i].GetXaxis().SetRangeUser(1200,2000)
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 5.8 fb^{-1}')
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
    hists_m5_b1[i].GetXaxis().SetRangeUser(1200,2000)
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 5.8 fb^{-1}')
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
    hists_m5_b9[i].GetXaxis().SetRangeUser(1200,2000)
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 5.8 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 5}{b-inclusive}')
leg.Draw()

c5 = ROOT.TCanvas('c5','c5',800,600)
foo.limit_2up.Draw('AFL')
c5.Modified()
#foo.limit_2up.GetXaxis().SetRangeUser(1200,2000)
#c5.Update()
foo.limit_1up.Draw('FL')
foo.limit_nom.Draw('FL')

foo.limit_1down.Draw('FL')
foo.limit_nom.Draw('L')
foo.limit_2down.Draw('FL')


foo.limit_2up.GetXaxis().SetLimits(800,1800)
foo.limit_2up.GetXaxis().SetTitle('m_{#tilde{g}} [GeV]')
foo.limit_2up.GetYaxis().SetTitle('m_{#tilde{#chi}} [GeV]')

ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 5.8 fb^{-1}')
lumiLatex.DrawLatexNDC(0.7,0.3,'#splitline{#splitline{n_{jet}#geq 5}{b-tag}}{MJ > 800 GeV}')
