#!/usr/bin/env python
import argparse,ROOT,sys,pprint,array
from pointDict import pointDict
from massDict import massDict
class interpolator:
    def __init__(self,fileName):
        #3 sigma
#        self.contours =[ 50.9, 40.9, 30.5, 24.2, 19.6, 85.6, 66.0, 48.0, 37.6, 30.4, 11.5, 9.4, 7.5, 6.1, 4.4, 17.3, 13.9, 11.5, 9.5, 7.7]

        #95% CL
        self.contours =[48.8,39.0,30.1,24.2,19.7,
                        98.9,69.9,51.5,38.1,29.1,
                        12.7,11.2,9.7,8.8,7.5,
                        17.5,14.8,12.8,11.2,9.9]

        self.rootFile = ROOT.TFile.Open(fileName)
        if not self.rootFile:
            print 'File not found',fileName
            sys.exit(1)
        self.lumi = 5800
#        self.lumi = 7140
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
    def trimGraph(self,g,p=False):
        pointArr = []
        #xarr = []
        #yarr = []
        for iPoint in range(g.GetN()):
            x,y=ROOT.Double(0),ROOT.Double(0)
            g.GetPoint(iPoint,x,y)
            if x < 1000 or (x < 1200 and y > 850) or (x < 1400 and y > 1050) or (x < 1600 and y >1250) or (x < 1800 and y > 1450):
                pass
            else:
                pointArr.append((x,y))
                if p:
                    print iPoint,x,y
                #xarr.append(x)
                #yarr.append(y)
        if len(pointArr) > 0:
            #pointArr.sort(key=lambda x:x[0])
            xarr = [ p[0] for p in pointArr ]
            yarr = [ p[1] for p in pointArr ]
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

    def makeLimitHist(self,index,uncert=0.3):
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
        self.limit_nom = self.trimGraph(self.limitConts_1sig[1],True)
        self.limit_1up = self.trimGraph(self.limitConts_1sig[0])
        self.limit_1down = self.trimGraph(self.limitConts_1sig[2])
        self.limit_2up = self.trimGraph(self.limitConts_2sig[0])
        self.limit_2down = self.trimGraph(self.limitConts_2sig[2])
        self.canvas.Clear()
        self.limit_2up.SetFillColor(ROOT.kGreen)

parser = argparse.ArgumentParser(add_help=False, description='Plot Yields')
parser.add_argument('input')
args = parser.parse_args()
foo = interpolator(args.input)    
foo.makeLimitHist(14,0.15)
foo.write()

ROOT.gROOT.LoadMacro('~/atlasstyle/AtlasStyle.C')
ROOT.gROOT.LoadMacro('~/atlasstyle/AtlasLabels.C')
ROOT.SetAtlasStyle()
ROOT.gStyle.SetPaintTextFormat('2.3f')

#hists_m4_b1 = foo.yieldHists[0:5]
#hists_m4_b9 = foo.yieldHists[5:10]
#hists_m5_b1 = foo.yieldHists[10:15]
#hists_m5_b9 = foo.yieldHists[15:20]

hists_m4_b1 = foo.yieldGraphs[0:5]
hists_m4_b9 = foo.yieldGraphs[5:10]
hists_m5_b1 = foo.yieldGraphs[10:15]
hists_m5_b9 = foo.yieldGraphs[15:20]
cols = [ROOT.kRed,ROOT.kGreen,ROOT.kBlue,ROOT.kViolet,ROOT.kCyan]
labs = ['MJ > 600 GeV','MJ > 650 GeV','MJ > 700 GeV','MJ > 750 GeV','MJ > 800 GeV']
leg = ROOT.TLegend(0.75,0.65,0.8,0.9)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.03)
lumiLatex=ROOT.TLatex()

for i in range(len(hists_m5_b1)):
    leg.AddEntry(hists_m5_b1[i],labs[i],'l')

c1=ROOT.TCanvas('c1','c1',800,600)
drawn=False
for i in range(len(hists_m4_b1)):
    if hists_m4_b1[i] == 0:
        continue

    hists_m4_b1[i].SetLineColor(cols[i])
    hists_m4_b1[i].SetLineWidth(2)
    if not drawn:
        drawn=True
        hists_m4_b1[i].Draw('AL')
        hists_m4_b1[i].GetXaxis().SetLimits(700,2200)
        hists_m4_b1[i].GetYaxis().SetRangeUser(0,1600)
        c1.Modified()
        c1.Update()
    else:
        hists_m4_b1[i].Draw('L')

ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 7.1 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 4}{b-tag}')
leg.Draw()

c2=ROOT.TCanvas('c2','c2',800,600)
drawn=False
for i in range(len(hists_m4_b9)):
    if hists_m4_b9[i] == 0:
        continue
    hists_m4_b9[i].SetLineColor(cols[i])
    hists_m4_b9[i].SetLineWidth(2)
    if not drawn:
        drawn=True
        hists_m4_b9[i].Draw('AL')
        hists_m4_b9[i].GetXaxis().SetLimits(700,2200)
        hists_m4_b9[i].GetYaxis().SetRangeUser(0,1600)
        c2.Modified()
        c2.Update()
    else:
        hists_m4_b9[i].Draw('L')

ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 7.1 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 4}{b-inclusive}')
leg.Draw()

c3=ROOT.TCanvas('c3','c3',800,600)
for i in range(len(hists_m5_b1)):
    if hists_m5_b1[i] == 0:
        continue
    hists_m5_b1[i].SetLineColor(cols[i])
    hists_m5_b1[i].SetLineWidth(2)
    if i == 0:
        hists_m5_b1[i].Draw('AL')
        hists_m5_b1[i].GetXaxis().SetLimits(800,2200)
        hists_m5_b1[i].GetYaxis().SetRangeUser(0,1600)
        c3.Modified()
        c3.Update()
    else:
        hists_m5_b1[i].Draw('L')

ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 7.1 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 5}{b-tag}')
leg.Draw()

c4=ROOT.TCanvas('c4','c4',800,600)
for i in range(len(hists_m5_b9)):
    if hists_m5_b9[i] == 0:
        continue
    hists_m5_b9[i].SetLineColor(cols[i])
    hists_m5_b9[i].SetLineWidth(2)
    if i == 0:
        hists_m5_b9[i].Draw('AL')
        hists_m5_b9[i].GetXaxis().SetLimits(700,2200)
        hists_m5_b9[i].GetYaxis().SetRangeUser(0,1600)
        c4.Modified()
        c4.Update()
    else:
        hists_m5_b9[i].Draw('L')
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 7.1 fb^{-1}')
lumiLatex.DrawLatexNDC(0.75,0.5,'#splitline{n_{jet}#geq 5}{b-inclusive}')
leg.Draw()

c5 = ROOT.TCanvas('c5','c5',800,600)
foo.limit_1up.SetFillColor(5)
foo.limit_1up.Draw('aF')
foo.limit_1up.GetXaxis().SetLimits(700,2200)
foo.limit_1up.GetYaxis().SetRangeUser(0,1600)
foo.limit_1down.SetFillColor(10)
foo.limit_1down.Draw('Fsame')
foo.limit_nom.Draw('lsame')
ROOT.ATLASLabel(0.2,0.85,'Internal')
lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 7.1 fb^{-1}')
lumiLatex.DrawLatexNDC(0.65,0.4,'#splitline{#splitline{n_{jet}#geq 5}{b-tag}}{MJ > 800 GeV}')
line = ROOT.TLine()
line.SetLineColor(ROOT.kWhite)
line.SetLineWidth(3)
#line.DrawLine(1000,200,1000,800)

#line.DrawLine(1000,268.882177216,1000,828.200665761)
#line.DrawLine(1000,828.200665761,1006.35299758,850)
#x1  =1006.35299758
#y1 = 850
#x3 = 1200.0
#y3 = 946.644554518
#x2 = 1190
x1=1000.0
y1=238.028
x3=1200.0
y3=975.051254841
x2=1186
y2=y1+(x2-x1)*(y3-y1)/(x3-x1)

line.SetLineColor(ROOT.kWhite)
line.DrawLine(x1,y1,x2,y2)

#3 sigma
#c1.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m4_b1_dy14_3sigma.pdf')
#c2.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m4_b9_dy14_3sigma.pdf')
#c3.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m5_b1_dy14_3sigma.pdf')
#c4.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m5_b9_dy14_3sigma.pdf')
#c5.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m5_b1_dy14_mj800_3sigma_1sigmaband.pdf')

#95% CL
#c1.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m4_b1_dy14_95CL.pdf')
#c2.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m4_b9_dy14_95CL.pdf')
#c3.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m5_b1_dy14_95CL.pdf')
#c4.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m5_b9_dy14_95CL.pdf')
#c5.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_5p8fb/reach_RPV10_m5_b1_dy14_mj800_95CL_1sigmaband.pdf')


#3 sigma
#c1.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m4_b1_dy14_3sigma.pdf')
#c2.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m4_b9_dy14_3sigma.pdf')
#c3.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m5_b1_dy14_3sigma.pdf')
#c4.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m5_b9_dy14_3sigma.pdf')
#c5.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m5_b1_dy14_mj800_3sigma_1sigmaband.pdf')

#95% CL
c1.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m4_b1_dy14_95CL.pdf')
c2.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m4_b9_dy14_95CL.pdf')
c3.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m5_b1_dy14_95CL.pdf')
c4.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m5_b9_dy14_95CL.pdf')
c5.Print('/global/project/projectdirs/atlas/www/multijet/RPV/btamadio/ReachPlots/06_28_7p14fb/reach_RPV10_m5_b1_dy14_mj800_95CL_1sigmaband.pdf')


# c5 = ROOT.TCanvas('c5','c5',800,600)
# foo.limit_2up.Draw('AFL')
# c5.Modified()
# #foo.limit_2up.GetXaxis().SetRangeUser(1200,2000)
# #c5.Update()
# foo.limit_1up.Draw('FL')
# foo.limit_nom.Draw('FL')

# foo.limit_1down.Draw('FL')
# foo.limit_nom.Draw('L')
# foo.limit_2down.Draw('FL')


# foo.limit_2up.GetXaxis().SetLimits(800,1800)
# foo.limit_2up.GetXaxis().SetTitle('m_{#tilde{g}} [GeV]')
# foo.limit_2up.GetYaxis().SetTitle('m_{#tilde{#chi}} [GeV]')

# ROOT.ATLASLabel(0.2,0.85,'Internal')
# lumiLatex.DrawLatexNDC(0.2,0.75,'#int L dt = 7.1 fb^{-1}')
# lumiLatex.DrawLatexNDC(0.7,0.3,'#splitline{#splitline{n_{jet}#geq 5}{b-tag}}{MJ > 800 GeV}')
