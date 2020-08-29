from ROOT import TFile, TH1D, TMultiGraph, TCanvas, TLegend

filename = [];
filename.append("e1MeV_pred_bterm_smear0mm.root")
filename.append("e1MeV_pred_bterm_smear100mm.root")
filename.append("e1MeV_pred_bterm_smear230mm.root")
filename.append("e1MeV_pred_bterm_smear500mm.root")

rsd_graph = [];

for i in filename:
    ff = TFile(i, "read")
    graph = ff.Get("pred")
    rsd_graph.append(graph)

color = [36, 38, 41, 43]
label= ["no smear", "smear 100mm", "smear 230mm", "smear 500mm" ]

mg = TMultiGraph()
led = TLegend()
for i in range(4):
    rsd_graph[i].SetLineColor(color[i])
    rsd_graph[i].SetFillStyle(3003)
    rsd_graph[i].SetLineWidth(2)
    rsd_graph[i].SetFillColor(30)
    mg.Add(rsd_graph[i])
    led.AddEntry(rsd_graph[i], label[i], "L")

cc = TCanvas()
mg.Draw("A3 L")
led.Draw("SAME")
cc.SaveAs("test.root")
