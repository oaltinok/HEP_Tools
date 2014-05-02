{
//=========Macro generated from canvas: c1/c1
//=========  (Fri May  2 13:07:08 2014) by ROOT version5.30/00
   TCanvas *c1 = new TCanvas("c1", "c1",429,73,700,500);
   c1->Range(-0.125,-283.2375,1.125,2549.138);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   THStack *hs = new THStack();
   hs->SetName("hs");
   hs->SetTitle("Proton Score");
   
   TH1F *hs_stack_1 = new TH1F("hs_stack_1","Proton Score",20,0,1);
   hs_stack_1->SetMinimum(0);
   hs_stack_1->SetMaximum(2265.9);
   hs_stack_1->SetDirectory(0);
   hs_stack_1->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   hs_stack_1->SetLineColor(ci);
   hs_stack_1->GetXaxis()->SetTitle("Proton Score");
   hs_stack_1->GetXaxis()->SetLabelFont(42);
   hs_stack_1->GetXaxis()->SetLabelSize(0.035);
   hs_stack_1->GetXaxis()->SetTitleSize(0.035);
   hs_stack_1->GetXaxis()->SetTitleFont(42);
   hs_stack_1->GetYaxis()->SetTitle("Candidates / 0.05 ");
   hs_stack_1->GetYaxis()->SetLabelFont(42);
   hs_stack_1->GetYaxis()->SetLabelSize(0.035);
   hs_stack_1->GetYaxis()->SetTitleSize(0.035);
   hs_stack_1->GetYaxis()->SetTitleFont(42);
   hs_stack_1->GetZaxis()->SetLabelFont(42);
   hs_stack_1->GetZaxis()->SetLabelSize(0.035);
   hs_stack_1->GetZaxis()->SetTitleSize(0.035);
   hs_stack_1->GetZaxis()->SetTitleFont(42);
   hs->SetHistogram(hs_stack_1);
   
   
   TH1F *pID_other = new TH1F("pID_other","Other",20,0,1);
   pID_other->SetBinContent(1,46);
   pID_other->SetBinContent(2,388);
   pID_other->SetBinContent(3,281);
   pID_other->SetBinContent(4,288);
   pID_other->SetBinContent(5,297);
   pID_other->SetBinContent(6,270);
   pID_other->SetBinContent(7,144);
   pID_other->SetBinContent(8,95);
   pID_other->SetBinContent(9,55);
   pID_other->SetBinContent(10,51);
   pID_other->SetBinContent(11,48);
   pID_other->SetBinContent(12,57);
   pID_other->SetBinContent(13,51);
   pID_other->SetBinContent(14,73);
   pID_other->SetBinContent(15,93);
   pID_other->SetBinContent(16,70);
   pID_other->SetBinContent(17,40);
   pID_other->SetBinContent(18,34);
   pID_other->SetBinContent(19,33);
   pID_other->SetBinContent(20,10);
   pID_other->SetBinContent(21,2);
   pID_other->SetEntries(2426);

   ci = TColor::GetColor("#ff0000");
   pID_other->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   pID_other->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   pID_other->SetMarkerColor(ci);
   pID_other->SetMarkerStyle(21);
   pID_other->GetXaxis()->SetTitle("Other");
   pID_other->GetXaxis()->SetLabelFont(42);
   pID_other->GetXaxis()->SetLabelSize(0.035);
   pID_other->GetXaxis()->SetTitleSize(0.035);
   pID_other->GetXaxis()->SetTitleFont(42);
   pID_other->GetYaxis()->SetTitle("Candidates / 0.1 ");
   pID_other->GetYaxis()->SetLabelFont(42);
   pID_other->GetYaxis()->SetLabelSize(0.035);
   pID_other->GetYaxis()->SetTitleSize(0.035);
   pID_other->GetYaxis()->SetTitleFont(42);
   pID_other->GetZaxis()->SetLabelFont(42);
   pID_other->GetZaxis()->SetLabelSize(0.035);
   pID_other->GetZaxis()->SetTitleSize(0.035);
   pID_other->GetZaxis()->SetTitleFont(42);
   hs->Add(pID_other,"");
   
   TH1F *pID_piminus = new TH1F("pID_piminus","Pi Minus",20,0,1);
   pID_piminus->SetBinContent(1,40);
   pID_piminus->SetBinContent(2,326);
   pID_piminus->SetBinContent(3,189);
   pID_piminus->SetBinContent(4,144);
   pID_piminus->SetBinContent(5,134);
   pID_piminus->SetBinContent(6,86);
   pID_piminus->SetBinContent(7,45);
   pID_piminus->SetBinContent(8,25);
   pID_piminus->SetBinContent(9,21);
   pID_piminus->SetBinContent(10,25);
   pID_piminus->SetBinContent(11,36);
   pID_piminus->SetBinContent(12,37);
   pID_piminus->SetBinContent(13,22);
   pID_piminus->SetBinContent(14,29);
   pID_piminus->SetBinContent(15,27);
   pID_piminus->SetBinContent(16,9);
   pID_piminus->SetBinContent(17,13);
   pID_piminus->SetBinContent(18,10);
   pID_piminus->SetBinContent(19,11);
   pID_piminus->SetBinContent(20,5);
   pID_piminus->SetEntries(1234);

   ci = TColor::GetColor("#ffff00");
   pID_piminus->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   pID_piminus->SetLineColor(ci);

   ci = TColor::GetColor("#ffff00");
   pID_piminus->SetMarkerColor(ci);
   pID_piminus->SetMarkerStyle(21);
   pID_piminus->GetXaxis()->SetTitle("Pi Minus");
   pID_piminus->GetXaxis()->SetLabelFont(42);
   pID_piminus->GetXaxis()->SetLabelSize(0.035);
   pID_piminus->GetXaxis()->SetTitleSize(0.035);
   pID_piminus->GetXaxis()->SetTitleFont(42);
   pID_piminus->GetYaxis()->SetTitle("Candidates / 0.1 ");
   pID_piminus->GetYaxis()->SetLabelFont(42);
   pID_piminus->GetYaxis()->SetLabelSize(0.035);
   pID_piminus->GetYaxis()->SetTitleSize(0.035);
   pID_piminus->GetYaxis()->SetTitleFont(42);
   pID_piminus->GetZaxis()->SetLabelFont(42);
   pID_piminus->GetZaxis()->SetLabelSize(0.035);
   pID_piminus->GetZaxis()->SetTitleSize(0.035);
   pID_piminus->GetZaxis()->SetTitleFont(42);
   hs->Add(pID_piminus,"");
   
   TH1F *pID_piplus = new TH1F("pID_piplus","Pi Plus",20,0,1);
   pID_piplus->SetBinContent(1,162);
   pID_piplus->SetBinContent(2,916);
   pID_piplus->SetBinContent(3,505);
   pID_piplus->SetBinContent(4,394);
   pID_piplus->SetBinContent(5,277);
   pID_piplus->SetBinContent(6,200);
   pID_piplus->SetBinContent(7,124);
   pID_piplus->SetBinContent(8,87);
   pID_piplus->SetBinContent(9,66);
   pID_piplus->SetBinContent(10,52);
   pID_piplus->SetBinContent(11,121);
   pID_piplus->SetBinContent(12,89);
   pID_piplus->SetBinContent(13,88);
   pID_piplus->SetBinContent(14,80);
   pID_piplus->SetBinContent(15,75);
   pID_piplus->SetBinContent(16,50);
   pID_piplus->SetBinContent(17,41);
   pID_piplus->SetBinContent(18,30);
   pID_piplus->SetBinContent(19,17);
   pID_piplus->SetBinContent(20,7);
   pID_piplus->SetEntries(3381);

   ci = TColor::GetColor("#0000ff");
   pID_piplus->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   pID_piplus->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   pID_piplus->SetMarkerColor(ci);
   pID_piplus->SetMarkerStyle(21);
   pID_piplus->GetXaxis()->SetTitle("Pi Plus");
   pID_piplus->GetXaxis()->SetLabelFont(42);
   pID_piplus->GetXaxis()->SetLabelSize(0.035);
   pID_piplus->GetXaxis()->SetTitleSize(0.035);
   pID_piplus->GetXaxis()->SetTitleFont(42);
   pID_piplus->GetYaxis()->SetTitle("Candidates / 0.1 ");
   pID_piplus->GetYaxis()->SetLabelFont(42);
   pID_piplus->GetYaxis()->SetLabelSize(0.035);
   pID_piplus->GetYaxis()->SetTitleSize(0.035);
   pID_piplus->GetYaxis()->SetTitleFont(42);
   pID_piplus->GetZaxis()->SetLabelFont(42);
   pID_piplus->GetZaxis()->SetLabelSize(0.035);
   pID_piplus->GetZaxis()->SetTitleSize(0.035);
   pID_piplus->GetZaxis()->SetTitleFont(42);
   hs->Add(pID_piplus,"");
   
   TH1F *pID_proton = new TH1F("pID_proton","Proton",20,0,1);
   pID_proton->SetBinContent(1,78);
   pID_proton->SetBinContent(2,528);
   pID_proton->SetBinContent(3,331);
   pID_proton->SetBinContent(4,298);
   pID_proton->SetBinContent(5,295);
   pID_proton->SetBinContent(6,276);
   pID_proton->SetBinContent(7,239);
   pID_proton->SetBinContent(8,182);
   pID_proton->SetBinContent(9,171);
   pID_proton->SetBinContent(10,172);
   pID_proton->SetBinContent(11,184);
   pID_proton->SetBinContent(12,199);
   pID_proton->SetBinContent(13,201);
   pID_proton->SetBinContent(14,238);
   pID_proton->SetBinContent(15,298);
   pID_proton->SetBinContent(16,332);
   pID_proton->SetBinContent(17,360);
   pID_proton->SetBinContent(18,403);
   pID_proton->SetBinContent(19,375);
   pID_proton->SetBinContent(20,112);
   pID_proton->SetEntries(5272);

   ci = TColor::GetColor("#00ff00");
   pID_proton->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   pID_proton->SetLineColor(ci);

   ci = TColor::GetColor("#00ff00");
   pID_proton->SetMarkerColor(ci);
   pID_proton->SetMarkerStyle(21);
   pID_proton->GetXaxis()->SetTitle("Proton");
   pID_proton->GetXaxis()->SetLabelFont(42);
   pID_proton->GetXaxis()->SetLabelSize(0.035);
   pID_proton->GetXaxis()->SetTitleSize(0.035);
   pID_proton->GetXaxis()->SetTitleFont(42);
   pID_proton->GetYaxis()->SetTitle("Candidates / 0.1 ");
   pID_proton->GetYaxis()->SetLabelFont(42);
   pID_proton->GetYaxis()->SetLabelSize(0.035);
   pID_proton->GetYaxis()->SetTitleSize(0.035);
   pID_proton->GetYaxis()->SetTitleFont(42);
   pID_proton->GetZaxis()->SetLabelFont(42);
   pID_proton->GetZaxis()->SetLabelSize(0.035);
   pID_proton->GetZaxis()->SetTitleSize(0.035);
   pID_proton->GetZaxis()->SetTitleFont(42);
   hs->Add(pID_proton,"");
   hs->Draw("");
   
   TPaveText *pt = new TPaveText(0.3822414,0.94,0.6177586,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("Proton Score");
   pt->Draw();
   
   TLegend *leg = new TLegend(0.5,0.67,0.88,0.88,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("pID_other","Other","lpf");

   ci = TColor::GetColor("#ff0000");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#000099");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("pID_piminus","Pi Minus","lpf");

   ci = TColor::GetColor("#ffff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#000099");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ffff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("pID_piplus","Pi Plus","lpf");

   ci = TColor::GetColor("#0000ff");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#000099");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("pID_proton","Proton","lpf");

   ci = TColor::GetColor("#00ff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#000099");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
