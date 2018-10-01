{
    mGraph = new TMultiGraph();
    mGraph->SetTitle("Ratios of Pre & Post Cuts; Number of Cuts; sqrt(Signal/Background)");


    g1 = new TGraph("./Ratios/ratiosM1000.txt");
    g1->SetLineWidth(2);
    mGraph->Add(g1);

    g2 = new TGraph("./Ratios/ratiosM1100.txt");
    g2->SetLineColor(kBlue);
    g2->SetLineWidth(2);
    mGraph->Add(g2);

    g3 = new TGraph("./Ratios/ratiosM1300.txt");
    g3->SetLineColor(kBlack);
    g3->SetLineWidth(2);
    mGraph->Add(g3);

    g4 = new TGraph("./Ratios/ratiosM1500.txt");
    g4->SetLineColor(kRed);
    g4->SetLineWidth(2);
    mGraph->Add(g4);

    g5 = new TGraph("./Ratios/ratiosM1600.txt");
    g5->SetLineColor(kYellow+2);
    g5->SetLineWidth(2);
    mGraph->Add(g5);

    g6 = new TGraph("./Ratios/ratiosM1800.txt");
    g6->SetLineColor(kPink+10);
    g6->SetLineWidth(2);
    mGraph->Add(g6);

    g7 = new TGraph("./Ratios/ratiosM2000.txt");
    g7->SetLineColor(kRed);
    g7->SetLineWidth(2);
    mGraph->Add(g7);

    g8 = new TGraph("./Ratios/ratiosM2250.txt");
    g8->SetLineColor(kBlue+2);
    g8->SetLineWidth(2);
    mGraph->Add(g8);

    g9 = new TGraph("./Ratios/ratiosM2500.txt");
    g9->SetLineColor(kPink);
    g9->SetLineWidth(2);
    mGraph->Add(g9);

    g10 = new TGraph("./Ratios/ratiosM2750.txt");
    g10->SetLineColor(kBlue+3);
    g10->SetLineWidth(2);
    mGraph->Add(g10);

    g11 = new TGraph("./Ratios/ratiosM3000.txt");
    g11->SetLineColor(kRed+3);
    g11->SetLineWidth(2);
    mGraph->Add(g11);

    mGraph->Draw("AL");

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("Mass Legend","C"); 
    legend->AddEntry(g1,"Mass 1000 GeV");
    legend->AddEntry(g2,"Mass 1100 GeV");
    legend->AddEntry(g3,"Mass 1300 GeV");
    legend->AddEntry(g4,"Mass 1500 GeV");
    legend->AddEntry(g5,"Mass 1600 GeV");
    legend->AddEntry(g6,"Mass 1800 GeV");
    legend->AddEntry(g7,"Mass 2000 GeV");
    legend->AddEntry(g8,"Mass 2250 GeV");
    legend->AddEntry(g9,"Mass 2500 GeV");
    legend->AddEntry(g10,"Mass 2750 GeV");
    legend->AddEntry(g11,"Mass 3000 GeV");

  auto legend1 = new TLegend(0.1,0.7,0.48,0.9);
   legend1->SetHeader("Cuts Legend","C"); 
   legend1->AddEntry("1","1 - PreSel_PostDEtaCut_DeltaR");
   legend1->AddEntry("2","2 - PreSel_PostDRCut_DeltaR");
   legend1->AddEntry("3","3 - SR_PostEllipseCut_DeltaR");
   legend1->AddEntry("4","4 - SR_PostNSubJetsCut_DeltaR");
   legend1->AddEntry("5","5 - SR_PostnBTagTJCut_DeltaR");

   legend->Draw();
   legend1->Draw("SAME");
}