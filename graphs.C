{
    file = new TFile("data.root");

    TH1* h3 = 0;
    file->GetObject("h3",h3);

    h3->Draw();
}