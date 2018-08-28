//Testing to use a fit with a macro graph

{
    //Declare Data
    g = new TGraph();
    double x,y;

    //Set Data
    for(int i=0; i<80; i++){
        x = 0.5*i;
        y = 4* i+ 2+ gRandom -> Gaus(0,1);
        g -> SetPoint(i, x, y);
    }

    //Fitting Function
    f1 = new TF1("f1","[0]* x+ 1", 0, 5);
    g -> Fit(f1);

    g -> SetMarkerStyle(20);
    g -> Draw("APC");
}