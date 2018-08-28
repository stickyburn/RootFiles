{
    c1 = new TCanvas("c1","Different Surfaces",1600,600);
    c1 -> Divide(2,1);

    c1 -> cd(1);
    f1 = new TF2("f1","10*cos(x)*sin(y)",-5,5,-5,5);
    f1 -> Draw(); 

    c1 -> cd(2);
    f1 -> Draw("SURF2");
}