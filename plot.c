void plot() {

//x axis:
float x[5] = {1.,2.,3.,4.,5.};
//y axis:
float y[5] = {20., 40., 60., 80., 100.};


Plot_Example = new TGraph(5, x, y);
Plot_Example ->Draw("AC*");

}
