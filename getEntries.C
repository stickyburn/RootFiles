{

	//Not written as a 'named' macro. 
	//Gets all entries from a combined data file which combines all the histograms.

  //Read the combined data file
	file = new TFile("./RootFiles/Combined.root");
  //Get a list of all histograms
  histList = (TList*)gDirectory->GetListOfKeys();
  //Create an array of dummy histograms
	TH1* hist0Tag[50000];
  TH1* hist2Tag[30000];
  //Loop all data
  for(int i=0; i<histList->GetSize(); i++){
		//Get all names one by one 
    TString mystring(histList->At(i)->GetName());
			if(mystring.Contains("0tag")){
			if(mystring.Contains("PreSel" ) && mystring.Contains("deltaR")){
					hist0Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist0Tag");
					cout<<mystring<<endl;         
					cout<<"Number of enties after PreSel and PostDEtaCut and FJDT2D: "<<hist0Tag[i]->GetEntries()<<endl<<endl;
			}
			}
	// 		if(mystring.Contains("PreSel" ) && mystring.Contains("PostDRCut") && mystring.Contains("FJDT2D")){
	// 							hist1Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist1Tag");
	// 							cout<<mystring<<endl;         
	// 							cout<<"Number of enties after PreSel and PostDRCut and FJDT2D: "<<hist1Tag[i]->GetEntries()<<endl<<endl;
	// 			}
	// 	}			
  //   //Search for histograms with names
  //   if(mystring.Contains("2tag")) {
  //     //Choose the histogram of the deltaR distribution:
  //     if( mystring.Contains("PreSel") &&  mystring.Contains("deltaR")){
	// 		//Saving a copy
	// 		hist2Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist2Tag");
	// 		//Print the total entries on the screen:
	// 		cout<<mystring<<endl;   
	// 		cout<<"Number of enties after PreSel: "<<hist2Tag[i]->GetEntries()<<endl<<endl;
	// 	}
  //     if( mystring.Contains("PreSel") && mystring.Contains("PostDEtaCut") &&  mystring.Contains("deltaR") ){
	// 			hist2Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist2Tag");
	// 			cout<<mystring<<endl;
	// 			cout<<"Number of enties after PreSel and PostDEtaCut: "<<hist2Tag[i]->GetEntries()<<endl<<endl;
	// 		}

  //     if( mystring.Contains("PreSel") && mystring.Contains("PostDR") &&  mystring.Contains("deltaR") ){
	// 			hist2Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist2Tag");
	// 			cout<<mystring<<endl;         
	// 			cout<<"Number of enties after PreSel and PostDR: "<<hist2Tag[i]->GetEntries()<<endl<<endl;
	// }

  //     if( mystring.Contains("PreSel") && mystring.Contains("PostDEtaCut") && mystring.Contains("PostDR") && mystring.Contains("deltaR")){
	// 			hist2Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist2Tag");
	// 			cout<<mystring<<endl;         
	// 			cout<<"Number of enties after PreSel and PostDEtaCut and PostDR: "<<hist2Tag[i]->GetEntries()<<endl<<endl;
  //     }

  //     //SR
  //     if( mystring.Contains("SR") &&  mystring.Contains("deltaR")){
	// 				//Clone the matching histogram into a new histogram (this may not be strickly necessary)
	// 				hist2Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist2Tag");
	// 				//Print the number of entries to the screen:
	// 				cout<<mystring<<endl;         
	// 				cout<<"Number of enties after SR: "<<hist2Tag[i]->GetEntries()<<endl<<endl;
	// 			}
  //     if( mystring.Contains("SR") && mystring.Contains("PostDEtaCut") &&  mystring.Contains("deltaR") ){
	// 			hist2Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist2Tag");
	// 			cout<<mystring<<endl;         
	// 			cout<<"Number of enties after SR and PostDEtaCut: "<<hist2Tag[i]->GetEntries()<<endl<<endl;
	// 			}

  //     if( mystring.Contains("SR") && mystring.Contains("PostDEtaCut") && mystring.Contains("PostDR") && mystring.Contains("deltaR")){
	// 				hist2Tag[i] = (TH1F*)gDirectory->Get(histList->At(i)->GetName())->Clone("hist2Tag");
	// 				cout<<mystring<<endl;         
	// 				cout<<"Number of enties after SR and PostDEtaCut and PostDR: "<<hist2Tag[i]->GetEntries()<<endl<<endl;
  //     	} 
  //   }      
  //   //next if statement goes here
    
  }
}
