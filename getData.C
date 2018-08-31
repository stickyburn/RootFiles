#include <iostream>
#include <fstream>
using namespace std;

void getData(){
    //Opening files sequencially due to gDirect change.
    TFile* fileS = new TFile("./RootFiles/hist-RSG_C10_M1300-0.root");
    TList* signalList = (TList*)gDirectory->GetListOfKeys();
    signalList->Sort();

    TFile* fileB = new TFile("./RootFiles/Combined.root");  
    TList* backList = (TList*)gDirectory->GetListOfKeys();
    backList->Sort();

    //Make histos to copy and use
    TH1* signal2Tag[1000];
    TH1* back2Tag[1000];

    //Get enteries from Data
    double signalEntries[1000];
    double backEntries[1000];

    //Put in array using index
    int signalIndex = 0;
    int backIndex = 0;

    //Looping over all data for SIGNAL
    for(int i=0; i<signalList->GetSize(); i++){
        TString signalString(signalList->At(i)->GetName());
        
        if(signalString.Contains("2tag") ){ 
            if(signalString.Contains("PreSel") && signalString.Contains("PostDRCut")  
            && signalString.Contains("deltaR")){
                signal2Tag[i] = (TH1*)fileS->Get(signalList->At(i)->GetName())->Clone("signal2Tag");
                signalEntries[signalIndex] = signal2Tag[i]->GetEntries();
                signalIndex++;
                }
            if(signalString.Contains("PreSel") && signalString.Contains("PostDEtaCut")  
            && signalString.Contains("deltaR")){
                signal2Tag[i] = (TH1*)fileS->Get(signalList->At(i)->GetName())->Clone("signal2Tag");
                signalEntries[signalIndex] = signal2Tag[i]->GetEntries();
                signalIndex++;
                }
            if(signalString.Contains("SR") && signalString.Contains("PostnBTagTJCut")  
            && signalString.Contains("deltaR")){
                signal2Tag[i] = (TH1*)fileS->Get(signalList->At(i)->GetName())->Clone("signal2Tag");
                signalEntries[signalIndex] = signal2Tag[i]->GetEntries();
                signalIndex++;
                }
            if(signalString.Contains("SR") && signalString.Contains("PostNSubJetsCut")  
            && signalString.Contains("deltaR")){
                signal2Tag[i] = (TH1*)fileS->Get(signalList->At(i)->GetName())->Clone("signal2Tag");
                signalEntries[signalIndex] = signal2Tag[i]->GetEntries();
                signalIndex++;
                }
            if(signalString.Contains("SR") && signalString.Contains("PostEllipseCut")  
            && signalString.Contains("deltaR")){
                signal2Tag[i] = (TH1*)fileS->Get(signalList->At(i)->GetName())->Clone("signal2Tag");
                signalEntries[signalIndex] = signal2Tag[i]->GetEntries();
                signalIndex++;
                }
        }

        }

    //Looping over all data for BACKGROUND
    for(int j=0; j<backList->GetSize(); j++){
        TString backString(backList->At(j)->GetName());

        if(backString.Contains("2tag")){
            if(backString.Contains("PreSel") && backString.Contains("PostDRCut")
            && backString.Contains("deltaR")){
                    back2Tag[j] = (TH1*)fileB->Get(backList->At(j)->GetName())->Clone("back2Tag");
                    backEntries[backIndex] = back2Tag[j]->GetEntries();
                    backIndex++;
                }
            if(backString.Contains("PreSel") && backString.Contains("PostDEtaCut")
            && backString.Contains("deltaR")){
                    back2Tag[j] = (TH1*)fileB->Get(backList->At(j)->GetName())->Clone("back2Tag");
                    backEntries[backIndex] = back2Tag[j]->GetEntries();
                    backIndex++;
                }
            if(backString.Contains("SR") && backString.Contains("PostnBTagTJCut")
            && backString.Contains("deltaR")){
                    back2Tag[j] = (TH1*)fileB->Get(backList->At(j)->GetName())->Clone("back2Tag");
                    backEntries[backIndex] = back2Tag[j]->GetEntries();
                    backIndex++;
                }
            if(backString.Contains("SR") && backString.Contains("PostNSubJetsCut")
            && backString.Contains("deltaR")){
                    back2Tag[j] = (TH1*)fileB->Get(backList->At(j)->GetName())->Clone("back2Tag");
                    backEntries[backIndex] = back2Tag[j]->GetEntries();
                    backIndex++;
                }
            if(backString.Contains("SR") && backString.Contains("PostEllipseCut")
            && backString.Contains("deltaR")){
                    back2Tag[j] = (TH1*)fileB->Get(backList->At(j)->GetName())->Clone("back2Tag");
                    backEntries[backIndex] = back2Tag[j]->GetEntries();
                    backIndex++;
                }
        }
    }

    cout<<endl<<"Comparing with FILE : "<<fileS->GetName()<<endl<<endl;
    cout<<"Ratio of PreSel_PostDRCut_deltaR : "<<signalEntries[0]/backEntries[0]<<endl;
    cout<<"Ratio of PreSel_PostDEtaCut_deltaR : "<<signalEntries[1]/backEntries[1]<<endl;
    cout<<"Ratio of SR_PostnBTagTJCut_deltaR : "<<signalEntries[2]/backEntries[2]<<endl;
    cout<<"Ratio of SR_PostNSubJetsCut_deltaR : "<<signalEntries[3]/backEntries[3]<<endl;
    cout<<"Ratio of SR_PostEllipseCut_deltaR : "<<signalEntries[4]/backEntries[4]<<endl;

    ofstream saveVals ("data.txt");
    if(saveVals){
    for(int s=0; s<backIndex; s++){
        saveVals << (signalEntries[s]/backEntries[s]) << endl;
        }
    }


        // if(signalString.Contains("0ptag") ){
        //     signal0pTag[i] = (TH1*)fileS->Get(signalList->At(i)->GetName())->Clone("signal0pTag");
        //     signalEntries[signalIndex] = signal0pTag[i]->GetEntries();
        //     signalIndex++;
        //     }
        // if(signalString.Contains("1tag")){
        //     signal1Tag[i] = (TH1*)fileS->Get(signalList->At(i)->GetName())->Clone("signal1Tag");
        //     signalEntries[signalIndex] = signal1Tag[i]->GetEntries();
        //     signalIndex++;
        // }
        // if(signalString.Contains("2tag") ){
        //     signal2Tag[i] = (TH1*)fileS->Get(signalList->At(i)->GetName())->Clone("signal2Tag");
        //     signalEntries[signalIndex] = signal2Tag[i]->GetEntries();
        //     signalIndex++;
        //     }

         // if(backString.Contains("0ptag")){
        //         back0pTag[j] = (TH1*)fileB->Get(backList->At(j)->GetName())->Clone("back0pTag");
        //         backEntries[backIndex] = back0pTag[j]->GetEntries();
        //         backIndex++;
        //     }

        // if(backString.Contains("1tag")){
        //     back1Tag[j] = (TH1*)fileB->Get(backList->At(j)->GetName())->Clone("back1Tag");
        //     backEntries[backIndex] = back1Tag[j]->GetEntries();
        //     backIndex++;
        // }
        // if(backString.Contains("2tag")){
        //         back2Tag[j] = (TH1*)fileB->Get(backList->At(j)->GetName())->Clone("back2Tag");
        //         backEntries[backIndex] = back2Tag[j]->GetEntries();
        //         backIndex++;
        //     }
}