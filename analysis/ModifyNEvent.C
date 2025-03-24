const char *rootfiles[] = { };

void ModifyNEvent()
{
  for(const char *rootfile : rootfiles) {
    cout << rootfile << endl;
    TFile *file = new TFile(rootfile, "UPDATE");
    TTree *meta = (TTree *)file->Get("meta");
    Double_t NEvent = 0;
    meta->SetBranchAddress("Meta.NEvent", &NEvent);
    TTree *meta_new = meta->CloneTree(0);
    meta_new->SetBranchAddress("Meta.NEvent", &NEvent);
    for(Long64_t i = 0; i < meta->GetEntries(); ++i) {
      meta->GetEntry(i);
      NEvent = 1e7;
      meta_new->Fill();
    }
    meta_new->Write("meta", meta_new->kOverwrite);
    file->Close();
  }
}
