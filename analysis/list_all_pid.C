void list_all_pid(const char *infile = "../build/root_file/poca_CryMu.root") {

  TFile *file = TFile::Open(infile);
  TTree *tree = (TTree *)file->Get("tree");

  std::map<Int_t, Float_t> pid_count;

  vector<Int_t> *Pid = nullptr;
  tree->SetBranchAddress("Pid", &Pid);

  Long64_t nentries = tree->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    tree->GetEntry(ientry);

    for (size_t i = 0; i < Pid->size(); ++i) {
      pid_count[(*Pid)[i]]++;
    }
  }

  std::cout << "PID Types and Frequencies:" << std::endl;
  for (const auto &entry : pid_count) {
    std::cout << "PID: " << entry.first << ", Count: " << entry.second / 4 << std::endl;
  }

  file->Close();
}
