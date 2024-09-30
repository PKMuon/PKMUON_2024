void count(TCut cut, TChain* t_air, TTree* t_const_0p005, TTree* t_const_0p05, TTree* t_const_0p5, TTree* t_const_0p2, TTree* t_const_0p1, TTree* t_const_1, TTree* t_const_10, TTree* t_const_100);
double punzi(double eff, double Nbkg, double a);

void CutScan(){

TFile* f_const_0p005 = new TFile("../root/DMmuon_const_0p005.root","");
TFile* f_const_0p05 = new TFile("../root/DMmuon_const_0p05.root","");
TFile* f_const_0p5 = new TFile("../root/DMmuon_const_0p5.root","");
//TFile* f_const_0p4 = new TFile("../root/DMmuon_const_0p4.root","");
//TFile* f_const_0p3 = new TFile("../root/DMmuon_const_0p3.root","");
TFile* f_const_0p2 = new TFile("../root/DMmuon_const_0p2.root","");
TFile* f_const_0p1 = new TFile("../root/DMmuon_const_0p1.root","");
TFile* f_const_1 = new TFile("../root/DMmuon_const_1.root","");
TFile* f_const_10 = new TFile("../root/DMmuon_const_10.root","");
TFile* f_const_100 = new TFile("../root/DMmuon_const_100.root","");
//TFile* f_const_1000 = new TFile("../root/DMmuon_const_1000.root","");

TTree* t_const_0p005 = (TTree*)f_const_0p005->Get("T1");
TTree* t_const_0p05 = (TTree*)f_const_0p05->Get("T1");
TTree* t_const_0p5 = (TTree*)f_const_0p5->Get("T1");
//TTree* t_const_0p4 = (TTree*)f_const_0p4->Get("T1");
//TTree* t_const_0p3 = (TTree*)f_const_0p3->Get("T1");
TTree* t_const_0p2 = (TTree*)f_const_0p2->Get("T1");
TTree* t_const_0p1 = (TTree*)f_const_0p1->Get("T1");
TTree* t_const_1 = (TTree*)f_const_1->Get("T1");
TTree* t_const_10 = (TTree*)f_const_10->Get("T1");
TTree* t_const_100 = (TTree*)f_const_100->Get("T1");
//TTree* t_const_1000 = (TTree*)f_const_1000->Get("T1");

TChain* t_air_new = new TChain("T1");
//t_air_new->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job/analysis/root/muCRY_air_1.root");
t_air_new->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job/analysis/root/muCRY_air_*.root");

TCut cut;

cut = "1";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);

cout<<"DCA cut"<<endl;
cut = "dca_smear<500";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
cut = "dca_smear<400";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
cut = "dca_smear<300";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
cut = "dca_smear<200";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
cut = "dca_smear<100";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);


cout<<"PoCA cut"<<endl;
cut = "fabs(poca_x_smear)<1000 && fabs(poca_y_smear)<1000 && fabs(poca_z_smear)<1000";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
cut = "fabs(poca_x_smear)<800 && fabs(poca_y_smear)<800 && fabs(poca_z_smear)<800";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
cut = "fabs(poca_x_smear)<600 && fabs(poca_y_smear)<600 && fabs(poca_z_smear)<600";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
cut = "fabs(poca_x_smear)<500 && fabs(poca_y_smear)<500 && fabs(poca_z_smear)<500";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
cut = "fabs(poca_x_smear)<400 && fabs(poca_y_smear)<400 && fabs(poca_z_smear)<400";
count(cut, t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);

} 

void count(TCut cut, TChain* t_air, TTree* t_const_0p005, TTree* t_const_0p05, TTree* t_const_0p5, TTree* t_const_0p2, TTree* t_const_0p1, TTree* t_const_1, TTree* t_const_10, TTree* t_const_100){

double Nair = t_air->GetEntries(cut);
double Nsig_0p005 = t_const_0p005->GetEntries(cut);
double Nsig_0p05 = t_const_0p05->GetEntries(cut);
double Nsig_0p5 = t_const_0p5->GetEntries(cut);
double Nsig_0p2 = t_const_0p2->GetEntries(cut);
double Nsig_0p1 = t_const_0p1->GetEntries(cut);
double Nsig_1 = t_const_1->GetEntries(cut);
double Nsig_10 = t_const_10->GetEntries(cut);
double Nsig_100 = t_const_100->GetEntries(cut);

double Nsig_tot = 1000000.;
double eff_0p005 = Nsig_0p005 / Nsig_tot;
double eff_0p05 = Nsig_0p05 / Nsig_tot;
double eff_0p5 = Nsig_0p5 / Nsig_tot;
double eff_0p2 = Nsig_0p2 / Nsig_tot;
double eff_0p1 = Nsig_0p1 / Nsig_tot;
double eff_1 = Nsig_1 / Nsig_tot;
double eff_10 = Nsig_10 / Nsig_tot;
double eff_100 = Nsig_100 / Nsig_tot;

cout<<cut.GetTitle()<<endl;
cout<<"Nbkg(1y) = "<<Nair<<endl;
cout<<"DM 0.005 GeV\t Eff = "<<eff_0p005<<" ,\t FOM = "<<punzi(eff_0p005,Nair,3)<<endl;
cout<<"DM 0.05 GeV\t Eff = "<<eff_0p05<<" ,\t FOM = "<<punzi(eff_0p05,Nair,3)<<endl;
cout<<"DM 0.5 GeV\t Eff = "<<eff_0p5<<" ,\t FOM = "<<punzi(eff_0p5,Nair,3)<<endl;
cout<<"DM 0.2 GeV\t Eff = "<<eff_0p2<<" ,\t FOM = "<<punzi(eff_0p2,Nair,3)<<endl;
cout<<"DM 0.1 GeV\t Eff = "<<eff_0p1<<" ,\t FOM = "<<punzi(eff_0p1,Nair,3)<<endl;
cout<<"DM 1 GeV\t Eff = "<<eff_1<<" ,\t FOM = "<<punzi(eff_1,Nair,3)<<endl;
cout<<"DM 10 GeV\t Eff = "<<eff_10<<" ,\t FOM = "<<punzi(eff_10,Nair,3)<<endl;
cout<<"DM 100 GeV\t Eff = "<<eff_100<<" ,\t FOM = "<<punzi(eff_100,Nair,3)<<endl;
cout<<endl;

}

double punzi(double eff, double Nbkg, double a)
{
	double punzi = eff / (a/2+sqrt(Nbkg));
	return punzi;
}
