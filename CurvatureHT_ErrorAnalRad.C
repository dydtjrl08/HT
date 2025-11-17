#include "../Curvature_ht.hpp"
#include "../ParticleMomentum.hpp"


void CurvatureHT_ErrorAnalRad(TString particle,double driftlength,double radian){

	TDatime now;
	
	int dateInt = now.GetDate();

	Double_t min_degree = 0.;
	Double_t max_degree = 180.;
	Double_t delta_degree = 1.;

	gStyle -> SetOptStat(0);
	
	
	std::vector<double> MeV = {2,4,6,8,10,12,14,16,18,20};	
//	std::vector<double> MeV = {20};	
	
	const int n_points = MeV.size();

	std::vector<double> RelativeError(n_points);	
	std::vector<double> RelativeErrorConfidence(n_points);	

	auto cvs1 = new TCanvas(Form("%d_%s_dl = %fmm_AzimuthalAngle=%.1fRadian",dateInt,particle.Data(),driftlength,radian),"cvs1",1200,700);

	for(auto energy : MeV)		
	{
	
//	TVirtualPad* pad;
	
	TString nameHist = Form("hist_%s_%dMeV_%fmm(dl)",particle.Data(),int(energy),driftlength);

	auto mc_result = new TFile(Form("/data/yongseok/%s/%dMeV/SimPadData%dmm_modified.root",particle.Data(),int(energy),int(driftlength)),"read");

	Double_t Radius_dot_Mag = MomentumCalculator::GetRadius(particle.Data(),energy,1.); 

	if(!mc_result->IsOpen() || mc_result -> IsZombie()){
		cout << "file error" << endl;
		return;
	}
	else 
		cout << Form("%s %d MeV (%dmm of driftlength) simulation file is opened. ",particle.Data(),int(energy),int(driftlength)) << endl;


	auto tree = dynamic_cast<TTree*>(mc_result -> Get("data"));

	Double_t HitPad[50][50];
	Double_t TimePad[50][50];
	Double_t PositionPad[50][50][2];
	Double_t ADCSumRow[50];


	tree -> SetBranchAddress("HitPad",&HitPad);
	tree -> SetBranchAddress("TimePad",&TimePad);
	tree -> SetBranchAddress("PositionPad",&PositionPad);
	tree -> SetBranchAddress("SumRowPad",&ADCSumRow);
	
	UInt_t events = tree -> GetEntries();



	TH1D* relative_error = new TH1D(nameHist,nameHist+"Relative error stats of pT/q(MeV / Columb);",100,-0.5,0.5);

	
	TH2D* HoughSpace = new TH2D("HoughSpace","(X): Incident Degree vs. Candidate (B*R)^(-1);Incident degree [rad];BR_reco^-1 [T * mm]^-1",(max_degree - min_degree)/delta_degree,min_degree * (TMath::Pi() / 180.0),max_degree * (TMath::Pi() / 180.0),(max_degree - min_degree)/delta_degree,-0.005,0.005);
	

	//HoughSpace -> Draw("COLZ");


	Double_t fx;
	Double_t fy;
	Double_t ADC;
	Double_t ht_curvature;

	Int_t binx, biny, binz;
	Int_t globalBin;

	Double_t cand_rad;
	Double_t cand_curvature;

	Double_t cand_Radius;

	double epsilon = 1e-3;

	Double_t rad_val_center;

	for(UInt_t event = 0; event < events; event++){

		tree -> GetEntry(event);

		for(UInt_t i = 0; i < 50; i++){
			for(UInt_t j = 0; j < 50; j++){
				

				fx = PositionPad[i][j][0] - 30.;
				if(TMath::Abs(fx) <= 0.1 && j == 0) continue;
				fy = PositionPad[i][j][1];
				ADC = HitPad[i][j];
				if(ADC <= 100.) continue;		
				
			


				for(Double_t deg = min_degree+1; deg <= max_degree; deg+= delta_degree){
					
					try{
						ht_curvature = MomentumCalculator::Curvature_ht(fx,fy,deg);		
						
						//if(epsilon > TMath::Abs(ht_curvature)) continue;

						rad_val_center = (deg+0.5) * (TMath::Pi() / 180.0) *delta_degree;

					

						//std::cout << "good curvature : " << ht_curvature << std::endl;	
						
						HoughSpace -> Fill(rad_val_center, ht_curvature,ADC);

						//HoughSpace -> Fill(deg * (TMath::Pi() / 180.0), ht_curvature,ADC);

						//cout << "ht_curvature : " << ht_curvature << endl;
					} catch(const std::runtime_error &e){

						//std::cerr << "fail to calculate curvatrue (skipping event) : " << e.what() << std::endl;
						//std::cout << "skipping calculating event No : " << event << std::endl;
						
						continue;

					}
					
				}
			}
		}

		globalBin = HoughSpace -> GetMaximumBin();

		HoughSpace -> GetBinXYZ(globalBin,binx,biny,binz);

		cand_rad = HoughSpace -> GetXaxis() -> GetBinCenter(binx);
		cand_curvature = HoughSpace -> GetYaxis() -> GetBinCenter(biny);

		cand_Radius = 1./cand_curvature;	
		
//		cout << "cand Radius* B : " << cand_Radius << "(p/q) " << endl;		
//		cout << "True Radius* b : " << Radius_dot_Mag << "(p/q)" << endl;
		TAxis* xaxis = HoughSpace -> GetXaxis();

		Int_t bin_x_start = xaxis -> FindBin(cand_rad-2*delta_degree*(TMath::Pi() / 180.0));
		Int_t bin_x_end = xaxis -> FindBin(cand_rad+2*delta_degree*(TMath::Pi() / 180.0));
		
		TH1D* h_py_slice = HoughSpace -> ProjectionY("h_py_slice", bin_x_start, bin_x_end);
		

//		cout << "nb : " << h_py_slice -> GetNbinsX() << endl;

		//Int_t maxBin = h_py_slice -> GetMaximumBin();
		//Double_t xValueAtMax = h_py_slice -> GetXaxis() -> GetBinCenter(maxBin);

//		h_py_slice -> Smooth(1);
		
		h_py_slice -> Rebin(2);



		TSpectrum spec(10);
		int npeaks = spec.Search(h_py_slice,5,"nobackground",0.05);


		Double_t* xpk = spec.GetPositionX();
		
		for(Int_t k = 0; k < npeaks; ++k){
			
//			cout << "(Candidate Radius * B)^-1 (just peak) = " << xpk[k] << "(p/q)" << endl;


			int ib = h_py_slice -> FindBin(xpk[k]);
		}

//		HoughSpace -> GetYaxis() -> SetRangeUser(0.01,0.02);
//		HoughSpace -> GetXaxis() -> SetRangeUser(0,3.14);

//		HoughSpace -> Draw("colz");

//		h_py_slice -> Draw();

		
		auto nb = h_py_slice -> GetNbinsX();
		auto bw = h_py_slice -> GetBinWidth(1);

//		cout << " nb : " << nb << ", bw : " << bw << endl;


		int best = -1, best_il = 0, best_ir = 0;
		double bestScore = -1, best_x = 0;


		for(int r = 0; r < npeaks; ++r){
			int ib = h_py_slice -> FindBin(xpk[r]);
			double peak = h_py_slice -> GetBinContent(ib);

			if(peak <= 0) continue;

			double level = 0.3 * peak;
			int il = ib, ir = ib;
			while(il > 1 && h_py_slice -> GetBinContent(il) > level) --il;
			while(ir < nb && h_py_slice -> GetBinContent(ir) > level) ++ir;

			int widthBins = ir - il;
			double area = h_py_slice -> Integral(il,ir);

			double cand_x = xpk[r];
			if(widthBins < 5) continue;

			if(area > bestScore) {
				
				bestScore = area;
				best = r; best_il = il; best_ir = ir; best_x = cand_x;
			}
		}
		

		double xL = h_py_slice -> GetBinLowEdge(best_il);
		double xR = h_py_slice -> GetBinLowEdge(best_ir) + h_py_slice -> GetBinWidth(best_ir);
		double sigmaEst = (xR - xL) / 2.355;

		TF1 g("g","gaus",xL,xR);
			
		//TF1 g("g","crystalball",0.002,0.008);


		g.SetParameters(h_py_slice -> GetBinContent(h_py_slice->FindBin(best_x)), best_x, sigmaEst);
		h_py_slice -> Fit(&g,"RQ");
		
//		h_py_slice -> Draw();

//		cout << "Candidate Radius * B = " << 1./g.GetParameter(1) << "(p/q)" << endl;
//		cout << "events : " << event << " " << "Candidate (Radius * B)^(-1) = " << g.GetParameter(1) << "(p/q)^(-1)" << endl;


	//	cout << "Candidate Radius * B = " << 1./xValueAtMax << "(p/q)" <<endl;

		

	//	relative_error -> Fill(1./g.GetParameter(1) - Radius_dot_Mag);
	//	relative_error -> Fill((g.GetParameter(1) - 1./Radius_dot_Mag) / 1./Radius_dot_Mag);
		relative_error -> Fill(((1./g.GetParameter(1)) - Radius_dot_Mag) / Radius_dot_Mag);
		
		HoughSpace -> Reset();
	}
	
	cout << "True (B*R)^(-1)  = " << 1./Radius_dot_Mag << "(p/q)^(-1) " << endl;

	relative_error -> Draw();


	relative_error -> Fit("gaus");

	TF1* FitResult = relative_error->GetFunction("gaus");

//	FitResult -> Draw("SAME");

	cout << "Bias of pT/Q = " << FitResult -> GetParameter(1) << ", 1Sigma = " << FitResult -> GetParameter(2) << endl;

	RelativeError[(int(energy)-2)/2] = FitResult -> GetParameter(1);
	RelativeErrorConfidence[(int(energy)-2)/2] = FitResult -> GetParameter(2);

	relative_error -> Clear();


	mc_result -> Close();


	}


	TGraphErrors *graph = new TGraphErrors(n_points, MeV.data(), RelativeError.data(), nullptr, RelativeErrorConfidence.data());
	graph->SetTitle(Form("%s : Bias / PT_true vs. Energy (dl = %dmm)",particle.Data(),int(driftlength)));
        graph->GetXaxis()->SetTitle("Energy (MeV)");
        graph->GetYaxis()->SetTitle("Relative Error (bias)");
        graph->SetMarkerStyle(kFullCircle);
        graph->SetMarkerSize(1.0);
        graph->SetMarkerColor(kBlue);
        graph->SetLineColor(kBlue);	

	TAxis* yaxis = graph->GetYaxis();
	yaxis -> SetRangeUser(-0.20,0.20);

	graph -> Draw("AP");
	
	gPad -> SetGrid(1);

	cvs1->SaveAs(Form("figures_HT_4/%s.png",cvs1->GetName()));


}
