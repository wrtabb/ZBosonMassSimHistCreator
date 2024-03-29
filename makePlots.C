TString canvas_name = "canvas";
int canvas_number = 0;
void Plot2D(TString tag,TH2D*hist);
void Plot1D(TString tag,TH1D*hReco,TH1D*hHard,TH1D*hDressed);
void PlotProjections(TString tag,TH2D*hMatrix,TH1D*hReco,TH1D*hTrue);

void makePlots(TString hist_location)
{
    gROOT->SetBatch(true);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    TString input_file = "output_data/";
    input_file += hist_location;
    input_file += "/unfolding_histograms.root";

    // File to open
    TFile*open_file = new TFile(input_file);
    
    // Get histograms from file
    TH2D*hMatrixHard    = (TH2D*)open_file->Get("hMatrixHard");
    TH2D*hMatrixDressed = (TH2D*)open_file->Get("hMatrixDressed");
    TH1D*hReco          = (TH1D*)open_file->Get("hInvMassReco");
    TH1D*hHard          = (TH1D*)open_file->Get("hInvMassHard");
    TH1D*hDressed       = (TH1D*)open_file->Get("hInvMassDressed");

    hReco->SetMarkerStyle(20);
    hReco->SetMarkerColor(kBlack);
    hHard->SetLineColor(kRed);
    hDressed->SetLineColor(kBlue);

    Plot2D(hist_location,hMatrixHard);
    Plot2D(hist_location,hMatrixDressed);
    Plot1D(hist_location,hReco,hHard,hDressed);
    PlotProjections(hist_location,hMatrixHard,hReco,hHard);
    PlotProjections(hist_location,hMatrixDressed,hReco,hDressed);
}

void Plot2D(TString tag,TH2D*hist)
{
    TString canName = canvas_name;
    canName += canvas_number;
    TCanvas*canvas = new TCanvas(canName,"",0,0,1000,1000);
    canvas->SetGrid();
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.15);

    hist->GetYaxis()->SetTitle("m^{reco} [GeV]");
    hist->GetXaxis()->SetTitle("m^{true} [GeV]");
    hist->Draw("colz");

    canvas_number++;

    TString savename = "plots/";
    savename += tag;
    savename += hist->GetName();
    savename += ".png";
    canvas->SaveAs(savename);
}

void Plot1D(TString tag,TH1D*hReco,TH1D*hHard,TH1D*hDressed)
{
    TString canName = canvas_name;
    canName += canvas_number;
    TCanvas*canvas = new TCanvas(canName,"",0,0,1000,1000);
    canvas->SetGrid();
    canvas->SetLogy();

    TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
    legend->SetTextSize(0.02);
    legend->AddEntry(hReco,"reco");
    legend->AddEntry(hHard,"hard");
    legend->AddEntry(hDressed,"dressed");

    hReco->SetMaximum(1e7);
    hReco->GetXaxis()->SetTitle("m [GeV]");
    hReco->Draw("pe");
    hHard->Draw("hist,same");
    hDressed->Draw("hist,same");
    legend->Draw("same");
    canvas_number++;

    TString savename = "plots/";
    savename += tag;
    savename += "1DPlots";
    savename += ".png";
    canvas->SaveAs(savename);
}

void PlotProjections(TString tag,TH2D*hMatrix,TH1D*hReco,TH1D*hTrue)
{
    // This is a simple check to make sure that the migration matrix
    // Matches the 1D distributions
    TString canName = canvas_name;
    canName += canvas_number;
    TCanvas*canvas = new TCanvas(canName,"",0,0,1000,1000);
    canvas->SetGrid();
    canvas->SetLogy();

    TH1D*projX = (TH1D*)hMatrix->ProjectionX();
    projX->SetMarkerStyle(20);
    projX->SetMarkerColor(kRed);
    TH1D*projY = (TH1D*)hMatrix->ProjectionY();
    projY->SetMarkerStyle(20);
    projY->SetMarkerColor(kBlue);
    
    hReco->SetMaximum(1e7);

    TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
    legend->SetTextSize(0.02);
    legend->AddEntry(hReco,"reco");
    legend->AddEntry(hTrue,"true");
    legend->AddEntry(projX,"matrix x-projection");
    legend->AddEntry(projY,"matrix y-projection");

    hReco->GetXaxis()->SetTitle("m [GeV]");
    hReco->Draw("hist");
    hTrue->Draw("hist,same");
    projX->Draw("pe,same");
    projY->Draw("pe,same");
    legend->Draw("same");
    canvas_number++;

    TString savename = "plots/";
    savename += tag;
    savename += "Projections";
    savename += ".png";
    canvas->SaveAs(savename);
}
