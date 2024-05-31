#include <vector>  
#include "TEfficiency.h"
#include "TGraph.h"
#include "TList.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TPaveText.h"
#include "TDirectory.h"
#include "TClass.h"
#include "TRegexp.h"
#include "TMath.h"
#include "TLatex.h"
#include <regex>
#include "TH1D.h"
#include "THStack.h"
#include <TH1F.h>
#include <TLegend.h>
#include "TPaveStats.h"
#include <algorithm> 
#include "duneanaobj/StandardRecord/StandardRecord.h" 

// This MINERvA matching code is adapting and refined using Richie's logic.
double calculateDotProductAndExtrapolation(double dirX, double dirY, double dirZ,
                                           double start_z, double end_z, double start_x, double end_x,
                                           double start_y, double end_y,
                                           TVector3 tpcTrk_end, double &deltaExtrapY, double &deltaExtrapX) {
    double dotProduct = -999;
    double extrapdZ = start_z - tpcTrk_end.Z();  
    deltaExtrapY = (dirY / dirZ) * extrapdZ + tpcTrk_end.Y() - start_y;  // Calculate extrapolated Y deviation
    deltaExtrapX = (dirX / dirZ) * extrapdZ + tpcTrk_end.X() - start_x;  // Calculate extrapolated X deviation
    double diffExtrap = TMath::Sqrt(TMath::Power(deltaExtrapY - start_y, 2));  // Calculate the extrapolation difference
    dotProduct = dirX * dirX + dirY * dirY + dirZ * dirZ;  // Calculate the dot product
    return dotProduct;
}

bool isMINERvAMatch(TVector3 reco_vertex, caf::StandardRecord* sr) {
    double maxDotProduct = -999;  
    double minDeltaExtrapY = 9999999;  
    double minDeltaExtrapX = 9999999;  
    int maxPartMinerva = -999;  
    int maxTypeMinerva = -999;  
    int maxIxnMinerva = -999;  

    // Iterate over all MINERvA interactions
    for (int i = 0; i < sr->nd.minerva.ixn.size(); i++) {
        // Iterate over all tracks within a MINERvA interaction
        for (int j = 0; j < sr->nd.minerva.ixn[i].ntracks; j++) {
            // Extract start and end positions for the current MINERvA track
            double start_z = sr->nd.minerva.ixn[i].tracks[j].start.z;
            double end_z = sr->nd.minerva.ixn[i].tracks[j].end.z;
            double start_x = sr->nd.minerva.ixn[i].tracks[j].start.x;
            double end_x = sr->nd.minerva.ixn[i].tracks[j].end.x;
            double start_y = sr->nd.minerva.ixn[i].tracks[j].start.y;
            double end_y = sr->nd.minerva.ixn[i].tracks[j].end.y;

            // Calculate the length and direction of the MINERvA track
            double dX = end_x - start_x;
            double dY = end_y - start_y;
            double dZ = end_z - start_z;
            double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
            double dirX = dX / length;
            double dirY = dY / length;
            double dirZ = dZ / length;

            // Skip the track if the start position is NaN or the length is less than 10 cm
            if (std::isnan(start_z) || length < 10) continue;

            // Flip the direction if dirZ is negative to ensure end_z > start_z
            if (dirZ < 0) { 
                dirZ = -dirZ; 
                dirX = -dirX; 
                dirY = -dirY; 
                std::swap(start_x, end_x);
                std::swap(start_y, end_y);
                std::swap(start_z, end_z);
            }

            // Check if the MINERvA track and TPC track meet the conditions for matching
            if (start_z > 0 && end_z > 0 && (reco_vertex.Z() > 64.538 || end_z > 64.538)) {
                double deltaExtrapY, deltaExtrapX;
                // Calculate the dot product and extrapolation values
                double dotProduct = calculateDotProductAndExtrapolation(dirX, dirY, dirZ, start_z, end_z, start_x, end_x, start_y, end_y, reco_vertex, deltaExtrapY, deltaExtrapX);

                // Update the maximum dot product and associated variables if the current track is a better match
                // Include the dot product sanity check
                if (dotProduct > maxDotProduct && std::abs(deltaExtrapY) < 10 && dotProduct > 0.9975) {
                    maxDotProduct = dotProduct;
                    minDeltaExtrapY = deltaExtrapY;
                    minDeltaExtrapX = deltaExtrapX;
                    maxPartMinerva = sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                    maxTypeMinerva = sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
                    maxIxnMinerva = sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
                }
            }
        }
    }
    
    // Return true if a valid match is found (assuming a positive dot product indicates a match)
    return maxDotProduct > 0;
}

bool isMINERvAMatch(TVector3 reco_vertex, caf::StandardRecord* sr);
void SaveHistogramAsPNG(TObject* histogram);
void DrawAndSaveHistogramWithWeightedStats_onehistogram(TH1* h_truth, const char* title, const char* fileName);
void DrawAndSaveStackedHistogram(THStack* stack, TH1D* h_QE, TH1D* h_DIS, TH1D* h_Res, TH1D* h_Coh, TH1D* h_MEC, const char* fileName);

int caf_plotter(std::string input_file_list, std::string output_rootfile) {
    
    int length_bins = 100;
    int length_bins_L = 0;
    int length_bins_U = 150;
    
    int energy_bins = 100;
    int energy_bins_L = 0;
    int energy_bins_U = 1;

    int start_bins = 50;
    int start_bins_L = -70;
    int start_bins_U = 70;

    int mult_bins_genie = 20;
    int mult_edge_L_genie = 1;
    int mult_edge_U_genie = 21;
    
    int cosTheta_bins = 100;
    int cosTheta_bins_L = -1;
    int cosTheta_bins_U = 1;
    
    int mult_bins_muon = 5;
    int mult_bins_L_muon = 1;
    int mult_bins_U_muon = 6;
    
    int mult_bins_pion = 5;
    int mult_bins_L_pion = 1;
    int mult_bins_U_pion = 6;
    
    int mult_bins_proton = 5;
    int mult_bins_L_proton = 1;
    int mult_bins_U_proton = 6;
    
    int mult_bins_total = 5;
    int mult_bins_L_total = 1;
    int mult_bins_U_total = 6;    
	
	int nu_energy_bins = 50;
	int nu_energy_bins_L = 0;
	int nu_energy_bins_U = 15;

	TH1D *true_nu_E_matched = new TH1D("true_nu_E_matched", "True Neutrino Energy (Reco-Matched);Neutrino Energy (GeV);Counts", nu_energy_bins, nu_energy_bins_L, nu_energy_bins_U);
    THStack *hs_modes_matched = new THStack("hs_modes_matched", "True Neutrino Energy (Reco-Matched);Energy (GeV);Counts");
    TH1D *h_QE_matched = new TH1D("h_QE_matched", "Quasi-Elastic;Energy (GeV);Counts", nu_energy_bins, nu_energy_bins_L, nu_energy_bins_U);
    TH1D *h_DIS_matched = new TH1D("h_DIS_matched", "Deep Inelastic Scattering;Energy (GeV);Counts", nu_energy_bins, nu_energy_bins_L, nu_energy_bins_U);
    TH1D *h_Res_matched = new TH1D("h_Res_matched", "Resonant;Energy (GeV);Counts", nu_energy_bins, nu_energy_bins_L, nu_energy_bins_U);
    TH1D *h_Coh_matched = new TH1D("h_Coh_matched", "Coherent;Energy (GeV);Counts", nu_energy_bins, nu_energy_bins_L, nu_energy_bins_U);
    TH1D *h_MEC_matched = new TH1D("h_MEC_matched", "Meson Exchange Currents;Energy (GeV);Counts", nu_energy_bins, nu_energy_bins_L, nu_energy_bins_U);

    h_QE_matched->SetFillColor(kRed);
    h_DIS_matched->SetFillColor(kBlue);
    h_Res_matched->SetFillColor(kGreen);
    h_Coh_matched->SetFillColor(kYellow);
    h_MEC_matched->SetFillColor(kMagenta);

    hs_modes_matched->Add(h_QE_matched);
    hs_modes_matched->Add(h_DIS_matched);
    hs_modes_matched->Add(h_Res_matched);
    hs_modes_matched->Add(h_Coh_matched);
    hs_modes_matched->Add(h_MEC_matched);

    TH1D *true_nu_vtx_x = new TH1D("true_nu_vtx_x", "Neutrino Vertex-x;Vertex-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *true_nu_vtx_y = new TH1D("true_nu_vtx_y", "Neutrino Vertex-y;Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *true_nu_vtx_z = new TH1D("true_nu_vtx_z", "Neutrino Vertex-z;Vertex-z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1D *reco_nu_vtx_x = new TH1D("reco_nu_vtx_x", "Reconstructed Neutrino Vertex-x;Vertex-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_y = new TH1D("reco_nu_vtx_y", "Reconstructed Neutrino Vertex-y;Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_z = new TH1D("reco_nu_vtx_z", "Reconstructed Neutrino Vertex-z;Vertex-z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1D *true_genie_multiplicity = new TH1D("true_genie_multiplicity", "True GENIE Multplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);
    TH1D *true_trackOnly_multiplicity = new TH1D("true_trackOnly_multiplicity", "True Track Multplicity;Multiplicity;Counts", mult_bins_total, mult_bins_L_total, mult_bins_U_total);

    TH1D *reco_mult_prim_muon = new TH1D("reco_mult_prim_muon", "Reconstructed Muons;Number of Tracks;Counts", mult_bins_muon, mult_bins_L_muon, mult_bins_U_muon);
    TH1D *reco_mult_prim_pion = new TH1D("reco_mult_prim_pion", "Reconstructed Charged Pions;Number of Tracks;Counts", mult_bins_pion, mult_bins_L_pion, mult_bins_U_pion);
    TH1D *reco_mult_prim_proton = new TH1D("reco_mult_prim_proton", "Reconstructed Protons;Number of Tracks;Counts", mult_bins_proton, mult_bins_L_proton, mult_bins_U_proton);
    TH1D *reco_mult_prim_total = new TH1D("reco_mult_prim_total", "Reconstructed Track Multplicity;Multiplicity;Counts", mult_bins_total, mult_bins_L_total, mult_bins_U_total);
    
    TH1F *reco_length_prim_muon = new TH1F("reco_length_prim_muon", "Reconstructed Muons;Track Length;Counts", length_bins, length_bins_L, length_bins_U);
    TH1F *reco_length_prim_pion = new TH1F("reco_length_prim_pion", "Reconstructed Charged Pions;Track Length;Counts", length_bins, length_bins_L, length_bins_U);
    TH1F *reco_length_prim_proton = new TH1F("reco_length_prim_proton", "Reconstructed Protons;Track Length;Counts", length_bins, length_bins_L, length_bins_U);

    TH1F *reco_cosTheta_prim_muon = new TH1F("reco_cosTheta_prim_muon", "Reconstructed Muons;Cosine Theta;Counts", cosTheta_bins, cosTheta_bins_L, cosTheta_bins_U);
    TH1F *reco_cosTheta_prim_pion = new TH1F("reco_cosTheta_prim_pion", "Reconstructed Charged Pions;Cosine Theta;Counts", cosTheta_bins, cosTheta_bins_L, cosTheta_bins_U);
    TH1F *reco_cosTheta_prim_proton = new TH1F("reco_cosTheta_prim_proton", "Reconstructed Protons;Cosine Theta;Counts", cosTheta_bins, cosTheta_bins_L, cosTheta_bins_U);

    TH1F *reco_startX_prim_muon = new TH1F("reco_startX_prim_muon", "Reconstructed Muons;Start-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1F *reco_startY_prim_muon = new TH1F("reco_startY_prim_muon", "Reconstructed Muons;Start-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1F *reco_startZ_prim_muon = new TH1F("reco_startZ_prim_muon", "Reconstructed Muons;Start-z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1F *reco_startX_prim_pion = new TH1F("reco_startX_prim_pion", "Reconstructed Charged Pions;Start-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1F *reco_startY_prim_pion = new TH1F("reco_startY_prim_pion", "Reconstructed Charged Pions;Start-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1F *reco_startZ_prim_pion = new TH1F("reco_startZ_prim_pion", "Reconstructed Charged Pions;Start-z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1F *reco_startX_prim_proton = new TH1F("reco_startX_prim_proton", "Reconstructed Protons;Start-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1F *reco_startY_prim_proton = new TH1F("reco_startY_prim_proton", "Reconstructed Protons;Start-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1F *reco_startZ_prim_proton = new TH1F("reco_startZ_prim_proton", "Reconstructed Protons;Start-z (cm);Counts", start_bins, start_bins_L, start_bins_U);
    
    TH1F *reco_energy_prim_pion = new TH1F("reco_energy_prim_pion", "Reconstructed Charged Pions;Kinetic Energy (GeV);Counts", energy_bins, energy_bins_L, energy_bins_U);
    TH1F *reco_energy_prim_proton = new TH1F("reco_energy_prim_proton", "Reconstructed Protons;Kinetic Energy (GeV);Counts", energy_bins, energy_bins_L, energy_bins_U);

    // Counters for the print statements
    int totalInteractions = 0;
    int interactionsInRock = 0;
    int ccInteractions = 0;
    int interactionsWithLAr = 0;
    int interactionsMatchedToMINERvA = 0;
    int interactionsMatchedWithinLArFV = 0;
    int interactionsMatchedToTruthWithMuon = 0;

    // Give an input list
    std::ifstream caf_list(input_file_list.c_str());

    // Check if input list is present
    if (!caf_list.is_open()) {
        std::cerr << Form("File %s not found", input_file_list.c_str()) << std::endl;
        return 1;
    }

    // Add files to CAF chain from input list
    std::string tmp;
    TChain *caf_chain = new TChain("cafTree");

    while (caf_list >> tmp) {
        caf_chain->Add(tmp.c_str());
        std::cout << Form("Adding File %s", tmp.c_str()) << std::endl;
    }

    // Check if CAF tree is present
    if (!caf_chain) {
        std::cerr << Form("There is no tree in %s", tmp.c_str()) << std::endl;
        return 1;
    }

    long Nentries = caf_chain->GetEntries();
    std::cout << Form("Total number of spills = %ld", Nentries) << std::endl;

    // Define Standard Record and link it to the CAF tree branch "rec"
    auto sr = new caf::StandardRecord;
    caf_chain->SetBranchAddress("rec", &sr);

    // Bounds on outer boundaries
    const double distance_from_wall = 5.0;
    const double minX = -63.931 + distance_from_wall;
    const double maxX = +63.931 - distance_from_wall;
    const double minY = -62.076 + distance_from_wall;
    const double maxY = +62.076 - distance_from_wall;
    const double minZ = -64.538 + distance_from_wall;
    const double maxZ = +64.538 - distance_from_wall;
    const double general_track_length = 5;
    
    // Bounds on inner boundaries
    const double X_average = (minX + maxX) / 2;
    const double Z_average = (minZ + maxZ) / 2;
    const double inner_epsilon = 5.0;
    const double Xinner_neg = X_average - inner_epsilon;
    const double Xinner_pos = X_average + inner_epsilon;
    const double Zinner_neg = Z_average - inner_epsilon;
    const double Zinner_pos = Z_average + inner_epsilon;

    // Right before entering the loop over spills/triggers
    int entryCounter = 1;
    
    // Loop over spills/triggers
    for (long n = 0; n < Nentries; ++n) { 
    
        if (n % 10000 == 0) std::cout << Form("Processing trigger %ld of %ld", n, Nentries) << std::endl;
        caf_chain->GetEntry(n); // Get spill from tree

        // Loop over truth interactions
        for (long unsigned ntrue = 0; ntrue < sr->mc.nu.size(); ntrue++) {
            auto vertex = sr->mc.nu[ntrue].vtx;
            double true_nu_vtxX = vertex.x;
            double true_nu_vtxY = vertex.y;
            double true_nu_vtxZ = vertex.z;

            auto truePrimary = sr->mc.nu[ntrue].prim;

            int trueGenieMultiplicity = 0;         // for true_genie_multiplicity
            int trueTrackMultiplicity = 0;     // true_trackOnly_multiplicity

            if (sr->mc.nu[ntrue].targetPDG != 1000180400) continue;
            if (sr->mc.nu[ntrue].iscc == false) continue;

            for (long unsigned primaries = 0; primaries < truePrimary.size(); primaries++) {
                int pdg = truePrimary[primaries].pdg;
                auto start_pos_true = sr->mc.nu[ntrue].prim[primaries].start_pos;
                auto end_pos_true = sr->mc.nu[ntrue].prim[primaries].end_pos;

                double dX_true = (end_pos_true.x - start_pos_true.x);
                double dY_true = (end_pos_true.y - start_pos_true.y);
                double dZ_true = (end_pos_true.z - start_pos_true.z);
                double length_true = TMath::Sqrt(dX_true * dX_true + dY_true * dY_true + dZ_true * dZ_true);
                double cosTheta_true = dZ_true / length_true;    

                double currentX_prim_true = start_pos_true.x;
                double currentY_prim_true = start_pos_true.y;
                double currentZ_prim_true = start_pos_true.z;

                bool isXOutOfRange_true = abs(currentX_prim_true) > maxX || abs(currentX_prim_true) < minX;
                bool isYOutOfRange_true = abs(currentY_prim_true) < minY || abs(currentY_prim_true) > maxY;
                bool isZOutOfRange_true = abs(currentZ_prim_true) < minZ || abs(currentZ_prim_true) > maxZ;

                // Check if true_nu_vtxX, true_nu_vtxY, true_nu_vtxZ are within the inner exclusion zone (out of range if inside)
                bool isXInRange_inner_true = true_nu_vtxX > Xinner_neg && true_nu_vtxX < Xinner_pos;
                bool isZInRange_inner_true = true_nu_vtxZ > Zinner_neg && true_nu_vtxZ < Zinner_pos;

                if (isXOutOfRange_true || isYOutOfRange_true || isZOutOfRange_true || isXInRange_inner_true || isZInRange_inner_true) {
                    continue; // Skip to the next iteration if any condition is true
                }

                if ((abs(pdg) == 2212 || pdg == 211) && length_true > general_track_length) {
                    trueTrackMultiplicity++;
                }

                trueGenieMultiplicity++;
            }

            true_genie_multiplicity->Fill(trueGenieMultiplicity);
            true_trackOnly_multiplicity->Fill(trueTrackMultiplicity);
        }

        // Loop over reconstructed interactions
        for (long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++) {
            totalInteractions++; // Increment total number of interactions

            // declare some variables to find the best match for the current interaction based on truth information
            double biggestMatch = -999;  
            int biggestMatchIndex = -999;

            // Variables to track if there's a muon in reco and truth interactions
            bool hasRecoMuon = false;
            bool hasTruthMuon = false;

            // This loop iterates over all "truth" entries associated with the nth interaction. 
            // Inside the loop, it checks if the current "truth" entry's overlap (measured by sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)) is larger than the largest overlap found so far (biggestMatch). 
            // If the condition is true, the code updates biggestMatch with the current "truth" entry's overlap value, indicating this is the largest overlap found so far.
            // It also records the index of this "truth" entry (ntruth) by storing it in biggestMatchIndex.         
            for (int ntruth = 0; ntruth < sr->common.ixn.dlp[nixn].truth.size(); ntruth++) {
                if (biggestMatch < sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)) {
                    biggestMatch = sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth);
                    biggestMatchIndex = sr->common.ixn.dlp[nixn].truth.at(ntruth);
                }
            }

            // Continue if there's no match
            if (biggestMatchIndex == -999) continue;

            // Checks if the neutrino ID at the biggestMatchIndex is greater than 1E9. If this condition is met, it sets the variable rock to 1.
            if (sr->mc.nu[biggestMatchIndex].id > 1E9) {
                interactionsInRock++; // Increment interactions in the rock
            }

            double reco_nu_vtxX = sr->common.ixn.dlp[nixn].vtx.x;
            double reco_nu_vtxY = sr->common.ixn.dlp[nixn].vtx.y;
            double reco_nu_vtxZ = sr->common.ixn.dlp[nixn].vtx.z;

            double true_nu_vtxX = sr->mc.nu[biggestMatchIndex].vtx.x;
            double true_nu_vtxY = sr->mc.nu[biggestMatchIndex].vtx.y;
            double true_nu_vtxZ = sr->mc.nu[biggestMatchIndex].vtx.z;
            bool   iscc         = sr->mc.nu[biggestMatchIndex].iscc;
            double targetPDG    = sr->mc.nu[biggestMatchIndex].targetPDG;
			double energy 		= sr->mc.nu[biggestMatchIndex].E; 
			double mode 		= sr->mc.nu[biggestMatchIndex].mode; 

            // Check if true neutrino vertices are within FV bounds
            bool isXOutOfRange_true = abs(true_nu_vtxX) < minX || abs(true_nu_vtxX) > maxX;
            bool isYOutOfRange_true = abs(true_nu_vtxY) < minY || abs(true_nu_vtxY) > maxY;
            bool isZOutOfRange_true = abs(true_nu_vtxZ) < minZ || abs(true_nu_vtxZ) > maxZ;
            bool isXInRange_true = true_nu_vtxX > Xinner_neg && true_nu_vtxX < Xinner_pos;
            bool isZInRange_true = true_nu_vtxZ > Zinner_neg && true_nu_vtxZ < Zinner_pos;

            // Check if reconstructed neutrino vertices are within FV bounds
            bool isXOutOfRange_reco = abs(reco_nu_vtxX) < minX || abs(reco_nu_vtxX) > maxX;          
            bool isYOutOfRange_reco = abs(reco_nu_vtxY) < minY || abs(reco_nu_vtxY) > maxY;
            bool isZOutOfRange_reco = abs(reco_nu_vtxZ) < minZ || abs(reco_nu_vtxZ) > maxZ;
            bool isXInRange_reco = reco_nu_vtxX > Xinner_neg && reco_nu_vtxX < Xinner_pos;
            bool isZInRange_reco = reco_nu_vtxZ > Zinner_neg && reco_nu_vtxZ < Zinner_pos;

            if (iscc != 1) continue; // Skip if not a CC interaction 
            ccInteractions++; // Increment CC interactions

            if (targetPDG != 1000180400) continue; // Skip if not interacting with Argon
            interactionsWithLAr++; // Increment interactions with LAr

            // MINERvA matching
            TVector3 reco_vertex = TVector3(reco_nu_vtxX, reco_nu_vtxY, reco_nu_vtxZ);
            if (isMINERvAMatch(reco_vertex, sr)) {
                interactionsMatchedToMINERvA++; // Increment interactions matched to MINERvA
                continue;
            }

            // If the true and reco points are out of outer bounds OR within the inner exclusion zone
            if (isXOutOfRange_true || isYOutOfRange_true || isZOutOfRange_true || isXInRange_true || isZInRange_true) {
                continue; 
            }

            if (isXOutOfRange_reco || isYOutOfRange_reco || isZOutOfRange_reco || isXInRange_reco || isZInRange_reco) {
                continue;
            }
			
            interactionsMatchedWithinLArFV++; // Increment interactions matched within LAr FV

            // Check for a reco muon in the current interaction
            int reco_muonCount = 0; 

            for (long unsigned npart = 0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++) {
                auto& particle = sr->common.ixn.dlp[nixn].part.dlp[npart];
                auto start_pos_reco = particle.start;
                auto end_pos_reco = particle.end;       

                double dX_reco = (end_pos_reco.x - start_pos_reco.x);
                double dY_reco = (end_pos_reco.y - start_pos_reco.y);
                double dZ_reco = (end_pos_reco.z - start_pos_reco.z);
                double length_reco = TMath::Sqrt(dX_reco * dX_reco + dY_reco * dY_reco + dZ_reco * dZ_reco);
                double cosTheta_reco = dZ_reco / length_reco;

                double currentX_reco = start_pos_reco.x;
                double currentY_reco = start_pos_reco.y;
                double currentZ_reco = start_pos_reco.z;

                // Check if starting point of reconstructed particle is within FV
                bool isXOutOfRange_reco = abs(currentX_reco) < minX || abs(currentX_reco) > maxX;
                bool isYOutOfRange_reco = abs(currentY_reco) < minY || abs(currentY_reco) > maxY;
                bool isZOutOfRange_reco = abs(currentZ_reco) < minZ || abs(currentZ_reco) > maxZ;
                bool isXInRange_reco = currentX_reco > Xinner_neg && currentX_reco < Xinner_pos;
                bool isZInRange_reco = currentZ_reco > Zinner_neg && currentZ_reco < Zinner_pos;

                // Skip particle if it is out of the specified range
                if (isXOutOfRange_reco || isYOutOfRange_reco || isZOutOfRange_reco || isXInRange_reco || isZInRange_reco) {
                    continue;
                }

                // If within bounds, check if particle is a muon and the track length exceeds the general threshold
                if (abs(particle.pdg) == 13 && length_reco > general_track_length) {
                    reco_muonCount++; 
                    if (reco_muonCount > 1) {
                        break; // Exit the loop if more than one muon is found
                    }
                    hasRecoMuon = true;
                }
            }

            // Skip the current interaction if there is more than one muon
            if (reco_muonCount > 1) {
                continue;
            }

            // Check for a truth muon in the current interaction         
            if (biggestMatchIndex != -999) { // Ensure there is a match
                for (int k = 0; k < sr->mc.nu[biggestMatchIndex].prim.size(); k++) {
                    if (sr->mc.nu[biggestMatchIndex].prim[k].interaction_id == sr->mc.nu[biggestMatchIndex].id) {
                        if (abs(sr->mc.nu[biggestMatchIndex].prim[k].pdg) == 13) {
                            hasTruthMuon = true;
                        }
                    }
                }
            }

            if (!hasTruthMuon) continue;
            interactionsMatchedToTruthWithMuon++; // Increment interactions matched to truth and having a truth muon
			
			true_nu_E_matched->Fill(energy);
			// Fill histograms based on mode
            switch (static_cast<int>(mode)) {
                case 1: h_QE_matched->Fill(energy); break; // QE
                case 3: h_DIS_matched->Fill(energy); break; // DIS
                case 4: h_Res_matched->Fill(energy); break; // Res
                case 5: h_Coh_matched->Fill(energy); break; // Coh
                case 10: h_MEC_matched->Fill(energy); break; // MEC
                default: break;
            }

            true_nu_vtx_x->Fill(true_nu_vtxX);
            true_nu_vtx_y->Fill(true_nu_vtxY);
            true_nu_vtx_z->Fill(true_nu_vtxZ);

            reco_nu_vtx_x->Fill(reco_nu_vtxX);
            reco_nu_vtx_y->Fill(reco_nu_vtxY);
            reco_nu_vtx_z->Fill(reco_nu_vtxZ);

            int recoPrimary_total = 1;
            int recoPrimary_muon = 0;
            int recoPrimary_pion = 1;
            int recoPrimary_proton = 1;

            // Loop over reconstructed particles
            for (long unsigned npart = 0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++) {
                auto& particle = sr->common.ixn.dlp[nixn].part.dlp[npart];
                auto start_pos_reco = particle.start;
                auto end_pos_reco = particle.end;

                double dX_reco = (end_pos_reco.x - start_pos_reco.x);
                double dY_reco = (end_pos_reco.y - start_pos_reco.y);
                double dZ_reco = (end_pos_reco.z - start_pos_reco.z);
                double length_reco = TMath::Sqrt(dX_reco * dX_reco + dY_reco * dY_reco + dZ_reco * dZ_reco);
                double cosTheta_reco = dZ_reco / length_reco;

                double currentX_reco = start_pos_reco.x;
                double currentY_reco = start_pos_reco.y;
                double currentZ_reco = start_pos_reco.z;

                bool isXOutOfRange_reco = abs(currentX_reco) < minX || abs(currentX_reco) > maxX;
                bool isYOutOfRange_reco = abs(currentY_reco) < minY || abs(currentY_reco) > maxY;
                bool isZOutOfRange_reco = abs(currentZ_reco) < minZ || abs(currentZ_reco) > maxZ;
                bool isXInRange_reco = currentX_reco > Xinner_neg && currentX_reco < Xinner_pos;
                bool isZInRange_reco = currentZ_reco > Zinner_neg && currentZ_reco < Zinner_pos;

                if (isXOutOfRange_reco || isYOutOfRange_reco || isZOutOfRange_reco || isXInRange_reco || isZInRange_reco) {
                    continue;
                }

                // Common logic for counting and filling histograms
                bool isMuon = abs(particle.pdg) == 13;
                bool isPion = abs(particle.pdg) == 211;
                bool isProton = abs(particle.pdg) == 2212;
                bool validLength = length_reco > general_track_length;

                if (isMuon && validLength) {
                    reco_length_prim_muon->Fill(length_reco);
                    reco_cosTheta_prim_muon->Fill(cosTheta_reco);
                    reco_startX_prim_muon->Fill(start_pos_reco.x);
                    reco_startY_prim_muon->Fill(start_pos_reco.y);
                    reco_startZ_prim_muon->Fill(start_pos_reco.z);
                    recoPrimary_muon++;
                }
                if (isPion && validLength) {
                    reco_length_prim_pion->Fill(length_reco);
                    reco_cosTheta_prim_pion->Fill(cosTheta_reco);
                    reco_startX_prim_pion->Fill(start_pos_reco.x);
                    reco_startY_prim_pion->Fill(start_pos_reco.y);
                    reco_startZ_prim_pion->Fill(start_pos_reco.z);
                    reco_energy_prim_pion->Fill(particle.E);
                    recoPrimary_pion++;
                }
                if (isProton && validLength) {
                    reco_length_prim_proton->Fill(length_reco);
                    reco_cosTheta_prim_proton->Fill(cosTheta_reco);
                    reco_startX_prim_proton->Fill(start_pos_reco.x);
                    reco_startY_prim_proton->Fill(start_pos_reco.y);
                    reco_startZ_prim_proton->Fill(start_pos_reco.z);
                    reco_energy_prim_proton->Fill(particle.E);
                    recoPrimary_proton++;
                }

                if ((isProton || isPion) && validLength) {
                    recoPrimary_total++;
                }
            }

            reco_mult_prim_muon->Fill(recoPrimary_muon);
            reco_mult_prim_pion->Fill(recoPrimary_pion);
            reco_mult_prim_proton->Fill(recoPrimary_proton);
            reco_mult_prim_total->Fill(recoPrimary_total);

        }
        entryCounter++;
    }

    // Print the final counts
    std::cout << "Total number of interactions: " << totalInteractions << std::endl;
    std::cout << "Total number of interactions in the rock: " << interactionsInRock << std::endl;
    std::cout << "Total number of CC interactions: " << ccInteractions << std::endl;
    std::cout << "Total number of interactions with LAr: " << interactionsWithLAr << std::endl;
    std::cout << "Total number of interactions matched to MINERvA: " << interactionsMatchedToMINERvA << std::endl;
    std::cout << "Total number of interactions matched within LAr FV: " << interactionsMatchedWithinLArFV << std::endl;
    std::cout << "Total number of interactions matched to truth and having a truth muon: " << interactionsMatchedToTruthWithMuon << std::endl;

    // Save histograms
    TFile* outputFile = new TFile(output_rootfile.c_str(), "RECREATE");
    true_nu_vtx_x->Write();
    true_nu_vtx_y->Write();
    true_nu_vtx_z->Write();
    reco_nu_vtx_x->Write();
    reco_nu_vtx_y->Write();
    reco_nu_vtx_z->Write();
    true_genie_multiplicity->Write();
    true_trackOnly_multiplicity->Write();
    reco_mult_prim_muon->Write();
    reco_mult_prim_pion->Write();
    reco_mult_prim_proton->Write();
    reco_mult_prim_total->Write();
    reco_length_prim_muon->Write();
    reco_length_prim_pion->Write();
    reco_length_prim_proton->Write();
    reco_cosTheta_prim_muon->Write();
    reco_cosTheta_prim_pion->Write();
    reco_cosTheta_prim_proton->Write();
    reco_startX_prim_muon->Write();
    reco_startY_prim_muon->Write();
    reco_startZ_prim_muon->Write();
    reco_startX_prim_pion->Write();
    reco_startY_prim_pion->Write();
    reco_startZ_prim_pion->Write();
    reco_startX_prim_proton->Write();
    reco_startY_prim_proton->Write();
    reco_startZ_prim_proton->Write();
    reco_energy_prim_pion->Write();
    reco_energy_prim_proton->Write();
	
	DrawAndSaveHistogramWithWeightedStats_onehistogram(true_genie_multiplicity,"Multiplicity;Number of Charged Tracks;Counts", "true_genie_multiplicity.png");    
    DrawAndSaveHistogramWithWeightedStats_onehistogram(true_trackOnly_multiplicity,"Multiplicity;Number of Charged Tracks;Counts", "true_trackOnly_multiplicity.png");    

    DrawAndSaveHistogramWithWeightedStats_onehistogram(reco_mult_prim_total,"Multiplicity;Number of Charged Tracks;Counts", "reco_mult_total.png");    
    DrawAndSaveHistogramWithWeightedStats_onehistogram(reco_mult_prim_pion, "Charged Pions;Number of Tracks;Counts", "reco_mult_pion.png");    
    DrawAndSaveHistogramWithWeightedStats_onehistogram(reco_mult_prim_muon, "Muons;Number of Tracks;Counts", "reco_mult_muon.png");    
    DrawAndSaveHistogramWithWeightedStats_onehistogram(reco_mult_prim_proton, "Protons;Number of Tracks;Counts", "reco_mult_proton.png");

    SaveHistogramAsPNG(reco_length_prim_pion);    
    SaveHistogramAsPNG(reco_length_prim_proton);    
    SaveHistogramAsPNG(reco_cosTheta_prim_pion);    
    SaveHistogramAsPNG(reco_cosTheta_prim_proton);    
    SaveHistogramAsPNG(reco_cosTheta_prim_muon);    
    SaveHistogramAsPNG(reco_startX_prim_muon);    
    SaveHistogramAsPNG(reco_startY_prim_muon);    
    SaveHistogramAsPNG(reco_startZ_prim_muon);    
    SaveHistogramAsPNG(reco_startX_prim_pion);    
    SaveHistogramAsPNG(reco_startY_prim_pion);    
    SaveHistogramAsPNG(reco_startZ_prim_pion);
    SaveHistogramAsPNG(reco_startX_prim_proton);    
    SaveHistogramAsPNG(reco_startY_prim_proton);    
    SaveHistogramAsPNG(reco_startZ_prim_proton);

    SaveHistogramAsPNG(reco_energy_prim_pion);
    SaveHistogramAsPNG(reco_energy_prim_proton);

    SaveHistogramAsPNG(reco_nu_vtx_x);
    SaveHistogramAsPNG(reco_nu_vtx_y);
    SaveHistogramAsPNG(reco_nu_vtx_z);
	
	true_nu_E_matched->Write();
	h_QE_matched->Write();
    h_DIS_matched->Write();
    h_Res_matched->Write();
    h_Coh_matched->Write();
    h_MEC_matched->Write();
    hs_modes_matched->Write();

    // Draw and save the stacked histogram
    DrawAndSaveStackedHistogram(hs_modes_matched, h_QE_matched, h_DIS_matched, h_Res_matched, h_Coh_matched, h_MEC_matched, "stacked_energy_histogram_matched.png");

    outputFile->Close();    
    return 1;
}   // caf_plotter(std::string input_file_list, std::string output_rootfile)


int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "\n USAGE: " << argv[0] << " input_caf_file_list output_root_file\n" << std::endl;
        return 1;
    }

    std::string input_file_list = argv[1];
    std::string output_rootfile = argv[2];

    caf_plotter(input_file_list, output_rootfile);

    return 0;
}

void SaveHistogramAsPNG(TObject* histogram) {
    if (!histogram) return; // If the histogram pointer is null, exit the function

    if (histogram->InheritsFrom(TH1::Class())) {
        TH1* hist = static_cast<TH1*>(histogram);

        hist->SetLineColor(kRed); 
        hist->SetLineWidth(3);
        hist->SetLineStyle(1);

        // Create a canvas to draw the histogram
        TCanvas canvas("canvas", "Canvas", 800, 600);
        gStyle->SetOptStat("emr"); // Set the option for the statistics box

        hist->Draw();

        // Save the canvas as a PNG file with the name of the histogram
        std::string fileName = std::string(hist->GetName()) + ".png";
        canvas.SaveAs(fileName.c_str());
    } else {
        std::cerr << "The provided object is not a histogram." << std::endl;
    }
}

void DrawAndSaveHistogramWithWeightedStats_onehistogram(TH1* h_truth, const char* title, const char* fileName) {
    if (!h_truth) return; // Check if h_truth is nullptr

    // Turn off the default stat box
    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c", title, 800, 600);
    c->SetMargin(0.1, 0.02, 0.1, 0.1); // Adjust canvas margins

    // Function to calculate custom statistics
    auto calculateCustomStats = [](TH1* h, double& customMean, double& customStdDev, int& totalEntries) {
        double weightedSum = 0;
        double weightedSumSq = 0;
        totalEntries = 0;

        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double binContent = h->GetBinContent(i);
            double binCenter = h->GetBinCenter(i);
            weightedSum += binCenter * binContent;
            weightedSumSq += binCenter * binCenter * binContent;
            totalEntries += binContent;
        }

        customMean = (totalEntries > 0) ? weightedSum / totalEntries : 0;
        double variance = (totalEntries > 0) ? (weightedSumSq / totalEntries) - (customMean * customMean) : 0;
        customStdDev = TMath::Sqrt(variance);
    };

    // Calculate custom statistics for the truth histogram
    double customMeanTruth, customStdDevTruth;
    int totalEntriesTruth;

    calculateCustomStats(h_truth, customMeanTruth, customStdDevTruth, totalEntriesTruth);

    // Prepare the histogram
    h_truth->SetLineColor(kRed);
    h_truth->SetLineWidth(3);
    h_truth->SetLineStyle(1); // Solid line

    // Find the maximum y value to adjust the scale
    double maxY = h_truth->GetMaximum();
    h_truth->SetMaximum(maxY * 1.2); // Scale up by x% for some extra space

    // Draw histogram
    h_truth->Draw("HIST");

    // Prepare legend text including custom statistics
    char legendTextTruth[256];
    snprintf(legendTextTruth, sizeof(legendTextTruth), "Reco:  %d, %.2f, %.2f", totalEntriesTruth, customMeanTruth, customStdDevTruth);

    // Create and fill the legend
    TLegend* legend = new TLegend(0.65, 0.75, 0.97, 0.85);
    legend->SetTextSize(0.03);
    legend->AddEntry(h_truth, legendTextTruth, "l");
    legend->Draw();

    // Save canvas
    c->SaveAs(fileName);

    // Cleanup
    delete c;
    delete legend;
}

void DrawAndSaveStackedHistogram(THStack* stack, TH1D* h_QE, TH1D* h_DIS, TH1D* h_Res, TH1D* h_Coh, TH1D* h_MEC, const char* fileName) {
    if (!stack || !h_QE || !h_DIS || !h_Res || !h_Coh || !h_MEC) return; // Check if any histogram is nullptr

    // Create a canvas to draw the histogram stack
    TCanvas* c = new TCanvas("c_modes_matched", "True Neutrino Energy", 800, 600);

    // Draw the histogram stack
    stack->Draw("hist");

    // Create a legend and add entries
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h_QE, "QE", "f");
    legend->AddEntry(h_DIS, "DIS", "f");
    legend->AddEntry(h_Res, "RES", "f");
    legend->AddEntry(h_Coh, "COH", "f");
    legend->AddEntry(h_MEC, "MEC", "f");
    legend->Draw();

    // Save the canvas as a PNG file
    c->SaveAs(fileName);

    // Cleanup
    delete c;
    delete legend;
}
