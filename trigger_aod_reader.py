import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Open the ROOT file
file = ROOT.TFile.Open("AOD_trig_2000.pool.root.1.root")

# Access the TTree
tree = file.Get("analysis")
bin_edges = np.array([0.0, 30.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0, 300.0, 350.0, 500.0])
n_bins = len(bin_edges) - 1

# Create histograms for the numerator and denominator
h_num = ROOT.TH1F("h_num", "Numerator;llp_pt;Entries", n_bins, bin_edges)
h_den = ROOT.TH1F("h_den", "Denominator;llp_pt;Entries", n_bins, bin_edges)
# Define the numerator cut
cut_Lxy_num = "(llp_Lxy[HLT_jet_minDR_LLPindex[0]] > 2000 && llp_Lxy[HLT_jet_minDR_LLPindex[0]] < 3800)"
cut_pt_num = "(llp_pt[HLT_jet_minDR_LLPindex[0]] > 30)"
cut_eta_num = "(abs(llp_eta[HLT_jet_minDR_LLPindex[0]]) < 1.4)"
cut_minDR_num ="(HLT_jet_minDR_LLP < 0.2)"

cut_numerator = cut_Lxy_num+ " && " + cut_pt_num + " && "+cut_eta_num + " && " + cut_minDR_num
# Define the denominator cut
cut_denominator = "llp_Lxy > 2000 && llp_Lxy < 3800 && llp_pt > 30 && abs(llp_eta) < 1.4"

# Draw the numerator histogram
tree.Draw("llp_pt[HLT_jet_minDR_LLPindex[0]] >> h_num", cut_numerator)
# Draw the denominator histogram
tree.Draw("llp_pt >> h_den", cut_denominator)
import uproot
#convert TH1F to numpy
num_values=uproot.from_pyroot(h_num).to_numpy()
den_values=uproot.from_pyroot(h_den).to_numpy()

# Calculate efficiency and binomial errors
efficiency = np.divide(num_values[0], den_values[0], out=np.zeros_like(num_values[0]), where=den_values[0] != 0)
efficiency_error = np.sqrt(efficiency * (1 - efficiency) / den_values[0])

plt.figure(figsize=(10, 6))
plt.hist(bin_edges[:-1], bins=bin_edges, weights=efficiency, histtype='step', label="Efficiency", color="green")
# Customize and show the efficiency plot
plt.xlabel("llp_pt")
plt.ylabel("Efficiency")
plt.legend()
# Plot the efficiency with binomial errors
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Calculate bin centers for plotting
plt.errorbar(bin_centers, efficiency, yerr=efficiency_error, fmt='o', color="green", label="Efficiency", capsize=5)
plt.xlabel("llp_pt")
plt.ylabel("Efficiency")
plt.legend()

# Save the efficiency plot to a file
plt.savefig("efficiency_emf_008.png")
