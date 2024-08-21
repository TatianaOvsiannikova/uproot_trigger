import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
# Open the ROOT file and access the tree
file2 = uproot.open("/data/AOD_trig_2000.pool.root.1.root")
tree2 = file2["analysis"]
# Read the branches into an Awkward Array
arrays = tree2.arrays(["llp_pt", "llp_Lxy", "llp_eta", "HLT_jet_minDR_LLPindex", "HLT_jet_minDR_LLP"], library="ak")

# Create a mask to check if HLT_jet_minDR_LLPindex is non-empty
non_empty_index_mask = ak.num(arrays["HLT_jet_minDR_LLPindex"][arrays["HLT_jet_minDR_LLP"] < 0.2]) > 0

# Filter arrays where HLT_jet_minDR_LLPindex is non-empty
filtered_arrays = arrays[non_empty_index_mask]

# Extract the indices from HLT_jet_minDR_LLPindex
indices = ak.to_numpy(filtered_arrays["HLT_jet_minDR_LLPindex"][:, 0])
# Apply the Lxy cuts 
llppt_cut_Lxy_num_up = ak.Array([
     (lxy[idx] > 2000) & (lxy[idx] < 3800)
     for lxy, idx in zip(filtered_arrays["llp_Lxy"], indices) ])
     # Gordon: This looks like loop on each Lxy - why not used index slicing?
     # this means `filtered_arrays.llp_Lxy[indicies]`
     # would that work?

# Apply the pt cut 
llppt_cut_pt_num_up = ak.Array([
  pt[idx] > 30
 for pt, idx in zip(filtered_arrays["llp_pt"], indices) ])

# Apply the eta cut 
llppt_cut_eta_num_up = ak.Array([
 abs(eta[idx]) < 1.4
 for eta, idx in zip(filtered_arrays["llp_eta"], indices) ])

# Combine all cuts
mask_num = llppt_cut_Lxy_num_up & llppt_cut_pt_num_up & llppt_cut_eta_num_up
# Filter the original arrays based on the combined cut
llp_pt_num = filtered_arrays["llp_pt"][mask_num]
# Extract the indices from HLT_jet_minDR_LLPindex
indices2= ak.to_numpy(filtered_arrays["HLT_jet_minDR_LLPindex"][mask_num][:, 0])

# flatten array
llp_pt_num_fl = ak.Array([pt[idx]
 for pt, idx in zip( llp_pt_num, indices2) ])

cut_Lxy_den = (arrays["llp_Lxy"] > 2000) & (arrays["llp_Lxy"] < 3800)
cut_pt_den = arrays["llp_pt"] > 30
cut_eta_den = abs(arrays["llp_eta"]) < 1.4

# Combine all denominator cuts
mask_den = cut_Lxy_den & cut_pt_den & cut_eta_den

# Apply the denominator mask 
llp_pt_den = arrays["llp_pt"][mask_den]

non_empty_index_den = ak.num(llp_pt_den) > 0
# flatten array
llp_pt_den_fl= ak.flatten(llp_pt_den[non_empty_index_den], axis=1)
# binning
bin_edges = np.array([0.0, 30.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0, 300.0, 350.0, 500.0])
n_bins = len(bin_edges) - 1
# Histogram for the nominator
hist_num, _ = np.histogram(ak.to_numpy(llp_pt_num_fl), bins=bin_edges)
# Histogram for the denominator
hist_den, _ = np.histogram(ak.to_numpy(llp_pt_den_fl), bins=bin_edges)

# Calculate efficiency and binomial errors
efficiency = np.divide(hist_num, hist_den, out=np.zeros_like(hist_num, dtype=float), where=hist_den != 0)
efficiency_error = np.sqrt(efficiency * (1 - efficiency) / hist_den)

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
plt.savefig("efficiency_emf_008_up.png") 
