# Purpose:

# Local c++ code:
- Read output `*.root` files from JETSCAPE. These contain:
  - four-momenta of final state particles (principally hadrons)
  - four-momenta of initiating parton (IP) from pp collision
     - only keep IP up to |\eta|<1.0
- Cluster final state particles into (``truth'') jets:
  - use anti-kT algorithm
  - jet resolution parameter = 0.4
  - only cluster partiles up to |\eta|<1.0
- Take the highest-pT (``leading'') IP and search within $|\Delta\R|<0.3$ for clustered jets
- Use the highest-pT jet within search
- Save the pT of the IP and the matched jet
- Generate a thermal background of particles up to $|\eta|<1.1$
- Embed the matched-jet constituents into the particle background and re-cluster (``reco'') jets
- Match the truth-jet to the highest-pT reco-jet within $|\DeltaR|<=0.3$ 
- Store the IP, truth-jet, and reco-jet to a TTree in a TFile
- Save rho median background density estimator from all the thermal particles (``rho_bkg_thermal'')
  - use the kT clustering algorithm
  - user R=0.4 resolution parameter
- Save the rho median background density estimator after embedding the JetScape final state particles:
  - same as for rho_bkg_thermal, but excluding the two highest-pT jets
