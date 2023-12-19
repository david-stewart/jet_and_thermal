# Update the code to mimic the input variables in arXiv:2303.00592v1 

 -- Based on these considerations, the following input features were selected: 
 * jet transverse momentum corrected by the standard area-based method
 * first radial moment of constituent momenta (jet angularity) 
    $\equiv g = \Sigma_{i\injet}\frac{p_{T,i const}}{p_{T,jet}} \DeltaR_{jet to i const}
 * number of constituents within the jet
 * transverse momenta of the eight leading (highest pT) particles within the jet
 
 -- Target is PYTHIA pt



# Purpose

This code embed Pythia8 jets into a random background isotropically distributed
in eta<|1|, phi\in[-pi,pi]

# How to build

  * Clone this code locally. Generate directories:
      `./bin`
      `./obj`

  * Make a virtual link to an installation of FastJet3
  * Update `./Makefile` for paths for libraries for `root6`, `Pythia8`, and `fastjet3`
  * run: `./make`

## updated to do, 11.08.2023
  - [ ] Add ability to do background subtration (i.e. calculate rhoxA)
  - [ ] Add ability to look into the clustering history
  
## 11.13.2023
  - pythia8 -- see [PYTHIA8 documentation](https://pythia.org/latest-manual/CrossSectionsAndWeights.html)
    - store event info ::weightSum() -- what about sigmaGen()?

