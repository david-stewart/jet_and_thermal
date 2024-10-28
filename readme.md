# Analysis Code:

## JETSCAPE

The input data was generated with [site](https://github.com/jetscape).
Download getscape and compile and run it according to JETSCAPE
documentation. 

JETSCAPE is run with input `xml` files. The input files for this analysis are
located in the `jetscape_input` directory. All these input `.xml` files are
located in the `jetscape_input` directory. These files are in the following 
directories:

  1. `jetscape_input/pp/xml/`: The input files for the non-quenched jets
  2. `jetscape_input/lbt_brick/`: input files for jets quenched in `bricks` of QGP
  3. `jetscape_input/hydro/`: input files for the full hydro: Au+Au events,
      0-5% centrality with hydrodynamically simulated QGP evolution.

Note that each input `xml` file has a set values for the corresponding output
file (`outputFilename`) that is set from where these files were intially run
and will need to be updated appropriately.

## Process JETSCAPE output files:

Unzip the `.dat.gz` output files from JETSCAPE for the non-hydro runs and write
the output into input `*.root` files by:

  `./scripts/run_gunziptree.py [input-dir-path]`

This will use the the `scripts/gunziptree.py` on each of the files in the
provided file path.

For the hydro events, run `./scripts/gunziptree_hydro.py [input-file-name]` for
each hydro file individually. Additionally, use the ROOT `hadd` utility to
combine all the hydro `.root` files together into a single input file that contians
all the hydro backgrounds that will be used as the background embedding for the 
`pp` and `lbt_brick` jets.

## Run the C++ code to cluster the jets together and do the jet matching. For
the hydro events, the backgrounds in the hydro file will be used. Otherwise, an
input file from the `hydro` events will be used.

To do this:
`cd jet_and_thermal`
`make`
`./bin/IP_match_trees [opts]`

The options used depend on the run, and are:
  1. input file name
  2. maximum number of events, (-1 if using all events)
  3. output file name
  4. string: either 'is_hydro' if this is a hydro event, otherwise it is the
     name of the file containing the background distributions from the hydro
     runs which was previoudly `hadd`-ed together
  5. 0 or 1, determining if a vector of the pT's of all the jets constitents
     should be written to the output tree. Generally this isn't needed
  6. The number of background samples used (generally can be 1)
  7. The seed used for the random generator which isn't used in the default use
     of the code

If runing a single file, it is advisable to use a shell script to run these. An
example file is provided as
`jet_and_thermal/test_IP.sh`

For running all the events in a directory, can use script
`scripts/run_IP_match_trees.py`. Note that this script must be modified in order
to point to the hydro background file location, as well as to determin how many
core can be used simultaneously.

for the `pp` and `lbt_brick` runs, use ROOT's hadd utility to join the output 
files together. Then use the `script/to_parquet.py` script to convert the ROOT
files into `parquet` files, which is the file format used in the ipynb in the 
following steps.

##

General process:

The rest of this analysis was run in jupyter notebooks.




 1. Run JETSCAPE to simulate events. These are of 2 (/3) flavors:
    - Au+Au events, 0-5% centrality with hydrodynamically simulated QGP
    evolution.
    These are the most sophisticated simulations available to us.
    - Simulations in which jets are quenched through fixed length ``bricks'' of
    QGP.
    - Note that events with brick length zero are comparable to no-QGP.

    Both cases are run using standard JETSCAPE downloaded from the JETSCAPE 
    [site](https://github.com/jetscape). While running, the input files are
    `*.xml` files, and the output files are `*.dat.gz` files.

  2. The output `*.dat.gz` files are processessed into `*.root` files using
  python scripts, saved here as `gunziptree_hydro.py` and `gunziptree.py`, with
  a bundling script `run_gunziptree.py`.

  These scripts separate these following outputs in the JETSCAPE `.dat.gz`
  output files:
    1. The `# Energy loss Shower Initating Parton: JetEnergyLoss` section --
    which gives the energies of the initial hard parton scattering.
    2. The `# Final State Hadrons` section: give all the final state hadrons
    *from the initiaing parton shower* (i.e. these are the truth level jets)
    3. (in the case of the hydro events) `# Final State Bulk Hadrons`: the
    background particles from the QGP.
    4. the cross section for each event ("Xsec")

  The python scripts store these files into `.root` files. 

  The trees from the hydro files are collected separately so that they can be
  read in, separately, so that their backgrounds can be used with the jets
  quenched in the QGP events using bricks of plasma.

  3. The C++ code `IP_match_trees` is run using the `*.root` files from the
  previous step. This enacts the following logic:

    - Check if there is a separate input `.root` files for the bulk background.
    If there is, then read the background particles from that file for each 
    event. If not (this is a hydro event) then read the background from 
    each event.

    - Per-event, read in:
       - highest-pT IP: (phi, eta, pT)
       - the collisions Xsec
       - a vector of final state hadrons from the IP (phi, eta, pt)
       - a vector of the background hadrons (also phi, eta, pt)

       - If the |IP_eta| > 1., then discard the event (warning to me: this may
       be a source of error in the Xsec weightings...)
       
       - Cluster all the constitents from the IP, and check for all jets within
       dR<=0.4 of it. If there are none, just return that this IP wasn't
       matched.  If there are, then use the highest pT matched jet as the "Truth
       Jet"

       - cluster the background particles with the kT algorithm to get the 
       background median estimation. This is the `rho_bkg_thermal` in the trees.
       (as opposed to the additionally included `rho_bkg_thermalandjet`, but
       unused later in the analayis, clusteres the constituents of the Truth Jet
       together with the background particles, removes the high two pT jets, and
       then takes the median pT/area of the jets)

       - cluster together the background particles with the constituents the
       Truth Jet into "Reco Jets". Look at all Reco Jets within R<=0.3. If there
       are non, then write to the output tree that there is not a match from
       Truth Jet to Reco Jet in the `matched_TtoR` branch, and proceed to the
       next branch. If there is a match, then consider the highest-pT of
       the matched jet as the "Reco Jet"

       - Write the output parameters of the truth jet, reco jet (including 
       angularity and number of constituents, and the pT of the ten
       highest constituents) to the output file
  
  4. Use the `ROOT` commandline utility to concatenate the output files together:
  `hadd hadd.root *.root`

  5. Convert `hadd.root` to `hadd.parquet`. (Done in this case with 
  `./scripts/to_parquet` sometimes in conjunction with
  `./scripts/run_to_parquet.py`)

  6. Run the analysis on the input `.parquet` files using IPython notebooks
  under the scripts in `./Ipynb`











The `*.dat.gz` files divide the output into three 



# Purpose:

# Local c++ code:
- Read output `*.root` files from JETSCAPE. These contain:
  - four-momenta of final state particles (principally hadrons)
  - four-momenta of initiating parton (IP) from pp collision
     - only keep IP up to |\eta|<1.0
- Cluster final state particles into (``truth'') jets:
  - use anti-kT algorithm
  - jet resolution parameter = 0.4
  - only cluster partiles up to $|\eta|<1.0$ and $p_\mathrm{T}>0.2$
- Take the highest-pT (``leading'') IP and search within $|\Delta{}R|<0.3$ for clustered jets
- Use the highest-pT jet within search
- Save the pT of the IP and the matched jet
- Generate a thermal background of particles up to $|\eta|<1.1$
- Embed the matched-jet constituents into the particle background and re-cluster (``reco'') jets
- Match the truth-jet to the highest-pT reco-jet within $|\Delta{}R|<=0.3$ 
- Store the IP, truth-jet, and reco-jet to a TTree in a TFile
- Save rho median background density estimator from all the thermal particles (``rho_bkg_thermal'')
  - use the kT clustering algorithm
  - user R=0.4 resolution parameter
- Save the rho median background density estimator after embedding the JetScape final state particles:
  - same as for rho_bkg_thermal, but excluding the two highest-pT jets
