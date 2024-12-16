# Analysis Code and Process:

## Run JETSCAPE

The input data was generated with [site](https://github.com/jetscape).
Download JETSCAPE, then compile and run it according to JETSCAPE
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

## Run the C++ Code to Cluster and Match Jets

Use the `.root` files from the step above as input for the C++ binary. For the hydro
events, the backgrounds in the hydro file will be used. Otherwise, an input
file from the `hydro` events will be used.

To do this:

`cd jet_and_thermal`
`make`
`./bin/IP_match_trees [options]`

The options used depend on the run, and are:
  1. input file name
  2. maximum number of events, (-1 if using all events)
  3. output file name
  4. string: either 'is_hydro' if this is a hydro event, otherwise it is the
     name of the file containing the background distributions from the hydro
     runs which was previoudly `hadd`-ed together
  5. 0 or 1, determining if a vector of the pT's of all the jets constitents
     should be written to the output tree. Generally this isn't required.
  6. The number of background samples used (generally can be 1)
  7. The seed used for the random number generator (note that this isn't used in 
          the code's default use case, or in the analysis as described here)

## Transform Output data in Input `.parquet` Files for Ipynb

The output files from `IP_match_tree` are all root files. For all the `pp` and
`lbt_brick` runs, use ROOT's hadd utility to join the output files together. 
Convert these new output files (the joined ones for `pp` and `lbt_brick` runs and 
individual block files one-by-one) to `.parquet` files. To do this, 
`script/to_parquet.py` may be used.

## Do the analyses with the `Ipynb` files

### Generate input files flat in pT truth jet:

In `./Ipynb/NoBrickInput`:

Here `NoBrickInput` is just misnomer for `pp` (no quenching) input.

 - Run `min_bias_input.ipynb`. Make sure that notebook has the input path to
   the parquet file containing all the `pp` jets from the `IP_match_trees` C++
   process.
   
   output: `min_bias_input.parquet`

   This is the output file used for training Neural Networks.

 - Run `build_Xsec.ipybn`
   input: the input files from `IP_match_trees` (same as for `min_bias_input.ipybn`)
   output: directory `Xsec_groups` with dataframes for each Xsec bin
  
 - Run `build_Xsec_clean.ipynb`
   input: output of `build_Xsec.ipynb`
   output: files in which outliers have been removed
   use of output: files required to make a weighted spectrum for the jets

### Train the Neural Network:

In `./NeuralNetwork`:

 - `mkdir train`

  Run the following notebooks:
  n.b.: make sure that the input path to the flattened input is present.
  - `train_reco_cpt.ipybn`
  - `train_reco_all.ipybn`
  - `train_angularity.ipybn`
  - `train_nconsts.ipybn`
  - `train_rhoA.ipybn`

  Run the predictions on all of the input. Note that these use a local
  python file `JetscapeFileGetter.py`. Note that this script may need to 
  be updated based on your file paths. This is done in each of 
  these notebooks:

  - Run the 5 NNs on the `pp` and `lbt_brick` runs:
    - Predict_angularity.ipynb
    - Predict_nconsts.ipynb
    - Predict_reco_all.ipynb
    - Predict_reco_cpt.ipynb
    - Predict_rhoAonly.ipynb

  - Run the NNs on the hydro runs:
    - Predict_reco_hydro_all.ipynb
    - Predict_reco_hydro_angularity.ipynb
    - Predict_reco_hydro_cpt.ipynb
    - Predict_reco_hydro_nconsts.ipynb
    - Predict_reco_hydro_rhoAonly.ipynb

  For each input file, an output `.json` file containing a few summary
  statistics are generated in addition to a few plots.

  Collect and plot a summary the results by running
  - plot_allscores_paper.ipynb. This notebook uses the python file
          `../JetMatchesSet.py`. 

### Analysize the Fragmentation Function:

  - Make the paper figure from the fragementation of the input jets with
    `./FragFunc/FIG_fragfn.ipynb`.
    This notebook uses the output files from the C++ program `IP_match_trees`, as
    it reads the pT of all the constituents of each jet.

### Make a `pp` and a quenched spectra:

  In `Ipynb/Threep5FermiInput`:

  - Run `min_bias_input.ipnb` (make sure the input file is to the `.parquet` file
  with all the 3.5 fm input data)

 - Run `build_Xsec.ipybn`
 - Run `build_Xsec_clean.ipynb`

 These last two files are entirely analogous to those in `NoBrickInput`.
 Manually check the directories so that there are comparable events in each
 cross section bin for both. In analysis, two of the three highest cross-section
 groups were represented in the NoBrickInput data (although this is all essentially
 min-bias data). If needed, remove them here, too.

  - Apply the neural network correction on the bricks, so run:
    - M_angularity_NoBrick/ApplyModel.ipynb
    - M_AB_method_NoBrick/ApplyModel.ipynb note: this isn't a nueral network
    - M_reco_all_NoBrick/ApplyModel.ipynb
    - M_reco_cpt_NoBrick/ApplyModel.ipynb
    - M_rhoAonly_NoBrick/ApplyModel.ipynb
    - M_nconsts_NoBrick/ApplyModel.ipynb

### Compare RAA like results

  In `RAA_like`:
  Make sure the input files are correct and run:
   -  `./UnfoldPapar.ipynb`
   -  `./plot_RAA_rat.ipynb`
   -  `./UnfoldPapar.ipynb`

   Note that the input looks must be modified to fun on all the desired inputs.
   The logic will collect and weight the spectra from each of the input files
   and add them together. It will also make a response matrix and do the unfolding.
   

#  Appendix: Some Additiona and Repeated Logic and Detail

## Jetscape:

   The JETSCAPE files are 2 (/3) flavors:

    - Au+Au events, 0-5% centrality with hydrodynamically simulated QGP
    evolution.
    These are the most sophisticated simulations available to us.
    - Simulations in which jets are quenched through fixed length ``bricks'' of
    QGP.
    - Note that events with brick length zero are comparable to no-QGP.

    Both cases are run using standard JETSCAPE downloaded from the JETSCAPE 
    [site](https://github.com/jetscape). While running, the input files are
    `*.xml` files, and the output files are `*.dat.gz` files.

## Jetscape output Processing Files

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

## C++ Code:

 `IP_match_trees` is run using the `.root` files from the
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
