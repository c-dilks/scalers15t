# scalers15t

See `doc_diagram.pdf` for a flowchart of data and scripts. Rectangles are scripts to
be executed and parallelograms are files. You need to have done the polarisation
analysis using `polar15t` as well. Note that this `README` file serves as a more
descriptive guide as to what each script does; see the section "Scripts and Files" below
for descriptions of the scripts and data files (somewhat in order of execution).

It is best to read through `doc_diagram.pdf` and search this document for details at each
step; most of the flowcharts are read from top-left corner to bottom-left corner

Caveat: if this is your first time running and you have cloned from git, you may need to
make some subdirectories such as `datfiles`; scan through `doc_diagram.pdf`


## Adding new runs to analysis

- On RCAS, follow directions in `~/scalers2015/README.md` to download new files from HPSS

- Download datfiles to local machine: `get_scaler_files` 

- Update run list `culled_run_list.txt` by running `update_runlist`

- Update spin pattern directory by running `getspinpats`

- Update `spinpat/*.spin` files by running `spin_cogging`

- Execute `read_scalers` to produce datfiles which haven't been produced yet

- Execute `run_qa` and the runs which do not pass the QA will be moved to `datfiles/badruns`;
  this script only QAs runs which have run number greater than that in `lastrun_qa`
  - if you change anything in `run_qa`, make sure you move all files in `datfiles/badruns` back
    to `datfiles/`

- You now have a set of `datfiles` and you are ready to update the rellum webpage; 
  see the "Updating Relative Luminosity Webpage" page of `doc_diagram.pdf`

  - first remove bad bXings and update `counts.root` by doing the following:

    - execute `remake_bXing_dists` to remake all bXing distributions
    - scan through `pdf_bXings_fills_{lin,log,lin_zoom}/*.pdf` and look for any anomalies; 
      if you find any, add them to both `pathologies.txt`  
    - if you added any new bXings to `patholigies.txt`, you will now need to re-execute
      `bunch_kicker` and `accumulate` to produce a new `counts.root` file; if no new bad
      bXings were found in the update, you don't need to do this as `remake_bXing_dists`
      will have already produce a new `counts.root` file for you

  - Now that you have an updated `counts.root` file, you may proceed with the full rellum
    analysis; you can do this automatically by running `update_htmlfiles`; all the relevant 
    scripts will execute and a new `htmlfiles.tar.gz` will be produced, ready to be uploaded
    to whereve you wish


## Scripts and Files
- `get_scaler_files`
  - uses `scp` to get files from `RCAS`

- `read_scalers`
  - makes a condor job for executing the scaler reader for all scaler board files in the
    subdirectory `sca2015`; the files are output to `datfiles/*.dat` and contain columns
    [bXing] [bbc bits 0-7] [zdc bits 0-7] [vpd bits 0-7]
  - the scaler bit reader is called `scaler2_reader_bit.c`; you need to compile it
    using `make` which produces the binary `scaler2_reader_bit.exe`
    - this file is the 32-bit reader for run15; see `bit_doc.txt` and `bit_map.pdf` 
      as well as the scaler router map on the trigger webpage for further info about
      the scaler bit maps

- `run_qa`
  - for all logical scaler bits (x0 w0 e1; x0 w1 e0; x1, w1, e1), we sum the scales 
    over all bXings
  - if all 3 logical scales for all 3 detectors each are greater than `$THRESHOLD`, then
    the run is ok; if not, it is moved to `datfiles/badruns`

- `update_runlist`
  - gets most recent `culled_run_list.txt` from `RCAS`

- `getspinpats`
  - gets all spin pattern files from `CDEV` for runs listed in `culled_run_list.txt` and 
    places them in `spinpat`
  - only spin patterns which you do not yet have are downloaded
  - `fill.txt` is regenerated which is a list of run numbers and fill numbers
  - fills for which no spin pattern was found are printed out at the end of execution and
    culled from `fills.txt`

- `spin_cogging`
  - creates `spinpat/*.spin` files which contain columns of STAR bXing, blue spin, yell spin; 
  - yellow beam is cogged so that blue bunch 0 collides with yellow bunch 80
  - since STAR spin is opposite that of CDEV spin, we implement a sign filp here

- `bunch_kicker` 
  - marks bad bXings; if you manually mark any as bad, you must do so in `pathologies.txt`
  - if `run_randomizer=1` then it kicks the minimum amount of additional bXings such that
    the number of spin bits for the fill are equalized (see comments in `bunch_kicker` 
    for more details)
  - output file is a list of kicked bXings is `kicked` with columns `fill, bXing, spinbit`

- `accumulate`
   - collects all the datfiles into `datfiles/acc.dat`, with columns:
     (** and filters out bad part of 17600)
     - run index
     - runnumber
     - fill
     - run time (seconds)
     - bunch crossing (bx)
     - BBC[1-8]
     - ZDC[1-8]
     - VPD[1-4]
     - total bXings
     - blue spin
     - yellow spin
 - `accumulate` then creates `counts.root,` which contains useful trees, using `mk_tree.C`
 - `mk_tree.C` reads `datfiles/acc.dat` file
   - BXING SHIFT CORRECTIONS ARE IMPLEMENTED HERE (FOR RUN 12 ONLY!!!)
   - acc tree: simply the `acc.dat` table converted into a tree
   - sca tree: restructured tree containing containing branches like 
     - bbc east, bbc west, bbc coincidence 
     - similar entries for zdc and vpd
     - run number, fill number, bunch crossing, spin bit
     - num_runs = number of runs in a fill
     - kicked bunches: 
       certain bunches which are empty according to scalers but filled according to
       cdev are labelled as 'kicked' in the output tree; bXings which are kicked
       are not projected to any distributions in the rellum4 output

- `draw_spin_patterns.C`
  - produces a pdf of spin patterns; vertical axes are bXing numbers and horizontal
    axis is the spin
  - there are 2 fills per page, one plot for blue beam and one for yellow beam for 
    each fill
  - red means opposite spin bXing, grean means same spin, cyan means no collisions

- `bit_combinations.C`
  - produces images in `bit_combos` which plot the number of each of the scaler `3bc's`

- `nbx_check_2.C`
  - produces plots of the number of bXings vs. bXing number for each run, normalized so that
    the max bin is equal to unity... this is an odd structure which no one fully understands

- `rellum4.C`
  - this is the analysis script
  - `var` is the independent variable
    - `i` is run index, `fi` is fill index, `bx` is bXing
  - objects in outputted `rdat{i,fi,bx}.root` file:
    - `c_spin_pat` -- `R_spin = same spin / diff spin`
    - `c_raw_{bbc,zdc,vpd}` = raw scaler counts vs var
    - `c_acc_{bbc,zdc,vpd}` = accidentals corrected scaler counts vs var
    - `c_mul_{bbc,zdc,vpd}` = multiples corrected scaler counts vs var
    - `c_fac_{bbc,zdc,vpd}` = correction factor (mult/raw) vs var
    - `c_R#_{bbc,zdc,vpd}` = relative luminosity vs. var
    - `c_mean_R#` = mean rellum over EWX (bbc,zdc,vpd on same canvas)
    - `c_R#_zdc_minus_vpd` = difference between zdc and vpd
    - `c_deviation_R#_{bbc,zdc,vpd}` = rellum minus mean rellum
    - `rate_dep_R#_{bbc*,zdc*,vpd*}` = rellum vs. multiples corrected rate
    - `c_rate_fac_{bbc*,zdc*,vpd*}` = correction factor vs. multiples corrected rate
                                     (tprofile from `rate_fac` for each spinbit)
  - there are more notes in the header comments of this script

- `rellum_all`
  - basically used to run `rellum4.C` for various independent
    variables etc. 
  - outputs pngs in `png_rellum`, ready to be copied to protected area to link
    to scalers web page

- `rellum_fills [drawLog] [zoomIn]`
  - runs `rellum4.C` for all fills separately and output pdfs in subdirectories of
    `pdf_bXings_fills`; this is for looking at fill dependence of bXing 
    distributions
  - execute `ghost_script_fills` afterward to combine all the pdfs into 
    `pdf_bXings_fills/*.pdf`
  - produces `matrix` trees (see matrix section below)

- `rellum_runs`
  - analagous to `rellum_fills` but for run-by-run bXing distributions
  - use `ghost_script_runs` for pdf concatenation
  - also produces `matrix` trees (see matrix section below)

- `sumTree.C`
  - builds `sums.root` which sums the counts in `counts.root` for each run
  - also recognizes the spin pattern which is used; see the script comments
    for the variable `subpattern` for further details
  - details about the pattern recognition are written out to `pattern_log.txt`

- `nbx_check.C`
  - compares total number of bunch crossings divided by clock frequency to 
    run time, which should be approximately equal
  - plots are output to the `nbx_check/` subdirectory

- `combineAll.C`
  - combines the `sums.root` file with the `rdat_i.root` file to produce a final
    relative luminosity data tree contained in `rtree.root`

- `make_run_table.C`
  - builds `run_table.txt` which is a text table of run numbers, run indices, 
    fill numbers, and fill indices, to be copied to the webpage

- `remake_bXing_dists`
  - clears `kicked` and reruns `accumulate` to produce a `counts.root` file
    with no kicked bXings
  - then it runs `rellum_fills` for three different drawing settings, merges 
    the pdfs with `ghost_script_fills` and moves them to three separate output 
    directories: `pdf_bXings_fills_{lin,lin_zoom,log}`
  - then it reruns `bunch_kicker` and `accumulate` to recreate `counts.root`
    with the original kicked bXings list in the tree
  - you may then go through these pdfs and hunt for bad bXings; if you find any
    new ones, add them to `bunch_kicker` and re-execute `bunch_kicker` and 
    `accumulate`

- `update_htmlfiles`
  - reruns all of the other scripts (besides `remake_bXing_dists`) to build all
    the relevant files for the webpage; it then executs `copyToHtmlfiles` to copy
    all relevant files into the `htmlfiles` subdirectory
  - a tarball of `htmlfiles` is then created


## Matrix Subdirectory: Bunch Fitting Algorithm

- Running `rellum4.C` with `var="bx"` and with `specificFill>0` XOR `specificRun>0` will
  produce `matx` tree files, found in `matrix/rootfiles/*.root`
  - this can be done for each fill or run using `rellum_fills` or `rellum_runs`
  - the `matx` tree contains scales, corrected scales, and correction factors for each 
    cbit, tbit, and bXing
- execute `hadd matrix/rootfiles/all.root matrix/rootfiles/matx*.root` to merge the matx trees
- `DrawMatrix.C` draws the desired matrix and produces `matrix.root` and `for_mathematica`
  - `matrix.root` contains the matx tree and the matrix `mat`
  - `for_mathematica` contains the matrix `mat` in text form for reading with mathematica
  - singular value decomposition (SVD) is then performed using `SVD.nb`
- bunch fitting
  - the main code for bunch fitting is `BF.C`, which requires `rootfiles/all.root`, 
    `../counts.root`, and `../sums.root`
  - you need to specify the ratio of scalers to bunch fit over
    - numerator tbit and cbit (see `rellum4.C` for definitions of tbit and cbit)
    - denominator tbit and cbit
    - evaluateChiSquare = true will try to draw Chi2 profiles (not working well yet...)
    - specificPattern != 0 will only consider the specified spin pattern
  - it's best to just use `RunPatterns`, which runs `BF.C` under various interesting
    conditions, for all spin patterns; the following files are produced:
    - `fit_result.[num].[den].root`: bunch fit results for all spin patterns, where the fit
      is done to `r^i=num/den`
    - `pats/fit_result.[num].[den].pat[pat].root`: bunch fit results for spin pattern `pat`
    - `colour.[num].[den].root`: bunch fit results, with colour code according to spin 
      patterns (see the TCanvas `legend` in the ROOT file)
