# scalers15t

See `doc_diagram.pdf` for a flowchart of data and scripts. Rectangles are scripts to
be executed and parallelograms are files. You need to have done the polarisation
analysis using `polar15t` as well. Note that this `README` file serves as a more
descriptive guide as to what each script does; see the section "Scripts" below
for descriptions of the scripts (somewhat in order of execution). The section 
"Files" documents what is contained in each file

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
  if you change anything in `run_qa`, make sure you move all files in `datfiles/badruns` back
  to `datfiles/`

- You now have a set of `datfiles` and you are ready to update the rellum webpage; 
  see the "Updating Relative Luminosity Webpage" page of `doc_diagram.pdf`

  - removing bad bXings and update `counts.root` by doing the following:

    - execute `remake_bXing_dists` to remake all bXing distributions
    - scan through `pdf_bXings_fills_{lin,log,lin_zoom}/*.pdf` and look for any anomalies; 
      if you find any, add them to both `pathologies.txt` and `bunch_kicker` 
      (in `bunch_kicker`, they are lines like `echo [fill] [bXing] >> kicked_manual`
    - if you added any new bXings to `bunch_kicker`, you will now need to re-execute
      `bunch_kicker` and `accumulate` to produce a new `counts.root` file; if no new bad
      bXings were found in the update, you don't need to do this as `remake_bXing_dists`
      will have already produce a new `counts.root` file for you

  - Now that you have an updated `counts.root` file, you may proceed with the full rellum
    analysis; you can do this automatically by running `update_htmlfiles`; all the relevant 
    scripts will execute and a new `htmlfiles.tar.gz` will be produced, ready to be uploaded
    to whereve you wish


## Scripts
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
  - marks bad bXings; if you manually mark any as bad, you must have this script echo them
    to `kicked_manual` during execution
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
