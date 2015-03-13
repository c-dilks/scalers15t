# scalers15t
run 15t relative luminosity


## Adding new runs to analysis

- On RCAS, follow directions in `~/scalers2015/README.md` to download new files from HPSS

- Download datfiles to local machine: `get_scaler_files` 

- Update run list `culled_run_list.txt` by running `update_runlist`

- Update spin pattern directory by running `getspinpats`

- Update `spinpat/*.spin` files by running `spin_cogging`

- Execute `read_scalers` to produce datfiles which haven't been produced yet
