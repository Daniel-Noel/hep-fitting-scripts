## HEP fitting scripts

Scripts to run fits with pyhf and produce outputs in cabinetry


To set up the pyhf and cabinetry environments
```
. setup.sh
```


For producing yield tables where we fit in a control region then extrapolate to a region not fit in
```
python yield_tables.py
```

For producing post-fit plots where we fit in a control region and apply fit to variables
```
python postfit_plots.py
```


For producing a systematic table for all systematics, individual systs or combined systs
```
python syst_table.py
```
