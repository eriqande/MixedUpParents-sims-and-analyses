# MixedUpParents-sims-and-analyses

## Notes on SEDNA

Something weird goes on with SEDNA whereby the slim_sim rule writes all its output
no problem, and then fails when returning control back to Snakemake.  It might
have something to do with the reticulate package and its conda environment.

At any rate, you can run it like this (which is assuming that we are running snakemake
on a node with 10 cores, 9 of which we will use for local jobs):
```sh
snakemake -np --use-envmodules --keep-incomplete --until slim_sim --profile hpcc-profiles/slurm/sedna
```
and once that is done, the outputs will still be there and you can continue with:
```sh
snakemake -np --use-envmodules --use-conda --ignore-incomplete --local-cores 9 --profile hpcc-profiles/slurm/sedna
```

Total PITA.  Gonna have to figure out why it is failing, but we can use the above
until that time.

Note, the --use-conda on the second line is to have bedtools for mup_logls.

You gotta define a MUP_CONDA variable, too. see the readme for MixedUpSlimSims for that.




