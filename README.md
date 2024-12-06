# MixedUpParents-sims-and-analyses

## Notes on SEDNA

Something weird goes on with SEDNA whereby the slim_sim rule writes all its output
no problem, and then fails when returning control back to Snakemake.  It might
have something to do with the reticulate package and its conda environment.

At any rate, you can run it like this:
```sh
snakemake -np --use-envmodules --keep-incomplete --until slim_sim --profile hpcc-profiles/slurm/sedna
```
and once that is done, the outputs will still be there and you can continue with:
```sh
snakemake -np --use-envmodules --use-conda --ignore-incomplete --profile hpcc-profiles/slurm/sedna
```

Total PITA.  Gonna have to figure out why it is failing, but we can use the above
until that time.

Note, the --use-conda on the second line is to have bedtools for mup_logls.

You gotta define a MUP_CONDA variable, too. see the readme for MixedUpSlimSims for that.


## To-Dos:

1. Write a separate ROC rule for sequoia using the tibble that Jared has already made.
2. In sequoia.R figure out how to use only the variable markers or both sets.



