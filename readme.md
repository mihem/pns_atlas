# Sural snRNAseq and Xenium analysis
This is the code of the sural snRNAseq and Xenium analysis of the corresponding manuscript. Link:
All scripts can be found in the scripts folder.

## Reproducibility
To ensure reproducibility, we used the *renv* package. To restore the environment from the *renv.lock* file, use:

```R
renv::restore()
```

If you have trouble with restoring the environment via *renv*, or want to also use the same operating system,
you can also use the *Docker* image.

```bash
docker pull mihem/pns_atlas:v1.0
```

To reproduce the figures, we created a self-contained quarto report: [Reproducing the figures](https://mihem.github.io/pns_atlas/)

## Questions?
If you have any questions, please contact us via [mheming.de](https://www.mheming.de/).