# Sural snRNAseq and Xenium analysis
This is the code of the sural snRNAseq and Xenium analysis of the corresponding manuscript. Link:
All scripts can be found in the scripts folder.

## Reproducibility

### renv
To ensure reproducibility, we used the *renv* package. To restore the environment from the *renv.lock* file, use:

```R
renv::restore()
```

## Docker
If you have trouble restoring the environment via *renv*, or if you want to be on the safe side regarding reproducibility,
you can also use the *Docker* image.

```bash
docker pull mihem/pns_atlas:v2.0
```

### Figures
To reproduce the figures, we created a self-contained *Quarto* document: [Reproducing the figures](https://mihem.github.io/pns_atlas/).
All required data are automatically downloaded from *Zenodo*  in the *Quarto* document.

**Steps to reproduce all figures:**
1. Download `docs/index.qmd`.
2. Pull the *Docker* image or restore the environment via *renv* as explained above.
3. Run the following command in the terminal:

```bash
quarto render index.qmd
```

## Questions?
f you have any questions, please contact us at [mheming.de](https://www.mheming.de/).