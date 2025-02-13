# Readme

This is a very simple package with a main object, an `R6` class called `hapBlock` with several methods. Its main purpose is to analyse groups of phased SNP markers and to convert them to haplotype codes.

Although this task is relatively simple, when the number of SNP markers increases (especially after having more than 10 SNPs) the number of "rare" alleles is quickly inflated. This is because a single mistake in a single marker produces a seemingly unique, new, haplotype. An easy solution to this situation is to compare haplotypes and group them based on similarity. This is exactly what `hapgroup` does. Find out more in the `vignette`.

To install `hapgroup` you can use the following code:

```{r}
devtools::install_github("https://github.com/alethere/hapgroup",build_vignettes = TRUE)
```

And to see the vignette you can go to this link or browse it from the package:

```{r}

```
