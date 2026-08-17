# Savont recall by observed pseudo-reference depth

For each of the 24 actual ONT subsamples, the workflow:

1. dereplicates reads with USEARCH `-fastx_uniques -minuniquesize 1 -sizeout`;
2. assigns those reads directly to the matching HiFi-derived UNOISE3
   pseudo-reference with USEARCH `-otutab`;
3. uses the resulting per-reference count as that HiFi ASV's observed depth in
   that particular ONT subsample; and
4. marks the HiFi ASV as recovered when a Savont ASV has an `NM:i:0`
   alignment to it (terminal clipping is permitted).

Each point in the figure is the fraction of eligible HiFi ASVs recovered in
one dataset and depth bin. Bar heights are unweighted means of those dataset
fractions. HiFi ASVs assigned fewer than 10 ONT reads are outside the plotted
bins.

Current outputs use `observed_depth` and `direct_mapping` in their names. The
`legacy_ont_unoise_bridge/` contains superseded tables from an outdated method.
