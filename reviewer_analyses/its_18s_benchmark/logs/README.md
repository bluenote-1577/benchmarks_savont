# Execution records from the manuscript run

These logs were retained from the analysis runs and cover all six nested
depths:

- `savont/`: 12 complete ASV calls (six each for 18S and ITS), including the
  command parameters, stage counts, timing, and successful-completion marker;
- `unoise/`: 72 USEARCH dereplication, denoising, and read-to-zOTU mapping
  records for default UNOISE3 and `-minsize 3`;
- `banana/`: six complete, losslessly gzip-compressed BaNaNA command/stage logs;
- `classification/`: 42 harmonized Savont classification records for every
  marker, method, and depth;
- `database/`: the two minimap2 index-construction logs; and
- `input_preparation/`: a compact summary of the original 298 MB nhmmer report,
  including its HMMER header, model, threshold, thread count, and extraction
  outcome.

Raw FASTQs, full nhmmer hit tables/reports, databases, PAFs, and full run trees
are omitted because they are too large for Git. See `../HOWTO.md` to reconstruct
the inputs and rerun the workflow.
