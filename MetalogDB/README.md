# Fetching microbiome profiles from Metalog as a `TreeSummarizedExperiment`

[Metalog](https://metalog.embl.de/) is a database of curated metagenomic taxonomic profiles and harmonised sample metadata, organised into
four collections: `"human"`, `"animal"`, `"ocean"` and `"environmental"`.

`fetchMetalogTSE()` downloads one or more of those collections and returns them
as a
[`TreeSummarizedExperiment`](https://bioconductor.org/packages/TreeSummarizedExperiment)
(TSE), the container used throughout the
[`mia`](https://bioconductor.org/packages/mia) microbiome-analysis ecosystem.
Profiles are MetaPhlAn4, at SGB (species-level genome bin) resolution.

This document walks through the function argument by argument. The outputs shown
were produced by actually running the code against the `2026-07-07` Metalog
snapshot, your own numbers will differ as the database is updated.

## Contents

- [Setup](#setup)
- [Quick start](#quick-start)
- [Anatomy of the returned object](#anatomy-of-the-returned-object)
- [`collection`: one or several collections](#collection-one-or-several-collections)
- [`meta.type`: how much sample metadata to fetch](#metatype-how-much-sample-metadata-to-fetch)
- [`samplelist`: restrict to a subset of samples](#samplelist-restrict-to-a-subset-of-samples)
- [`use.cache`: reusing downloads, and forcing a refresh](#usecache-reusing-downloads-and-forcing-a-refresh)
- [`make.dense`: dense vs sparse assays](#makedense-dense-vs-sparse-assays)
- [Provenance](#provenance)
- [A short downstream example](#a-short-downstream-example)
- [Caching and troubleshooting](#caching-and-troubleshooting)

## Setup

The function lives in a single, standalone-sourceable script. Install its
dependencies, then source it:

```r
# mia pulls in TreeSummarizedExperiment and S4Vectors through its Depends chain
# BiocManager::install("mia")
# install.packages(c("data.table", "ape", "Matrix", "httr2", "jsonlite"))

source("fetchMetalogTSE.R")
```

The script attaches `mia` for you and stops with an informative error if any of
the six packages is missing. An internet connection is required for the first
call; later calls can run entirely from the cache.

## Quick start

The simplest call downloads a single collection together with its core metadata:

```r
tse <- fetchMetalogTSE("ocean")
tse
#> class: TreeSummarizedExperiment 
#> dim: 7861 5086 
#> metadata(1): metalog
#> assays(1): relabundance
#> rownames(7861): t__SGB100018 t__SGB100073 ... t__SGB99616 t__SGB99643
#> rowData names(8): Kingdom Phylum ... Species SGB
#> colnames(5086): Alneberg_2020_marine_Baltic.P4201_101
#>   Alneberg_2020_marine_Baltic.P4201_102 ...
#>   Zhou_2022_deepsea_sediment.T3L8-5 Zhou_2022_deepsea_sediment.T3L8-7
#> colData names(28): artificial collection_date ... temperature
#>   collection
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowLinks: a LinkDataFrame (7861 rows)
#> rowTree: 1 phylo tree(s) (7861 leaves)
#> colLinks: NULL
#> colTree: NULL
```

`"ocean"` is one of the smaller collections and makes a good first call.
`"human"` is by far the largest — its profile table alone is a 327 MB gzip!

The first call prints download progress; subsequent calls report what they took
from the cache instead:

```
#> Loaded cached file: /home/you/.cache/R/metalog/metaphlan4_clades.tsv.gz
#> Loaded cached file: /home/you/.cache/R/metalog/ocean_metaphlan4_2026-07-07.tsv.gz
#> Loaded cached file: /home/you/.cache/R/metalog/ocean_core_long_2026-06-07.tsv.gz
#> Loaded cached file: /home/you/.cache/R/metalog/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk
```

## Anatomy of the returned object

Everything is reachable through the standard `SummarizedExperiment` / `mia`
accessors.

```r
assayNames(tse)
#> [1] "relabundance"

assay(tse, "relabundance")[1:5, 1:3]
#>                   TARA_Oceans_protists.sample_TARA_N000001485
#> t__SGB50875                                           0.11359
#> t__SGB3824                                            0.45890
#> t__SGB50962                                           0.09359
#> t__SGB19048_group                                     1.36188
#> t__SGB68547                                           0.11506
```

The single assay holds MetaPhlAn4 relative abundances as percentages. Only
`t__SGB` rows are kept from the profile, so every feature is one SGB; the higher
ranks are recovered from the taxonomy rather than from separate rows.

```r
head(rowData(tse), 3)
#> DataFrame with 3 rows and 8 columns
#>                  Kingdom            Phylum                  Class
#>              <character>       <character>            <character>
#> t__SGB100018 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
#> t__SGB100073 k__Bacteria p__Proteobacteria c__Alphaproteobacteria
#> t__SGB10009  k__Bacteria p__Proteobacteria c__Gammaproteobacteria
#>                            Order                 Family               Genus
#>                      <character>            <character>         <character>
#> t__SGB100018     o__Moraxellales       f__Moraxellaceae    g__Acinetobacter
#> t__SGB100073 o__Hyphomicrobiales f__Methylobacteriaceae g__Methylobacterium
#> t__SGB10009  o__Enterobacterales      f__Morganellaceae      g__Providencia
#>                            Species          SGB
#>                        <character>  <character>
#> t__SGB100018 s__Acinetobacter_mod.. t__SGB100018
#> t__SGB100073 s__Methylobacterium_.. t__SGB100073
#> t__SGB10009  s__Providencia_heimb..  t__SGB10009
```

`rowData` carries eight ranks, Kingdom through SGB, mapped from the MetaPhlAn4
clade table. `colData` holds the sample metadata, pivoted from the long format
Metalog serves into one column per metadata item:

```r
colData(tse)[1:3, c("collection_date", "depth_meters", "temperature")]
#> DataFrame with 3 rows and 3 columns
#>                                       collection_date depth_meters temperature
#>                                           <character>  <character> <character>
#> Alneberg_2020_marine_Baltic.P4201_101      2013-01-23            2         3.1
#> Alneberg_2020_marine_Baltic.P4201_102      2013-04-03            2         2.2
#> Alneberg_2020_marine_Baltic.P4201_103      2013-04-16            2         3.1
```

Note that metadata values arrive as **character** columns regardless of their
content. Convert the ones you intend to compute with.

Finally, the MetaPhlAn4 SGB phylogeny tree is attached, pruned to the features
present:

```r
rowTree(tse)
#> Phylogenetic tree with 7861 tips and 7860 internal nodes.
#> 
#> Tip labels:
#>   139512, 139514, 90923, 19, 49495, 49494, ...
#> 
#> Rooted; includes branch length(s).
```

## `collection`: one or several collections

Pass a character vector to download and merge several collections into one
object. Features are unioned, samples concatenated, and a `collection` column
records each sample's source:

```r
tse <- fetchMetalogTSE(c("ocean", "animal"))

dim(tse)
#> [1] 19241 14846

table(tse$collection)
#> animal  ocean 
#>   9760   5086
```

Metadata columns are unioned too, and samples from a collection that does not
define a given column are filled with `NA`. That is easy to see in the merged
object — the oceanographic variables are populated only for ocean samples, the
host variables only for animal samples:

```r
cd <- as.data.frame(colData(tse))
sapply(split(cd, cd$collection), function(x) colSums(!is.na(x)))[
    c("depth_meters", "salinity", "host_tax_id", "sex"), ]
#>              animal ocean
#> depth_meters      0  5086
#> salinity          0  2190
#> host_tax_id    9760     0
#> sex            4604     0
```

Merging ocean and animal core metadata therefore widens `colData` from 28
columns to 46.

Sample aliases must be unique across collections. A collision is reported as an
error rather than silently merged:

```
#> Error: Duplicate sample aliases across collections: ...
```

## `meta.type`: how much sample metadata to fetch

One of `"core"` (the default), `"extended"`, `"all"` or `"none"`. The levels
differ only in how many metadata columns you get; the assay is unaffected.

| `meta.type`  | What it contains                                 | `colData` columns (ocean) |
| ------------ | ------------------------------------------------ | ------------------------- |
| `"core"`     | Curated fields available for most or all samples | 28                        |
| `"extended"` | The partially harmonised set                     | 40                        |
| `"all"`      | Everything, down to study-specific variables     | 286                       |
| `"none"`     | No metadata download at all                      | 0                         |

(Those counts are what `ncol(colData())` returns, so they include the added
`collection` column.)

```r
tse_core <- fetchMetalogTSE("ocean", meta.type = "core")
ncol(colData(tse_core))
#> [1] 28

tse_ext <- fetchMetalogTSE("ocean", meta.type = "extended")
ncol(colData(tse_ext))
#> [1] 40
```

`"all"` is a different kind of object to work with. Study-specific variables are
included verbatim, so column names can be very long and most columns are sparsely
populated:

```r
tse_all <- fetchMetalogTSE("ocean", meta.type = "all")
ncol(colData(tse_all))
#> [1] 286

head(colnames(colData(tse_all)), 3)
#> [1] "AOU"                                                                 
#> [2] "Atom_ratio_of_dissolved_C_isotopes_in_dissolved_inorganic_carbon_ex..
#> [3] "Atom_ratio_of_dissolved_Fe_isotopes_expressed_in_conventional_DELTA..
```

Use `"none"` to skip the metadata download entirely when you only need
abundances. For a single collection this leaves `colData` genuinely empty:

```r
tse_min <- fetchMetalogTSE("ocean", meta.type = "none")
ncol(colData(tse_min))
#> [1] 0
```

The `collection` column is added in every other case, including single-collection
calls. `meta.type = "none"` with exactly one collection is the only combination
that produces a zero-column `colData`.

## `samplelist`: restrict to a subset of samples

The [Metalog web UI](https://metalog.embl.de/) lets you explore samples
interactively and export your selection as a file. Pass that file to
`samplelist` and the subsetting is applied *while the object is built*, rather
than after the fact.

That ordering matters for what the object costs you. The profile is still
downloaded and parsed in full, so you need enough memory to hold the whole
collection sparsely — `samplelist` will not let you squeeze `"human"` onto a
small machine. But the subset is taken before the metadata pivot, the tree
pruning, the densification and any merging, so everything downstream of the
filter works on your slice alone.

The file may be `.csv`, `.tsv` or `.json`, and **must contain a `sample_alias`
column**:

```csv
"sample_alias"
"Alneberg_2020_marine_Baltic.P4201_124"
"Alneberg_2020_marine_Baltic.P6071_515"
"Biller_2018_GEOTRACES_marine_pelagic.S0004"
```

```r
tse_sub <- fetchMetalogTSE("ocean", samplelist = "my_samples.csv")
#> Filtering to requested samples in 'ocean'...

dim(tse_sub)
#> [1] 3558   50
```

Only the listed samples are kept, and features that are zero across all of them
are dropped as well — 50 ocean samples retain 3558 of the collection's 7861
features. Aliases in the list that are absent from the collection are ignored, so
a list that only partly matches still works.

Two things to watch for:

- **The list must actually name a `sample_alias` column.** A file with the
  aliases under a different name fails with the same error as a list that matches
  nothing: `None of the samples in 'samplelist' were found in the dataset.`

- **Filtering is applied per collection.** Each requested collection is filtered
  independently, and any collection left with no matching samples aborts the
  whole call. A sample list drawn only from `"ocean"` therefore *fails* when
  combined with `c("ocean", "animal")`, rather than returning just the ocean
  samples:
  
  ```r
  fetchMetalogTSE(c("ocean", "animal"), samplelist = "ocean_only.csv")
  #> Filtering to requested samples in 'ocean'...
  #> Filtering to requested samples in 'animal'...
  #> Error: None of the samples in 'samplelist' were found in the dataset.
  ```
  
  Request only the collections your list actually covers.

## `use.cache`: reusing downloads, and forcing a refresh

Downloaded files are kept under `tools::R_user_dir("metalog", "cache")` — on
Linux that is `~/.cache/R/metalog`. Check the location on your own machine with:

```r
tools::R_user_dir("metalog", "cache")
#> [1] "/home/you/.cache/R/metalog"
```

Metalog publishes each file under a dated name and serves `..._latest...` as a
redirect to the current one. The function follows that redirect and saves the
file under its **resolved, dated** name, for example
`ocean_metaphlan4_2026-07-07.tsv.gz`.

**With `use.cache = TRUE` (the default) the server is never contacted if a
matching file is already cached.** The cache lookup happens before any network
request: the function globs the cache directory for a dated file matching the
requested name and takes the newest one it finds *locally*. A newer release on
the server is therefore **not** picked up automatically, however long the cached
copy has been sitting there.

Pass `use.cache = FALSE` to bypass that check and pull whatever `latest`
currently resolves to:

```r
tse <- fetchMetalogTSE("ocean", use.cache = FALSE)
#> Fetching file from: https://metalog.embl.de/static/download/profiles/ocean_metaphlan4_latest.tsv.gz
#> Intercepted redirect. Forcing HTTPS: https://metalog.embl.de//static/download/profiles/ocean_metaphlan4_2026-07-07.tsv.gz
#> Downloading to: /home/you/.cache/R/metalog/ocean_metaphlan4_2026-07-07.tsv.gz
```

The refresh sticks: the new dated file is added *alongside* the old one, and
because the cache lookup picks the newest local match, later default calls use
the refreshed copy.

```r
list.files(tools::R_user_dir("metalog", "cache"), pattern = "ocean_metaphlan4")
#> [1] "ocean_metaphlan4_2026-06-07.tsv.gz" "ocean_metaphlan4_2026-07-07.tsv.gz"
```

Nothing is ever evicted, so refreshing repeatedly accumulates snapshots. Delete
files you no longer need by hand.

One consequence worth knowing: because profiles and metadata are cached
independently, a partial refresh can leave you combining a profile from one
snapshot with metadata from another. This is legitimate, and recorded, but not
what you usually want. The provenance entry is where you would notice it:

```r
metadata(tse)$metalog$date_profile
#>        ocean 
#> "2026-07-07" 
metadata(tse)$metalog$date_metadata
#>        ocean 
#> "2026-06-07"
```

If those disagree and you would rather they did not, re-run the call with
`use.cache = FALSE`.

## `make.dense`: dense vs sparse assays

The assay is built sparsely to keep memory down while collections are loaded and
merged. By default it is converted to a dense base R matrix as a final step, for
maximum compatibility with downstream tools. Set `make.dense = FALSE` to keep the
sparse representation:

```r
tse_sp <- fetchMetalogTSE("ocean", make.dense = FALSE)
class(assay(tse_sp, "relabundance"))
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```

Microbiome profiles are mostly zeros, so the difference is substantial even for a
small collection — 306.1 Mb dense against 17.1 Mb sparse for ocean, roughly an
eighteenfold saving:

```r
format(object.size(assay(tse, "relabundance")), units = "MB")     # dense
#> [1] "306.1 Mb"
format(object.size(assay(tse_sp, "relabundance")), units = "MB")  # sparse
#> [1] "17.1 Mb"
```

Keep the default unless memory is tight or you are working with a large
collection; a dense matrix costs `nrow * ncol * 8` bytes whatever its contents,
which grows quickly for `"human"` or for several merged collections.

## Provenance

Every object records how it was built, under `metadata(tse)$metalog`:

```r
metadata(tse)$metalog
#> $source
#> [1] "https://metalog.embl.de/"
#> 
#> $license
#> [1] "Open Database License (ODbL) v1.0"
#> 
#> $collection
#> [1] "ocean"
#> 
#> $meta.type
#> [1] "core"
#> 
#> $date_profile
#>        ocean 
#> "2026-07-07" 
#> 
#> $date_metadata
#>        ocean 
#> "2026-06-07" 
#> 
#> $date_fetched
#> [1] "2026-08-13"
```

`date_profile` and `date_metadata` are parsed from the dated file names and are
named per collection, so they pin the exact snapshot behind every part of the
object. `date_fetched` is simply when the call ran. Report these alongside any
analysis: Metalog is a living database, and the snapshot dates are what make a
result reproducible.

Metalog data is released under the Open Database License
([ODbL](https://opendatacommons.org/licenses/odbl/)) v1.0.

## A short downstream example

Once you have a TSE the whole `mia` toolkit applies. For example, collapse
features to genus level:

```r
library(mia)

genus <- agglomerateByRank(tse, rank = "Genus")
dim(genus)
#> [1] 3786 5086

head(rownames(genus), 5)
#> [1] "g__Abiotrophia"  "g__Abyssibacter" "g__Abyssicoccus" "g__Achromatium" 
#> [5] "g__Achromobacter"
```

See the [Orchestrating Microbiome Analysis
(OMA)](https://microbiome.github.io/OMA/) book for a comprehensive guide to
working with these objects.

## Caching and troubleshooting

**What gets cached.** Profiles, metadata, the MetaPhlAn4 clade table
(`metaphlan4_clades.tsv.gz`) and the SGB tree
(`mpa_vJun23_CHOCOPhlAnSGB_202307.nwk`) all land in the same directory. The clade
table and tree are shared across collections, so they are downloaded once. The
tree comes from the MetaPhlAn database host rather than from Metalog.

**Reclaiming space, or starting clean.** Delete the cache directory; the next
call refetches what it needs.

```r
unlink(tools::R_user_dir("metalog", "cache"), recursive = TRUE)
```

**Redirect messages.** `Intercepted redirect. Forcing HTTPS: ...` is expected.
The server downgrades some redirects to plain HTTP, and the function rewrites
them back to HTTPS before downloading.

**Warnings you are likely to see.** These are informational, not failures:

```
#> Warning messages:
#> 1: 35 of 7896 taxa could not be matched to a tip in the SGB tree. 
#> 2: 35 row(s) couldn't be matched to the tree and are/is removed.
```

Features with no tip in the MetaPhlAn4 SGB tree are **dropped from the object**,
which is why the ocean collection returns 7861 rows rather than the 7896 read
from the profile. Two related warnings can also appear: `n of N taxa could not be
mapped to a full lineage` (the SGB is absent from the clade table, and its
`rowData` ranks are `NA`) and `n sample(s) present in assay data but missing from
metadata` (those samples get `NA` metadata but are kept).

**Errors.**

| Message                                                          | Cause                                                                                                          |
| ---------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `None of the samples in 'samplelist' were found in the dataset.` | The list has no `sample_alias` column, or one of the requested collections contains none of the listed samples |
| `Duplicate sample aliases across collections: ...`               | Two collections share a sample alias; they cannot be merged                                                    |
| `No taxa matched any tip in the SGB tree.`                       | The tree failed to load usefully, e.g. a truncated cached download — delete it and refetch                     |
| `Package 'x' is required.`                                       | A dependency is missing; see [Setup](#setup)                                                                   |

## Reference

Data is downloaded from the [Metalog database](https://metalog.embl.de/) and made
available under the Open Database License (ODbL) v1.0. Please cite Metalog when
using these data.
