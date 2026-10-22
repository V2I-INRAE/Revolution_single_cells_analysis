# Local cell explorer

The app opens the cleaned parent from
`data/job-44542861/clean_annotated_object.rds`. It displays numeric clusters
from the saved `harmony_snn_res.0.7`, without appending annotation labels.
All 234,183 cells and 21,123 RNA features are retained without refitting.
The cleaned object has no dimensional reductions. The app uses the original `umap_harmony` from
`data/job-44542861/annotated_clustered.rds`, matched by exact cell barcode.
RNA expression and metadata come only from the cleaned object. Cells removed
during cleanup are not displayed, and no UMAP is recomputed.

The cleaned parent is the default view. Seven subset views are also enabled,
using the following sources under `data/job-44542861/subsets/`:

| Subset | Resolution | Source RDS |
|---|---:|---|
| DC | 0.2 | `DC/job-44831038/subclustered.rds` |
| Endothelial | 0.3 | `Endothelial/job-44831031/subclustered.rds` |
| Epithelial | 0.3 | `Epithelial/job-44831032/subclustered.rds` |
| Granulocytes | 0.2 | `Granulocytes/job-44831033/subclustered.rds` |
| Macrophages | 0.4 | `Macrophage_Monocytes/job-44831034/subclustered.rds` |
| NK cells | 0.3 | `NK_cells/job-44831035/subclustered.rds` |
| Stromal | 0.3 | `Stromal/job-44831036/subclustered.rds` |

Each subset uses its own saved `umap_harmony`, normalized RNA and requested
cluster partition, with numeric cluster labels. Nothing is refitted or
annotated in the app. Refit job-44890936 views are no longer configured;
their source files and previous caches remain unchanged.
The parent cache is `app/cache/clean-44542861-original-umap-res0.7/`.

## Prepare once and launch locally

On a Genotoul compute node, run:

```sh
bash /work/project/revo-pig-sc/analysis/app/run.sh
```

The launcher loads GCC 12.2.0 and R 4.6.1 and selects the project directory.
The current views are already prepared. Keep the terminal open while
using the app, and press Ctrl+C to stop it. See the
[cluster connection instructions](../README.md#cell-explorer-on-the-cluster)
for the SSH tunnel from your computer. Only one app can use port 3838 on a node.

When running on the same computer as your browser, open
`http://127.0.0.1:3838`. In an R environment with the project packages,
`Rscript app/run.R` or `shiny::runApp("app")` from the project root also works.

Preparation requires enough memory to read one source Seurat object at a time.
For the parent, it extracts and releases the original UMAP source before
opening the cleaned RNA source. Later launches read the small metadata file and
retrieve one gene at a time from sparse HDF5 storage. Browser memory is
additional.
Each configured view has a separate directory under `app/cache/`. Keep its
`metadata.rds`, `expression.h5`, and `.INFO` files together. Copy the caches
with the code to avoid preparing the Seurat objects again on another local
computer. Source paths stored in
metadata are provenance, not runtime dependencies. Navigation opens only the
selected view and reads one expression vector at a time.

If a cache is missing, load the same GCC/R modules and run
`Rscript app/prepare.R` from the project root to prepare missing active views.
Pass a view ID from `config.R` to prepare only that view.
Preparation keeps completed caches after checking their source identifier,
UMAP source, reduction, partition and cell/gene order. It refuses incomplete
caches and does not overwrite files. If source contents or gene labels change, move the
old cache aside and prepare again explicitly. Preparation status and
saved-data checks are recorded in each cache's `.INFO` file.

Dependencies are already listed in the project's `renv.lock`. Preparation
uses `Seurat`, `Matrix`, and `HDF5Array`. The app uses `shiny`, `bslib`,
`plotly`, `htmlwidgets`, `viridisLite`, and `HDF5Array`; PDF export also uses
the existing `ggplot2` and `ggrastr` packages.

## Navigate

Choose the cleaned parent or a subset under **Dataset / view**.
Labels include the selected resolution. The parent uses the original
resolution 0.7 partition, not the refit clusters. Every view uses numeric
cluster labels only; no biological names are appended to cluster numbers.

All cells initially appear black. Check groups to color them, while
unchecked cells remain black. Changing the view or grouping clears checked
groups. A selected gene is kept across views if it exists in the new view.

Choose **Gene expression** and search a symbol or Ensembl ID. Search covers
all retained features, including the existing reviewed gene-name aliases.
Original feature IDs remain unchanged. Tick **Color checked groups only**
to restrict expression coloring; otherwise every cell receives a color.

Expression is the existing RNA `data` layer, already LogNormalize with
scale factor 10000. The app does not log-transform or normalize again.
Zero appears at the purple end of viridis, and higher values appear yellow.
The color scale uses all cells in the active view for that gene and stays
fixed when groups are selected. It is recalculated when the view changes,
so compare the displayed ranges before comparing colors between views.
High-expression cells are drawn last. A notice appears
if a feature was absent from any original sample layer, since joining layers
fills those cells with zero as in the existing FeaturePlots.

Drag to pan and scroll to zoom. Gene and group changes preserve the view.
Changing datasets resets zoom because the embeddings differ. **Reset view**
restores the active embedding's full extent. All cells in each saved source
are retained; the app does not downsample.
In subset views, dots grow gradually from 2 to at most 8 pixels as you zoom in.
Gene and group changes preserve that size; zooming out or resetting the view
shrinks them again. The full parent keeps fixed 2-pixel dots. Sizing runs in
the browser and does not reload cells or expression from R.

## Save a PDF for editing

Click **Save PDF** next to **Reset view** for the current cluster or
gene-expression plot. The download keeps the dataset, checked groups,
expression restriction, gene color range, zoom and dot size. Reset the view
first if you want the full embedding instead of a close-up.

The PDF uses the same method as the existing marker-panel PDFs: all cell dots
form one 600-dpi raster layer, while titles, axes, cluster legends and gene
color scales remain editable vectors. Open the PDF in Illustrator or Affinity
to edit those elements. Individual cell dots are not separate editable objects.
The page follows the browser plot dimensions; the PDF layout is redrawn from
the same cached values rather than taking a screenshot. The filename includes
the dataset, saved clustering resolution and gene ID or `clusters`.

The existing Plotly camera button still saves PNGs. PDF generation runs only
when requested and does not reopen the source Seurat object or rebuild caches.

## Code organization

* `config.R` lists registered views and the ordered active selection.
* `prepare.R` and `R/prepare_data.R` build and verify the viewing data.
* `R/data.R` opens the cache and retrieves individual expression vectors.
* `R/controls.R` owns the Shiny controls module.
* `R/plot.R` builds WebGL traces and updates the plot without rebuilding it.
* `R/export.R` saves the current plot as a hybrid PDF with editable frame elements.
* `app.R` connects the modules, layout, and status text.
* `www/explorer.css` contains layout styling.

Preparation reuses `viz_palette` from `R/utils/save_plots.R` and
`marker_feature_renames` from `R/clusters/markers.R`. It preserves original
IDs rather than modifying the Seurat object to apply display aliases.

To add another compartment, register its approved RDS, resolution and a new
cache directory in `config.R`, add its ID to `explorer_active_views`, then
prepare that view.
Only the parent specifies a separate `embedding_input`; its provenance is
stored alongside the RNA source and displayed in the plot footer. Preparation
rejects missing or duplicate UMAP barcodes before matching coordinates.
The shared controls and plotting modules do not need another implementation.
Existing lower or higher resolutions are not chosen automatically.

## Browser checks

The browser test uses the prepared real dataset, not generated data. With
the app running in a separate terminal, run from the project root:

```sh
npm ci --prefix app/tests
npx --prefix app/tests playwright install chromium
npm test --prefix app/tests
```

The checks cover exactly the cleaned parent and seven approved subset views,
their sources and resolutions, numeric cluster labels without annotation grouping, gene-symbol and
Ensembl-alias search, expression values, restricted and empty selections,
fixed color scales, zoom preservation and reset, and subset-only dot sizing.
They also check saved coordinates and expression-to-cell correspondence.
Screenshots are saved under `.amp/in/artifacts/explorer-clean-sources/`.
Parent UMAP and zoomed DC feature-plot PDF downloads are checked and saved
under `.amp/in/artifacts/explorer-pdf/`.
Set `EXPLORER_TEST_URL` when checking a separate local test instance instead
of the default `http://127.0.0.1:3838`.
