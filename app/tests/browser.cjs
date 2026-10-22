// Run against the local app after preparing the real, approved datasets.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { execFileSync } = require('node:child_process');
const { chromium } = require('playwright');

process.chdir(path.resolve(__dirname, '../..'));
const expected = JSON.parse(execFileSync('Rscript', ['-e', `
  source("app/config.R")
  source("app/R/data.R")
  configs <- explorer_config(".")
  genes <- c(parent="CD3E", dc="FCER1A", endothelial="PECAM1",
    epithelial="EPCAM", granulocytes="S100A8", macrophages="LYZ",
    nk="NKG7", stromal="COL1A1")
  views <- lapply(explorer_active_views, function(id) {
    config <- configs[[id]]
    if (config$parent) {
      wrong <- config
      wrong$embedding_input <- wrong$input
      stopifnot(inherits(try(open_explorer_data(wrong), silent=TRUE), "try-error"))
    }
    m <- readRDS(file.path(config$cache, "metadata.rds"))
    x <- HDF5Array::TENxMatrix(file.path(config$cache, "expression.h5"), group="expression")
    counts <- table(m$cells$cluster)
    selected <- names(counts)[c(1, length(counts))]
    checked <- m$cells$cluster %in% selected
    ids <- unique(c(genes[[id]], "PECAM1", "ENSSSCG00000039214"))
    describe <- function(gene, checked) {
      v <- as.numeric(x[, match(gene, colnames(x))])
      anchors <- unique(c(which.max(v), head(which(v > 0), 2), head(which(v == 0), 1)))
      list(id=gene, maximum=max(v), total=sum(v), positive=sum(v > 0),
        restricted_total=sum(v[checked]), restricted_positive=sum(v[checked] > 0),
        anchors=lapply(anchors, function(i) list(x=m$cells$umap_1[i],
          y=m$cells$umap_2[i], value=v[i], checked=checked[i])))
    }
    values <- lapply(ids, describe, checked=checked)
    names(values) <- ids
    anchors <- c(1L, nrow(m$cells) %/% 2L, nrow(m$cells))
    list(id=id, n=nrow(m$cells), label=config$label, gene=genes[[id]], parent=config$parent,
      clusters=names(counts), selected=selected, counts=as.list(counts), genes=values,
      ranges=lapply(m$cells[c("umap_1", "umap_2")], range),
      anchors=lapply(anchors, function(i) list(index=i-1L,
        x=m$cells$umap_1[i], y=m$cells$umap_2[i])))
  })
  cat(jsonlite::toJSON(list(registered=I(names(configs)), active=I(explorer_active_views),
    views=views), auto_unbox=TRUE, digits=NA))
`], { env: { ...process.env, RENV_CONFIG_NAMESPACES_CHECK: 'FALSE',
  OPENBLAS_NUM_THREADS: '1', OMP_NUM_THREADS: '1' }, encoding: 'utf8' }));

// The owner approved the cleaned parent at resolution 0.7 with its original UMAP.
const approved = [
  ['parent', '44827588', '0.7', 'clean_annotated_object.rds'],
  ['dc', '44831038', '0.2', 'subclustered.rds'],
  ['endothelial', '44831031', '0.3', 'subclustered.rds'],
  ['epithelial', '44831032', '0.3', 'subclustered.rds'],
  ['granulocytes', '44831033', '0.2', 'subclustered.rds'],
  ['macrophages', '44831034', '0.4', 'subclustered.rds'],
  ['nk', '44831035', '0.3', 'subclustered.rds'],
  ['stromal', '44831036', '0.3', 'subclustered.rds'],
];
assert.deepEqual(expected.active, approved.map(v => v[0]));
assert.equal(expected.views[0].n, 234183);
assert.equal(expected.views[0].clusters.length, 24);
assert.deepEqual(expected.registered, approved.map(v => v[0]),
  'Only the cleaned parent and seven approved subset views must be registered');
const artifacts = '.amp/in/artifacts/explorer-clean-sources';
fs.mkdirSync(artifacts, { recursive: true });

(async () => {
  const browser = await chromium.launch({ headless: true,
    args: ['--no-sandbox', '--use-gl=angle', '--use-angle=swiftshader',
      '--enable-unsafe-swiftshader'] });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.goto(process.env.EXPLORER_TEST_URL || 'http://127.0.0.1:3838',
      { waitUntil: 'networkidle' });
    const foregroundCount = n => page.waitForFunction(n =>
      document.getElementById('umap')?.data?.[1]?.x?.length === n, n);
    const snapshot = name => page.screenshot({ path: `${artifacts}/${name}.png`, fullPage: true });
    const view = () => page.evaluate(() => {
      const p = document.getElementById('umap');
      return { x: p._fullLayout.xaxis.range, y: p._fullLayout.yaxis.range };
    });
    const checkSize = size => page.waitForFunction(size =>
      document.getElementById('umap').data.every(trace =>
        Math.abs(trace.marker.size - size) < 1e-8), size);
    const checkView = async before => {
      const after = await view();
      for (const axis of ['x', 'y']) {
        assert.ok(after[axis].every((v, i) => Math.abs(v - before[axis][i]) < 1e-10),
          JSON.stringify({ axis, before, after }));
      }
    };
    const selectGene = async (query, id, n) => {
      const input = page.locator('#controls-gene-selectized');
      await input.press('Backspace');
      await input.pressSequentially(query, { delay: 20 });
      await page.locator(`.selectize-dropdown-content [data-value="${id}"]`).click();
      await page.waitForFunction(id =>
        document.getElementById('umap').data[1].text?.[0]?.includes(id), id);
      await foregroundCount(n);
    };
    const checkExpression = async (gene, restricted = false) => {
      const actual = await page.evaluate(() => {
        const m = document.getElementById('umap').data[1].marker;
        const values = Array.from(m.color);
        return { total: values.reduce((a, b) => a + b, 0),
          positive: values.filter(v => v > 0).length, minimum: m.cmin,
          maximum: m.cmax, ordered: values.every((v, i) => i === 0 || v >= values[i - 1]) };
      });
      const total = restricted ? gene.restricted_total : gene.total;
      assert.ok(Math.abs(actual.total - total) < 1e-7 * Math.max(total, 1));
      assert.equal(actual.positive, restricted ? gene.restricted_positive : gene.positive);
      assert.equal(actual.minimum, 0);
      assert.ok(Math.abs(actual.maximum - gene.maximum) < 1e-12);
      assert.ok(actual.ordered);
      assert.ok(await page.evaluate(anchors => {
        const p = document.getElementById('umap').data[1];
        return anchors.every(a => Array.from(p.x).some((x, i) =>
          Math.abs(x - a.x) < 1e-10 && Math.abs(p.y[i] - a.y) < 1e-10 &&
          Math.abs(p.marker.color[i] - a.value) < 1e-10));
      }, gene.anchors.filter(a => !restricted || a.checked)),
      'Expression must stay attached to the correct cells after sorting');
    };

    for (const [i, fixture] of expected.views.entries()) {
      if (i > 0) {
        await page.locator('#controls-view').locator('..').locator('.selectize-input').click();
        await page.locator(`.selectize-dropdown-content [data-value="${fixture.id}"]`).click();
      }
      await page.waitForFunction(({ id, n }) => {
        const p = document.getElementById('umap');
        return p?.layout?.uirevision === id && p.data[0].x.length === n;
      }, fixture, { timeout: 120000 });
      await foregroundCount(0);
      await checkSize(2);
      assert.equal(await page.locator('input[name="controls-groups"]:checked').count(), 0);
      assert.ok(await page.evaluate(anchors => {
        const p = document.getElementById('umap').data[0];
        return anchors.every(a => Math.abs(p.x[a.index] - a.x) < 1e-10 &&
          Math.abs(p.y[a.index] - a.y) < 1e-10);
      }, fixture.anchors), 'Each view must display its own saved UMAP');
      const fullView = await view();
      for (const [axis, column] of [['x', 'umap_1'], ['y', 'umap_2']]) {
        assert.ok(fullView[axis][0] <= fixture.ranges[column][0] &&
          fullView[axis][1] >= fixture.ranges[column][1], 'View changes must reset zoom');
      }
      const labels = await page.locator('#controls-groups .group-name').allTextContents();
      assert.deepEqual(labels, fixture.clusters, 'Cluster labels must remain numeric without annotations');
      assert.deepEqual(await page.locator('#controls-grouping').evaluate(
        el => Object.keys(el.selectize.options)), ['cluster'], 'Annotation grouping must not be offered');
      const source = await page.locator('#source').textContent();
      assert.ok(source.includes(approved[i][1]));
      assert.ok(source.includes(`harmony_snn_res.${approved[i][2]}`));
      assert.ok(source.includes(approved[i][3]));
      assert.ok(await page.locator('#dataset').textContent().then(s => s.includes(fixture.label)));
      assert.equal(await page.locator('#controls-grouping').locator('..').isVisible(), fixture.parent);
      if (i === 0) {
        assert.equal(await page.locator('#controls-view').inputValue(), 'parent');
        assert.deepEqual(await page.locator('#controls-view').evaluate(
          el => Object.keys(el.selectize.options)), expected.active);
        assert.ok(source.includes('UMAP source: job-44542861 annotated_clustered.rds (barcode matched)'));
        assert.equal(await page.evaluate(() =>
          document.getElementById('umap').data[0].marker.color), '#000000');
        await snapshot('00-default-reference');
      } else {
        // The preceding view leaves gene mode and restriction active.
        assert.equal(await page.locator('#controls-gene').inputValue(), 'PECAM1');
        await page.locator('#controls-restrict').uncheck();
        await foregroundCount(fixture.n);
        await checkExpression(fixture.genes.PECAM1);
      }
      await page.locator('input[name="controls-mode"][value="cluster"]').check();
      await foregroundCount(0);
      for (const cluster of fixture.selected) {
        await page.locator(`input[name="controls-groups"][value="${cluster}"]`).check();
      }
      const selectedCount = fixture.selected.reduce((n, id) => n + fixture.counts[id], 0);
      await foregroundCount(selectedCount);
      assert.equal(await page.evaluate(() => new Set(
        document.getElementById('umap').data[1].marker.color).size), 2);
      await page.locator('#controls-view').locator('..').scrollIntoViewIfNeeded();
      if (i === 0) await snapshot('01-parent-clusters');
      await page.locator('input[name="controls-mode"][value="gene"]').check();
      await selectGene(fixture.gene, fixture.gene, fixture.n);
      await checkExpression(fixture.genes[fixture.gene]);
      await snapshot(`${fixture.id}-gene`);
      await page.locator('#controls-restrict').check();
      await foregroundCount(selectedCount);
      await checkExpression(fixture.genes[fixture.gene], true);
      if (i === 0) await snapshot('02-parent-restricted');
      await page.getByRole('button', { name: 'Clear', exact: true }).click();
      await foregroundCount(0);
      await page.getByRole('button', { name: 'Select all', exact: true }).click();
      await foregroundCount(fixture.n);
      await checkExpression(fixture.genes[fixture.gene]);
      await page.locator('#controls-restrict').uncheck();
      await selectGene('PECAM1', 'PECAM1', fixture.n);
      await checkExpression(fixture.genes.PECAM1);
      const plot = await page.locator('#umap .nsewdrag').boundingBox();
      await page.mouse.move(plot.x + plot.width / 2, plot.y + plot.height / 2);
      await page.mouse.wheel(0, -1200);
      await page.waitForTimeout(500);
      const zoomed = await view();
      assert.ok(zoomed.x[1] - zoomed.x[0] < fullView.x[1] - fullView.x[0]);
      if (!fixture.parent) await page.waitForFunction(() =>
        document.getElementById('umap').data[0].marker.size > 2);
      const zoomSize = await page.evaluate(() => document.getElementById('umap').data[0].marker.size);
      if (fixture.parent) assert.equal(zoomSize, 2);
      else assert.ok(zoomSize > 2 && zoomSize < 8, 'Subset dots must grow gradually on zoom');
      await checkSize(zoomSize);
      // Test a user-set zoom, not the transient aspect padding of a new plot.
      await selectGene('IL1B', 'ENSSSCG00000039214', fixture.n);
      await checkExpression(fixture.genes.ENSSSCG00000039214);
      await checkView(zoomed);
      await checkSize(zoomSize);
      await selectGene('PECAM1', 'PECAM1', fixture.n);
      await checkExpression(fixture.genes.PECAM1);
      await checkView(zoomed);
      await checkSize(zoomSize);
      await page.getByRole('button', { name: 'Clear', exact: true }).click();
      await page.locator('#controls-restrict').check();
      await foregroundCount(0);
      await checkView(zoomed);
      await checkSize(zoomSize);
      // A much tighter view distinguishes a bounded size from unchecked growth.
      await page.evaluate(async full => {
        const inset = range => {
          const center = (range[0] + range[1]) / 2;
          const half = (range[1] - range[0]) / 200;
          return [center - half, center + half];
        };
        await Plotly.relayout(document.getElementById('umap'), {
          'xaxis.range': inset(full.x), 'yaxis.range': inset(full.y),
        });
      }, fullView);
      await checkSize(fixture.parent ? 2 : 8);
      await page.getByRole('button', { name: 'Reset view', exact: true }).click();
      // Plotly may add different aspect-ratio padding after an explicit reset.
      await page.waitForFunction(ranges => {
        const layout = document.getElementById('umap')._fullLayout;
        return layout.xaxis.range[0] <= ranges.umap_1[0] &&
          layout.xaxis.range[1] >= ranges.umap_1[1] &&
          layout.yaxis.range[0] <= ranges.umap_2[0] &&
          layout.yaxis.range[1] >= ranges.umap_2[1];
      }, fixture.ranges);
      await checkSize(2);
      console.log(`PASS: ${fixture.id}, ${fixture.n} cells, saved UMAP, cluster IDs, expression and zoom.`);
    }
    // A subset callback must not survive a return to the full parent.
    await page.locator('#controls-view').locator('..').locator('.selectize-input').click();
    await page.locator('.selectize-dropdown-content [data-value="parent"]').click();
    await page.waitForFunction(() => document.getElementById('umap').layout.uirevision === 'parent');
    await checkSize(2);
    assert.equal(await page.locator('.shiny-output-error').count(), 0);
    assert.deepEqual(errors, []);
    console.log('PASS: exactly the cleaned parent and seven approved subset views are registered and selectable.');
    console.log('PASS: numeric-only clusters, gene/alias search, restricted and empty selections, fixed scales, subset-only bounded dot sizing, zoom/reset, no UI errors.');
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
