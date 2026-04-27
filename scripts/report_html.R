create_html_render <- function(pipeline_name,
                               pipeline_version,
                               qc_sc_stats,
                               qc_sc_rna_plot,
                               qc_sc_atac_plot,
                               qc_tf_stats,
                               qc_tf_cutoff_plot,
                               qc_tf_scatter_plot,
                               qc_tf_boxplot_plot,
                               out_render_context)
{
    pipeline_info <- paste0(pipeline_name, " v", pipeline_version)
    rmd_render_context <- glue(r"(
---
title: "Single-cell Multi-omics QC Report"
subtitle: "{pipeline_info}"
date: "`r format(Sys.time(), '%d %B %Y -- %A -- %X')`"
output:
    html_document:
        toc: true
        toc_depth: 4
        theme: united
        highlight: tango
---

```{{r setup, include = FALSE}}
knitr::opts_chunk$set(echo = TRUE, fig.align = "center")
library(reactable)
library(sparkline)
library(UpSetR)
```

```{{js, echo = FALSE}}
function formatNumber(num, precision = 1) {{
    const map = [
        {{ suffix: 'T', threshold: 1e12 }},
        {{ suffix: 'B', threshold: 1e9 }},
        {{ suffix: 'M', threshold: 1e6 }},
        {{ suffix: 'K', threshold: 1e3 }},
        {{ suffix: '', threshold: 1 }},
    ];
    const found = map.find((x) => Math.abs(num) >= x.threshold);
    if (found) {{
        const formatted = (num / found.threshold).toFixed(precision) + found.suffix;
        return formatted;
    }}
    return num;
}}

function rangeMore(column, state) {{
    let min = Infinity
    let max = 0
    state.data.forEach(function(row) {{
        const value = row[column.id]
        if (value < min) {{
            min = Math.floor(value)
        }}
        if (value > max) {{
            max = Math.ceil(value)
        }}
    }})

    const filterValue = column.filterValue || min
    const input = React.createElement('input', {{
        type: 'range',
        value: filterValue,
        min: min,
        max: max,
        onChange: function(event) {{
            column.setFilter(event.target.value || undefined)
        }},
        style: {{ width: '100%', marginRight: '8px' }},
        'aria-label': 'Filter ' + column.name
    }})

    return React.createElement(
        'div',
        {{ style: {{ display: 'flex', alignItems: 'center', height: '100%' }} }},
        [input, formatNumber(filterValue)]
    )
}}

function filterMinValue(rows, columnId, filterValue) {{
    return rows.filter(function(row) {{
        return row.values[columnId] >= filterValue
    }})
}}

function rangeLess(column, state) {{
    let min = Infinity
    let max = 0
    state.data.forEach(function(row) {{
        const value = row[column.id]
        if (value < min) {{
            min = Math.floor(value)
        }}
        if (value > max) {{
            max = Math.ceil(value)
        }}
    }})

    const filterValue = column.filterValue || max
    const input = React.createElement('input', {{
        type: 'range',
        value: filterValue,
        min: min,
        max: max,
        onChange: function(event) {{
            column.setFilter(event.target.value || undefined)
        }},
        style: {{ width: '100%', marginRight: '8px' }},
        'aria-label': 'Filter ' + column.name
    }})

    return React.createElement(
        'div',
        {{ style: {{ display: 'flex', alignItems: 'center', height: '100%' }} }},
        [input, formatNumber(filterValue)]
    )
}}

function filterMaxValue(rows, columnId, filterValue) {{
    return rows.filter(function(row) {{
        return row.values[columnId] <= filterValue
    }})
}}
```

---

## 1. Introduction
**Pipeline:** {pipeline_name}

**Version:** {pipeline_version}

**Homepage:** https://github.com/wtsi-hgi/sc_analysis

Synthetic Lineage Project

To learn sequence to expression rules, we need to figure out how TFs interact with DNA and affect gene expression.

---

## 2. Single-Cell Data QC
This section summarises the basic QC metrics of the single-cell multi-omics data.

### 2.1 Basic QC metrics

```{{r, echo = FALSE}}
df <- as.data.frame(read.table("{qc_sc_stats}", header = TRUE, sep = "\t", check.names = FALSE))
min_row <- ifelse(nrow(df) > 10, 10, nrow(df))
reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,
          minRows = min_row, defaultColDef = colDef(minWidth = 150, align = "left"))
```
<br>

### 2.2 QC plots of RNA data

```{{r, echo = FALSE, fig.show = "hold", fig.align = "center", out.height = "50%", out.width = "50%"}}
knitr::include_graphics("{qc_sc_rna_plot}", rel_path = FALSE)
```
<br>

### 2.3 QC plots of ATAC data
```{{r, echo = FALSE, fig.show = "hold", fig.align = "center", out.height = "50%", out.width = "50%"}}
knitr::include_graphics("{qc_sc_atac_plot}", rel_path = FALSE)
```
<br>

---

## 3. TF QC Metrics
This section summarises the basic QC metrics of the TF data.

### 3.1 Basic QC metrics

```{{r, echo = FALSE}}
df <- as.data.frame(read.table("{qc_tf_stats}", header = TRUE, sep = "\t", check.names = FALSE))
min_row <- ifelse(nrow(df) > 10, 10, nrow(df))
reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,
          minRows = min_row, defaultColDef = colDef(minWidth = 150, align = "left"))
```
<br>

### 3.2 QC plots of TF cutoff
```{{r, echo = FALSE, fig.show = "hold", fig.align = "center", out.height = "100%", out.width = "100%"}}
knitr::include_graphics("{qc_tf_cutoff_plot}", rel_path = FALSE)
```
<br>

```{{r, echo = FALSE, fig.show = "hold", fig.align = "center", out.height = "80%", out.width = "80%"}}
knitr::include_graphics("{qc_tf_scatter_plot}", rel_path = FALSE)
```
<br>

```{{r, echo = FALSE, fig.show = "hold", fig.align = "center", out.height = "80%", out.width = "80%"}}
knitr::include_graphics("{qc_tf_boxplot_plot}", rel_path = FALSE)
```
<br>

---

<br>
    )")

    writeLines(rmd_render_context, out_render_context)
}