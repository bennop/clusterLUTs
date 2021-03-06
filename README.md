---
output:
  html_document: 
    highlight: tango
  pdf_document: default
---
<!-- badges: start -->
[![Maintainability](https://api.codeclimate.com/v1/badges/fc9b78585c15a6d0ed30/maintainability)](https://codeclimate.com/github/bennop/clusterLUTs/maintainability)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis build status](https://travis-ci.com/github/bennop/clusterLUTs.svg?branch=master)](https://travis-ci.com/github/bennop/clusterLUTs)
<!-- badges: end -->

# Cluster Lookup Tables - cluster LUTs

Coloring for a cluster tree with 145 leaves and 28 levels (10 (top), 15, ..., 145 clusters (bottom)):


## Concepts


A *hue range* specifies a range or set of ranges in hue space ([0..1]). In sets the hue ranges need not be disjunct but usually are in the context of this package.

A *tree range* is a hierarchical collection of hue ranges where each level corresponds to one level of clustering. 

A *color matrix* is a 3x*n* matrix specifying RGB colors (like the output of `col2rgb`)

### Generation
- `hue.range.split`
- `tree.ranges`

### Output
- `hue.range.lines`
- `tree.range.plot`

more details later

Benno Pütz 
