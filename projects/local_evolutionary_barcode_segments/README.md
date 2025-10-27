# Local Evolutionary Barcode Segments (LBS)  Demo

Minimal, **synthetic** demonstration of detecting *local* co-evolutionary signal segments
along the species axis (barcode-like windows), separate from global similarity.

## Overview
Implements a tiny windowed scoring over a toy matrix to illustrate how local signals
can recover a known seed set better than global correlation.

## Scripts
- scripts/window_consensus_demo.py  builds a tiny synthetic matrix and scores windows.

## Input & Output
Input: synthetic matrix (created in code).
Output: esults/top_windows_demo.csv (toy results), optional plot.

## Purpose
Explain the **idea** and show a runnable demo, without any real NPP/GRID data.

 All data here are synthetic and for demonstration only. Original datasets remain confidential.

Version: 2025-10-27 | Status: Demo skeleton
