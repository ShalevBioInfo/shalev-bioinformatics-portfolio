# Portfolio Guide — Computational Biology Projects  
**Author:** Shalev Yaacov  
**Purpose:** This document defines the complete workflow, structure, and rules for building, updating, and maintaining my GitHub portfolio.

---

## 1. Core Philosophy (What and Why)
The portfolio is **not a live workspace** — it’s a clean, curated archive showing stable, well-documented versions of computational biology projects.  
Each upload represents a **snapshot of maturity**, not work in progress.

**Why:**  
- To demonstrate professional coding and scientific thinking.  
- To protect confidential lab data.  
- To show a clear developmental path over time.

---

## 2. Repository Structure
```

shalev-bioinformatics-portfolio/
│
├── README.md               → General overview of the entire portfolio
├── LICENSE                 → MIT License (for code only)
├── .gitignore              → Ignore cache/temp/private files
├── VERSION_LOG.md          → Manual changelog for updates
├── PORTFOLIO_GUIDE.md      → This document (the rulebook)
│
├── projects/
│   └── project_name/
│       ├── scripts/        → Python/R scripts (clean, runnable)
│       ├── results/        → Synthetic demo outputs (png, csv)
│       └── README.md       → Project-specific documentation
│
└── notebooks/              → Small Jupyter/R notebooks for demonstrations

````

**Why:**  
A clear hierarchy makes it easy to review, reproduce, and extend projects without confusion.  
It separates logic (scripts), output (results), and explanations (readme).

---

## 3. Data Policy
**Allowed:**  
- Small synthetic or randomly generated datasets.  
- Simplified toy examples demonstrating methods.  

**Forbidden:**  
- Any real research data from the lab (e.g., NPP/GRID matrices, gene lists, patient data).  

Every project must clearly state:  
> “All data in this repository are synthetic and for demonstration purposes only.”

**Why:**  
To maintain full research confidentiality and prevent data leaks.

---

## 4. Definition of “Ready for Upload”
A project version is considered *ready to upload* only if all the following are true:
1. The code runs end-to-end without errors (on demo data).  
2. The results folder contains at least one example output (PNG, CSV, etc.).  
3. The README.md is complete and explains inputs, outputs, and purpose.  
4. Filenames are clear, consistent, and English-only.  
5. The data are fully synthetic (no sensitive content).  

**Why:**  
To ensure every upload represents something complete, understandable, and safe to share.

---

## 5. Workflow — Creating a New Project
1. **Create structure:**  
   `projects/NEW_PROJECT/{scripts, results, README.md}`  
2. **Write code:**  
   Place clean `.py` (or `.R`) scripts in `scripts/`.  
3. **Generate demo data:**  
   Use random/synthetic values to mimic the expected input.  
4. **Run & Save results:**  
   Store demo outputs (CSV, PNG, etc.) in `results/`.  
5. **Document:**  
   Complete the project `README.md` using the provided template below.  
6. **Update the version log:**  
   Add a line in `VERSION_LOG.md` describing the change.  
7. **Commit manually:**  
   Upload manually via GitHub (no automation).  

**Why:**  
This manual flow ensures that only final, verified, and documented code versions are exposed.

---

## 6. Workflow — Updating an Existing Project
1. Add new script(s) or modify existing ones in `scripts/`.  
2. Generate and save new demo results in `results/`.  
3. Update the project README.md (explain what changed).  
4. Append a new line to `VERSION_LOG.md` (with date + short summary).  
5. Commit the update manually.

**Why:**  
This creates an auditable history and keeps every project self-contained and reproducible.

---

## 7. Upload Checklist (2-minute pre-upload review)
Before every commit, confirm:  
- [ ] No real or sensitive data.  
- [ ] The code runs cleanly.  
- [ ] At least one example result exists.  
- [ ] README.md is complete and formatted.  
- [ ] `VERSION_LOG.md` was updated.  
- [ ] Filenames are clean and English-only.  

**Why:**  
Prevents mistakes and ensures portfolio consistency.

---

## 8. Project README Template
Every project must follow the same structure:

```markdown
# Project Title (Demo)

Brief description of what this project demonstrates using **synthetic data**.

## Overview
Explain what this analysis shows — biological idea, computational logic, or example workflow.

## Scripts
List and explain all included scripts:
- `main_analysis.py` — generates a small demo matrix and summary.
- `heatmap_generator.R` — plots a synthetic heatmap.

## Input & Output
Input: small synthetic dataset or data created in code.  
Output: example files saved under `results/` (e.g., `example_heatmap.png`, `demo_matrix.csv`).

## Example Usage
(Optional) short example command or notebook snippet.

## Purpose
Why this code exists and what biological/computational concept it demonstrates.

⚠️ All data here are synthetic and for demonstration purposes only. Original datasets remain confidential.

Version: YYYY-MM-DD | Status: Demo
````

**Why:**
Consistency makes your portfolio easier to read and maintain — every project “feels” the same.

---

## 9. Commit Message Style

Keep messages short, English-only, and action-based:

Examples:

```
setup: add initial project skeleton
analysis: add main_analysis.py demo
docs: update project README
results: add synthetic heatmap example
```

**Why:**
Clear commit messages tell the story of development at a glance.

---

## 10. File Naming Rules

* Lowercase English only.
* Use underscores (`_`) or hyphens (`-`) for clarity.
* No spaces, no Hebrew, no special characters.
  Examples:
  `main_analysis.py`, `demo_matrix.csv`, `example_heatmap.png`

**Why:**
Prevents errors across operating systems and keeps naming uniform.

---

## 11. When to Start a New Project

Start a new project folder in `projects/` when:

* The analysis or idea represents a distinct research concept (e.g., new method, new data type).
* The code and outputs can stand alone as a self-contained example.

**Why:**
Keeps each topic modular, readable, and shareable.

---

## 12. Moving from Private to Public

When the portfolio is ready to share:

1. Review every README and results folder for sensitive data.
2. Ensure clarity and professional documentation.
3. Go to **Settings → Change visibility → Public**.
4. Confirm the MIT License is in place.

**Why:**
Allows public sharing while keeping full data protection.

---

## 13. Maintenance Routine (Once a Month)

* Quickly review all projects — improve clarity, update outputs if needed.
* Add new projects only when they reach “ready for upload” status.
* Keep `VERSION_LOG.md` updated for all changes.

**Why:**
Regular review keeps the portfolio alive and professional.

---

## 14. Core Files Summary

| File                 | Purpose                                |
| -------------------- | -------------------------------------- |
| `README.md`          | Overall portfolio description          |
| `LICENSE`            | Legal permission for code reuse        |
| `.gitignore`         | Ignore temp or private files           |
| `VERSION_LOG.md`     | Manual changelog                       |
| `PORTFOLIO_GUIDE.md` | This rulebook (main reference)         |
| `projects/`          | Main projects (scripts + demo results) |
| `notebooks/`         | Small demos or testing notebooks       |

---

### In short:

Each project =
**scripts (code)** + **results (demo outputs)** + **README (documentation)**.
Each update =
**new code/output + updated README + new line in VERSION_LOG.md + clean commit.**

```

