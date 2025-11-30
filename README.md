# BLAST101 — CLI & Unittest Suite

This repository is based on a university assignment. The brief required designing a **simple, well‑documented command‑line interface (CLI)** and an **automatic `unittest` suite** for a BLAST/Smith–Waterman/Stats tool. The CLI acts as a gatekeeper: it reads sensible defaults from `settings.ini`, allows safe overrides, and validates inputs early so failures are clear and fast. The test suite (Python `unittest`) covers boundary cases, error paths, and end‑to‑end runs to ensure reproducibility. The full write‑up is available under **`docs/Final_Report.pdf`** (and the original brief under **`docs/Assignment_Brief.pdf`**). This submission was awarded a **Distinction**.

---

## Project goals
- Provide a **gatekeeper CLI** with clear, actionable errors
- Make the tool **reproducible** (deterministic defaults via a settings file)
- Prove correctness with **automated tests** (Python `unittest`)
- Demonstrate **end‑to‑end** runs that produce result files

---


## CLI design & rationale

* **Defaults → Parser → Overrides.** `parse_args()` pulls defaults from `settings.ini`, then safely applies user overrides.
* **Best-practice UX.** Options mirror familiar NCBI BLAST parameters; `--help` explains each flag with examples.
* **Validation as first-class logic.** Each argument has a type + dedicated validator; the CLI transforms and checks values *before* handing them to the pipeline (see *§1.2–1.3*). 
* **Reproducibility.** `current_arguments()` prints the resolved settings at run-time so results are traceable (example output on *p.3*). 


---


## CLI overview

**Entry point:** `main.py`

**Typical usage**

```bash
# 1) BLAST using defaults from settings.ini
python3 main.py

# 2) BLAST overriding defaults
python3 main.py --db custom_db.fasta --query custom_query.fasta --method blast
```

**Core arguments (typical set)**
- `--method {blast, sw, stats}` – choose pipeline
- `--db PATH` / `--query PATH` – FASTA files
- `--out PATH` – JSON or TSV, inferred by extension
- Optional knobs (examples): `--blosum`, `--gap`, `--word_size`, `--max_alignments`, `--max_scores`, `--max_sw_scores`
- Defaults can be provided via a settings file (e.g., `settings.ini`) so you can run with minimal flags

**Validation highlights**
- Ensures **files exist** and contain valid FASTA sequences
- Rejects **DNA‑only** queries for protein search
- Enforces sensible numeric ranges (BLOSUM, gap, word size, maximums); fails **early** with clear messages and non‑zero exit codes

> Note: exact ranges/choices are enforced in `cli.py`; see tests below for boundary cases.


---


## Tests (unittest)

Run everything:
```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

### What’s covered

| File | Focus | Highlights |
|------|------|------------|
| `tests/test_cli.py` | **CLI & validation** | `--help` exits 0; invalid BLOSUM/gap/word size/maximums; missing/empty/invalid FASTA; DNA‑only query; clear error + non‑zero exit |
| `tests/test_blast.py` | **BLAST internals** | extension logic across matrices; left/right choices; stop conditions; diagonal vs. non‑diagonal matches |
| `tests/test_sw_p.py` | **Smith–Waterman scoring** | controlled matches/mismatches, indels; gap penalties; matrix effects on raw scores |
| `tests/test_blast_search.py` | **End‑to‑end BLAST** | runs full pipeline; asserts exit code 0; writes `logs/BLsearch.csv` with content |
| `tests/test_smith_w_search.py` | **End‑to‑end SW** | runs SW pipeline; asserts exit code 0; writes `logs/SWsearch.csv` with content |
| `tests/test_stats.py` | **Stats workflow** | runs `--method stats`; asserts exit code 0 (plots disabled in test if applicable) |


---


## Tests: strategy, coverage & examples

* **Framework:** Python `unittest` with discovery; uses `unittest.mock` for isolated checks. (*p.4–6 of Report*) 
* **Scope:** **6 test scripts / 27 tests** total (unit + integration). Three are end-to-end: BLAST, SW, and Stats. (*p.4–6 of Report*) 
* **Integration pattern:** simulate full runs via `sys.argv` and assert **exit code 0** plus **non-empty CSV** outputs (e.g., `logs/BLsearch.csv`, `logs/SWsearch.csv`). Code example on *p.6 of Report*. 
* **Unit highlights:**

  * **CLI** validators: files, alphabet, DNA-only query guard, numeric ranges. (*p.7 of Report*)
  * **BLAST**: `extend_diagonal` scores across all BLOSUM tables; `process_blast` logic (e.g., 0 score when <2 matches). (*p.8 of Report*)
  * **Smith–Waterman**: `perform_smith_waterman` scores for match/indel scenarios across matrices and gaps. (*p.9 of Report*) 


---


### Coverage snapshot

Average coverage ≈ **81%**. Per-file examples: `blast_101_search.py 98%`, `cli.py 92%`, `smith_waterman_search.py 100%`, `main.py 56%` (driver code). Full table on *p.10 of Report*. 


---

## Repository layout (example)

```
.
├─ Blast101_Engine/
│  ├─ main.py
│  ├─ smith_waterman_p.py
│  ├─ smith_waterman_search.py
│  └─ ...               # other engine modules
├─ cli/
│  └─ cli.py            # CLI entry point (argument parsing + validation)
├─ docs/
│  ├─ Assignment_Brief.pdf
│  └─ Report.pdf
├─ tests/
│  ├─ test_cli.py
│  ├─ test_blast.py
│  ├─ test_blast_search.py
│  ├─ test_smith_w_search.py
│  ├─ test_stats.py
│  └─ test_sw_p.py
├─ .gitignore
├─ CITATION.cff
├─ LICENSE
├─ README.md
└─ settings       # runtime defaults

```

---

## Notes / decisions
- **Fail‑fast CLI**: strict validation saves runtime and produces actionable errors
- **Separation of concerns**: CLI parses/validates; algorithms live in modules
- **Tests first**: unit tests for edge cases + full pipeline checks
- **Reproducibility**: settings file + deterministic defaults

---

## License
MIT © 2025