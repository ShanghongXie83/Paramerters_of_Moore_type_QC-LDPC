# Supplementary Proof and Computation for Moore-Type QC-LDPC Codes

This repository contains the supplementary proof and reproducible computations for **Exact Parameter Classes and Bounds for Moore-Type QC-LDPC Codes with Odd Prime Size**.

The files support the structural result for codeword supports of weight `2j+2` and the finite computations used in the neighboring-orbit-order minimum-distance bound.

Repository: [ShanghongXie83/Paramerters_of_Moore_type_QC-LDPC](https://github.com/ShanghongXie83/Paramerters_of_Moore_type_QC-LDPC).

**Files**

| File | Purpose |
| --- | --- |
| [Supplementary_Material_for_Proposition_25.pdf](Supplementary_Material_for_Proposition_25.pdf) | Complete proof of the structural proposition (Proposition 26 in manuscript v3.2); the final section explains the computation output and its role in the proofs. The PDF filename is retained across versions. |
| [verify.py](verify.py) | The single program performing all symbolic checks and finite computations. |
| [results.json](results.json) | The supplied program-generated results, including a summary and complete finite certificates. |
| [README.md](README.md) | Installation and reproduction instructions. |

**Installation and reproduction**

Download the files from the same repository commit, keep them together, and open a terminal in that directory.

Install the two Python dependencies:

```bash
python -m pip install numpy sympy
```

To recompute everything and compare with the supplied `results.json`, run:

```bash
python verify.py --check
```

The comparison file is used only to check the regenerated output, not as mathematical input. A mismatch stops the program and preserves the supplied file.

To generate results without comparing against the supplied file, run:

```bash
python verify.py
```

This writes or replaces `results.json` next to the program. All computations start from the definitions and parameters in the program.

No other script, local data file, Magma installation, or network connection is needed after the dependencies have been installed. Use an ordinary Python interpreter without `-O` or `-OO`, since the assertions are mathematical verification checks. NumPy must provide the `bitwise_count` operation; the program reports an actionable error if it is absent.

The full run enumerates all 536870912 words of the boundary code. It uses bounded integer-array batches and does not store all these words at once. Runtime depends on the machine. Progress appears in the terminal, while all mathematical results are collected in the one JSON file. Successful completion prints `All checks passed`. In comparison mode, any mismatch stops the program and preserves the supplied output file.

**Output and verification scope**

The output contains:

- `status` and `summary`: completion status and the principal results.
- `symbolic_checks`: the Hessian identities with an indeterminate order, the nonzero coefficient minor, and the conic restriction identities.
- `moment_screen`: 69 parameter pairs, 19532 integer multisets, and the 43 surviving candidates.
- `polynomial_screen`: the three complete normalized polynomial searches and their matching coefficients.
- `boundary_code`: the reconstructed binary rank, kernel basis, complete weight distribution, and a minimum-weight word.
- `coefficient_check`: the exact coefficient obstruction and its unique admissible prime.
- `p211_ball`: the full radius-five tree, with node labels, distances, and parent relations.

Read the final section of the PDF for the mathematical meaning and proof role of these fields. The structural proposition has a complete written proof; the symbolic calculations check its algebraic identities, and the finite computations close the neighboring-order distance proof. Numerical examples, numerical tables, and conjecture scans are outside this package's verification scope.

**Referencing this material**

When referring to these files in a manuscript, use a link to the specific repository commit containing the corresponding supplementary material. Use the same commit for the PDF, program, and result file so that the cited proof and computations remain tied to one version. The repository homepage provides access to the latest version.
