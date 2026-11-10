# Supplementary proof and computation

Upload these four files together:

| File | Purpose |
| --- | --- |
| `Supplementary_Material_for_Proposition_25.pdf` | Complete supplementary proof; the final section explains every part of the computation output and how it is used in the proofs. |
| `verify.py` | The single program performing all symbolic checks and finite computations. |
| `results.json` | The single program-generated result file, including a summary and complete finite certificates. |
| `README.md` | These installation and reproduction instructions. |

Install the two Python dependencies:

```text
python -m pip install numpy sympy
```

Run the computation:

```text
python verify.py
```

The program writes `results.json` next to itself. It computes from its definitions and parameters and does not need the supplied result file as input. To recompute everything and compare with the supplied result, run:

```text
python verify.py --check
```

No other script, local data file, Magma installation, or network connection is needed after the dependencies have been installed. Use an ordinary Python interpreter without `-O` or `-OO`, since the assertions are mathematical verification checks. NumPy must provide the `bitwise_count` operation; the program reports an actionable error if it is absent.

The full run enumerates all 536870912 words of the boundary code. It uses bounded integer-array batches and does not store all these words at once. Runtime depends on the machine. Progress appears in the terminal, while all mathematical results are collected in the one JSON file. Successful completion prints `All checks passed`. In comparison mode, any mismatch stops the program and preserves the supplied output file.

The output contains:

- `status` and `summary`: completion status and the principal results.
- `symbolic_checks`: the Hessian identities with an indeterminate order, the nonzero coefficient minor, and the conic restriction identities.
- `moment_screen`: 69 parameter pairs, 19532 integer multisets, and the 43 surviving candidates.
- `polynomial_screen`: the three complete normalized polynomial searches and their matching coefficients.
- `boundary_code`: the reconstructed binary rank, kernel basis, complete weight distribution, and a minimum-weight word.
- `coefficient_check`: the exact coefficient obstruction and its unique admissible prime.
- `p211_ball`: the full radius-five tree, with node labels, distances, and parent relations.

Read the final section of the PDF for the mathematical meaning and proof role of these fields. The structural proposition has a complete written proof; the symbolic calculations check its algebraic identities, and the finite computations close the neighboring-order distance proof. Numerical examples, numerical tables, and conjecture scans are outside this package's verification scope.

Once a GitHub repository is chosen, use its URL in place of `xxxx` in the main manuscript. The upload consists of the four files listed above; the ZIP file is only a delivery container.
