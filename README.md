# Supplementary Material for Proposition 26

This package accompanies Proposition 26 of *Exact Parameter Classes and Bounds for Tanner (j,k) QC-LDPC Codes*. It contains the complete proof of the proposition and the computer certificates used in Theorem 30.

## Files

| File | Contents |
| --- | --- |
| `Supplementary_Material_for_Proposition_26.pdf` | Proposition 26, its complete proof, and a section explaining the computation output. |
| `verify.py` | The single verification program. |
| `results.json` | The complete program output, including the numerical certificates. |
| `README.md` | Package description and reproduction instructions. |

## Correspondence with the paper

- Lemma 8: equality of finite product-set cardinalities and the resulting coset structure.
- Lemma 9: the even-direction bound and the associated Redei-polynomial construction.
- Proposition 26, parts (i)--(iii): odd column occupancy, the quadratic quotient, and the two root polynomials.
- Lemma 27: common root-multiplicity and moment constraints, including the square-sum bound.
- Theorem 30: the minimum-distance bound for neighboring row and column orders.
- Proposition 16: the girth-based distance bound used at the exceptional prime.

The PDF recalls the statement of Lemma 9 and its polynomial notation. Its proof is given in the paper. Proposition 26 and its part labels agree with the paper.

## Running the verification

Install the dependencies:

```text
python -m pip install numpy sympy
```

Generate the output:

```text
python verify.py
```

The program writes `results.json` next to itself. To recompute all certificates and compare them with the supplied output, run:

```text
python verify.py --check
```

Use the ordinary Python execution mode so that assertions are enabled. NumPy provides the `bitwise_count` operation used in the full kernel enumeration. The program prints progress and reports successful completion after every check has passed.

The boundary-code calculation traverses all 536870912 words in bounded integer-array batches. Runtime depends on the machine.

## Output and proof roles

| Field | Result and use |
| --- | --- |
| `status` | Successful completion of the verification. |
| `proof_references` | The paper's proposition, lemma, and theorem numbers and the supplementary PDF filename. |
| `summary` | Principal counts and conclusions. |
| `symbolic_checks` | Hessian, coefficient-minor, and conic-restriction identities for Proposition 26, together with the square-sum identity of Lemma 27. |
| `moment_screen` | 69 parameter pairs, 19532 integer multisets, and 43 candidates passing the prime and moment conditions. |
| `polynomial_screen` | Complete normalized polynomial searches over the three remaining parameter triples. |
| `boundary_code` | The length-78 binary matrix rank, a kernel basis, the complete weight distribution, and a weight-16 codeword. |
| `coefficient_check` | The exact exceptional-prime obstruction 1055 = 5*211. |
| `p211_ball` | All nodes and parent relations of the radius-five tree at (p,j,k) = (211,3,5). |

The normalized polynomial searches test 13, 361, and 923521 polynomials for (p,j,k) = (13,4,6), (19,6,9), and (31,10,15), respectively. Only the first case gives matches, namely q(U) = plus or minus 2(U^2-U). Full enumeration of its binary kernel gives minimum distance 16. At the exceptional prime, the tree layers have sizes 1, 3, 12, 24, 96, 192, with zero non-tree edges. These calculations complete the finite branches of Theorem 30.

The final section of the PDF explains each output field, the binary-coordinate encoding, and the deductions used in the proofs. All computations use exact symbolic expressions, integers, rational numbers, and finite-field residues.

Repository: [Supplementary Material for Proposition 26](https://github.com/ShanghongXie83/Paramerters_of_Moore_type_QC-LDPC).
