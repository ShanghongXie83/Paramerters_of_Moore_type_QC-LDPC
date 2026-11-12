"""Computer certificates for Proposition 26, Lemma 27, and Theorem 30.

Run python verify.py to generate results.json.
Run python verify.py --check to recompute and compare all certificates.
Dependencies: numpy and sympy.

The program uses h=j+1 and delta=k-j. The integers m_u are the root-
multiplicity differences of Lemma 27, represented modulo p by q(U).
The PDF proof and result interpretation are in
Supplementary_Material_for_Proposition_26.pdf.
"""
from collections import Counter, deque
from functools import lru_cache
from fractions import Fraction
from math import comb, gcd, lcm
from pathlib import Path
import argparse
import json
import numpy as np
import sympy as sp


def symbolic_checks():
    z, y = sp.symbols('z y', nonzero=True)
    j = sp.symbols('j', integer=True, positive=True)
    f = z**j*y-y**(j+1)
    hessian = sp.diff(f,z,2)*sp.diff(f,y,2)-sp.diff(f,z,y)**2
    expected = -j**2*z**(j-2)*(z**j+(j**2-1)*y**j)
    assert sp.simplify(hessian-expected) == 0
    f = z**(j+1)-y**(j+1)
    hessian = sp.diff(f,z,2)*sp.diff(f,y,2)-sp.diff(f,z,y)**2
    assert sp.simplify(hessian+j**2*(j+1)**2*z**(j-1)*y**(j-1)) == 0
    h, lam = sp.symbols('h lam', nonzero=True)
    matrix = sp.Matrix([[1,0,0],[-lam/h,0,-1/h],[0,-1/h,lam]])
    assert sp.factor(matrix.det()) == -1/h**2

    a,b,c,d,e,f,u,x = sp.symbols('a b c d e f u x')
    form = sp.Matrix([[a,b,d],[b,c,e],[d,e,f]])
    aa = a*u*u+2*b*u+c
    bb = a*u*x+b*x+d*u+e
    cc = a*x*x+2*d*x+f
    line = sp.Matrix([1,-u,-x])
    dual = (line.T*form.adjugate()*line)[0]
    assert sp.expand(bb*bb-aa*cc+dual) == 0
    r1,r2,s1,s2 = sp.symbols('r1 r2 s1 s2')
    assert sp.expand((r1*s2+r2*s1)**2-4*r1*r2*s1*s2-(r1*s2-r2*s1)**2) == 0
    J, r = sp.symbols('J r', integer=True)
    merged_square_sum = (J-r+1)**2+r-1+(r+1)**2+J-r-1
    assert sp.expand(J**2+J+2-merged_square_sum-2*(r-1)*(J-r-1)) == 0
    return {
        'common_square_bound_gap': '2*(r-1)*(J-r-1)',
        'Hessian_of_z_power_j_times_y_minus_y_power_j_plus_one': 'proved symbolically with j indeterminate',
        'Hessian_of_root_form': 'proved symbolically with j indeterminate',
        'catalectic_minor': '-1/h**2',
        'tangency_discriminant_equals_negative_dual_form': True,
        'rank_two_restriction_discriminant': True,
    }



def prime_divisors(n):
    out = []
    q = 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out


def primitive(p):
    factors = prime_divisors(p - 1)
    return next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in factors))


def subgroup(p, k):
    a = pow(primitive(p), (p - 1) // k, p)
    return [pow(a, t, p) for t in range(k)]


@lru_cache(None)
def partitions(n, lower, cap):
    if n == 0:
        return ((),)
    out = []
    for first in range(lower, n + 1):
        for tail in partitions(n - first, first, cap):
            if tail.count(first) < cap:
                out.append((first,) + tail)
    return tuple(out)


def moment_screen():
    pair_count = multiset_count = 0
    by_delta = Counter()
    survivors = []
    for delta in range(1, 6):
        for j in range(max(3, delta + 1), 6 * delta):
            pair_count += 1
            k, h = j + delta, j + 1
            for S in range(1, h + 1):
                parts = partitions(S, 1, delta)
                for pos in parts:
                    for neg in parts:
                        if not h <= len(pos) + len(neg) <= k:
                            continue
                        multiset_count += 1
                        by_delta[delta] += 1
                        values = pos + tuple(-x for x in neg) + (0,) * (k - len(pos) - len(neg))
                        # Lemma 27: integer square-sum bound and parity.
                        square_sum = sum(x*x for x in values)
                        assert 0 < square_sum <= j*j+3*j+4
                        assert square_sum % 2 == 0
                        G = 0
                        for v in range(2, (k - 1) // delta + 1):
                            G = gcd(G, abs(sum(x ** v for x in values)))
                        for p in prime_divisors(G):
                            if p > 2 * k and (p - 1) % lcm(j, k) == 0:
                                survivors.append((p, j, k, values))
    return {
        'parameter_pairs': pair_count,
        'integer_multisets': multiset_count,
        'multisets_by_delta': dict(by_delta),
        'moment_survivors': survivors,
    }


def polynomial_screen(survivors):
    groups = {}
    for p, j, k, values in survivors:
        groups.setdefault((p, j, k), set()).add(tuple(sorted(x % p for x in values)))
    output = []
    for (p, j, k), targets in sorted(groups.items()):
        delta = k - j
        assert all(0 in target for target in targets)
        A = subgroup(p, k)
        # q(1)=0: q(X)=sum_{i=2}^delta c_i (X^i-X).
        # Domain multiplication by an A-root permutes evaluations; it does
        # not rescale their integer values or alter their multiplicities.
        evaluation_basis = np.array([
            [(pow(u, i, p) - u) % p for u in A]
            for i in range(2, delta + 1)
        ], dtype=np.int64)
        tested = p ** (delta - 1)
        matches = []
        batch_size = 32768
        for start in range(0, tested, batch_size):
            ids = np.arange(start, min(start + batch_size, tested), dtype=np.int64)
            remaining = ids.copy()
            coefficients = np.empty((len(ids), delta - 1), dtype=np.int64)
            for i in range(delta - 1):
                coefficients[:, i] = remaining % p
                remaining //= p
            evaluations = (coefficients @ evaluation_basis) % p
            evaluations.sort(axis=1)
            keep = np.zeros(len(ids), dtype=bool)
            for target in targets:
                keep |= np.all(evaluations == np.array(target, dtype=np.int64), axis=1)
            for coefficients_tail in coefficients[keep]:
                tail = [int(x) for x in coefficients_tail]
                matches.append([(-sum(tail)) % p] + tail)
        output.append({
            'p': p, 'j': j, 'k': k, 'target_multisets': len(targets),
            'all_targets_contain_zero': True, 'normalized_polynomials_tested': tested,
            'matching_coefficients_q1_to_qdelta': matches,
        })
        print('polynomial screen', output[-1], flush=True)
    return output


def binary_kernel(p, j, k):
    A, B = subgroup(p, k), subgroup(p, j)
    n = p * k
    rows = []
    for lam in B:
        for z in range(p):
            rows.append(sum(1 << (t * p + ((z - lam * u) % p)) for t, u in enumerate(A)))
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    free = [i for i in range(n) if i not in pivots]
    basis = []
    for f in free:
        word = 1 << f
        for pivot in sorted(pivots):
            if (pivots[pivot] & word).bit_count() & 1:
                word ^= 1 << pivot
        assert all((row & word).bit_count() % 2 == 0 for row in rows)
        basis.append(word)
    return n, len(pivots), basis, rows


def enumerate_boundary():
    n, rank, basis, rows = binary_kernel(13, 4, 6)
    assert n == 78 and len(basis) == 29
    split = 20
    size = 1 << split
    low = np.empty(size, dtype=np.uint64)
    high = np.empty(size, dtype=np.uint64)
    low[0] = high[0] = 0
    mask = (1 << 64) - 1
    for i, word in enumerate(basis[:split]):
        count = 1 << i
        low[count:2 * count] = low[:count] ^ np.uint64(word & mask)
        high[count:2 * count] = high[:count] ^ np.uint64(word >> 64)
    enumerator = np.zeros(n + 1, dtype=np.int64)
    witness = None
    for block in range(1 << (len(basis) - split)):
        shift = 0
        for i, word in enumerate(basis[split:]):
            if (block >> i) & 1:
                shift ^= word
        weights = np.bitwise_count(low ^ np.uint64(shift & mask))
        weights += np.bitwise_count(high ^ np.uint64(shift >> 64))
        enumerator += np.bincount(weights, minlength=n + 1)
        if witness is None and np.any(weights == 16):
            index = int(np.flatnonzero(weights == 16)[0])
            witness = int(low[index]) | (int(high[index]) << 64)
            witness ^= shift
        if block % 128 == 0:
            print('boundary blocks', block, '/ 512', flush=True)
    enumerated = int(enumerator.sum())
    distance = next(i for i in range(1, n + 1) if enumerator[i])
    assert enumerated == 1 << len(basis)
    assert witness is not None and witness.bit_count() == 16
    assert all((row & witness).bit_count() % 2 == 0 for row in rows)
    return {
        'p': 13, 'j': 4, 'k': 6, 'n': n, 'rank': rank, 'dimension': len(basis),
        'enumerated_words': enumerated, 'minimum_distance': distance,
        'kernel_basis_hex': [hex(word) for word in basis],
        'weight16_witness_hex': hex(witness),
        'weight_enumerator': {str(i): int(value) for i, value in enumerate(enumerator) if value},
    }


def ball_check():
    p, j, k = 211, 3, 5
    A, B = subgroup(p, k), subgroup(p, j)
    root = ('v', 1, 0)
    distance, parent = {root: 0}, {root: None}
    queue = deque([root])
    non_tree_edges = set()
    while queue:
        node = queue.popleft()
        if distance[node] == 5:
            continue
        kind, label, location = node
        if kind == 'v':
            neighbors = [('c', lam, (location + lam * label) % p) for lam in B]
        else:
            neighbors = [('v', u, (location - label * u) % p) for u in A]
        for neighbor in neighbors:
            if neighbor not in distance:
                distance[neighbor], parent[neighbor] = distance[node] + 1, node
                queue.append(neighbor)
            elif parent[node] != neighbor and parent[neighbor] != node:
                non_tree_edges.add(tuple(sorted((node, neighbor))))
    return {
        'p': p, 'j': j, 'k': k, 'A': A, 'B': B,
        'layer_sizes': dict(sorted(Counter(distance.values()).items())),
        'non_tree_edges': len(non_tree_edges),
        'nodes': [{'node': list(node), 'distance': distance[node],
                   'parent': list(parent[node]) if parent[node] else None}
                  for node in distance],
    }


def coefficient_check():
    # [X^5](1-X^5)/sqrt((1-X)(1-wX)), reduced modulo Phi_5(w).
    central = [Fraction(comb(2 * n, n), 4 ** n) for n in range(6)]
    coefficients = [central[5 - r] * central[r] for r in range(6)]
    coefficients[0] -= 1
    coefficients[0] += coefficients[5]  # w^5=1
    for r in range(4):
        coefficients[r] -= coefficients[4]  # w^4=-(1+w+w^2+w^3)
    reduced = coefficients[:4]
    assert reduced == [Fraction(-165, 256), 0, Fraction(-5, 256), Fraction(-5, 256)]
    # t=w^2+w^3 satisfies t^2+t-1=0.  Vanishing requires t=-33.
    obstruction = (-33) ** 2 - 33 - 1
    assert obstruction == 5 * 211
    return {'reduced_coefficient_w0_to_w3': [str(x) for x in reduced],
            'integer_obstruction': obstruction,
            'prime_factors': prime_divisors(obstruction),
            'only_legal_prime_p_greater_than_10': 211}



def generate():
    result = {'scope': 'Computer certificates for Proposition 26, Lemma 27, and Theorem 30. Complete proof and result interpretation: Supplementary_Material_for_Proposition_26.pdf.', 'proof_references': {'support_proposition': 26, 'coset_lemma': 8, 'even_direction_lemma': 9, 'multiplicity_moment_lemma': 27, 'neighboring_order_theorem': 30, 'girth_bound_proposition': 16, 'critical_support_theorem': 24, 'supplement': 'Supplementary_Material_for_Proposition_26.pdf'}, 'symbolic_checks': symbolic_checks()}
    result['moment_screen'] = moment_screen()
    print('moments', result['moment_screen']['parameter_pairs'],
          result['moment_screen']['integer_multisets'],
          len(result['moment_screen']['moment_survivors']), flush=True)
    assert result['moment_screen']['parameter_pairs'] == 69
    assert result['moment_screen']['integer_multisets'] == 19532
    assert len(result['moment_screen']['moment_survivors']) == 43
    result['polynomial_screen'] = polynomial_screen(result['moment_screen']['moment_survivors'])
    assert [r['target_multisets'] for r in result['polynomial_screen']] == [2, 5, 36]
    assert [r['normalized_polynomials_tested'] for r in result['polynomial_screen']] == [13, 361, 923521]
    assert sorted(result['polynomial_screen'][0]['matching_coefficients_q1_to_qdelta']) == [[2, 11], [11, 2]]
    assert all(not r['matching_coefficients_q1_to_qdelta'] for r in result['polynomial_screen'][1:])
    result['boundary_code'] = enumerate_boundary()
    print('boundary distance', result['boundary_code']['minimum_distance'], flush=True)
    result['p211_ball'] = ball_check()
    print('ball', result['p211_ball']['layer_sizes'], 'extra edges', result['p211_ball']['non_tree_edges'], flush=True)
    assert result['boundary_code']['rank'] == 49 and result['boundary_code']['minimum_distance'] == 16
    assert list(result['p211_ball']['layer_sizes'].values()) == [1, 3, 12, 24, 96, 192]
    assert result['p211_ball']['non_tree_edges'] == 0
    result['coefficient_check'] = coefficient_check()
    result['status'] = 'passed'
    result['summary'] = {
        'parameter_pairs': result['moment_screen']['parameter_pairs'],
        'integer_multisets': result['moment_screen']['integer_multisets'],
        'moment_survivors': len(result['moment_screen']['moment_survivors']),
        'normalized_polynomials_tested': [r['normalized_polynomials_tested'] for r in result['polynomial_screen']],
        'boundary_length': result['boundary_code']['n'],
        'boundary_rank': result['boundary_code']['rank'],
        'boundary_dimension': result['boundary_code']['dimension'],
        'boundary_words_enumerated': result['boundary_code']['enumerated_words'],
        'boundary_minimum_distance': result['boundary_code']['minimum_distance'],
        'tree_layer_sizes': list(result['p211_ball']['layer_sizes'].values()),
        'tree_non_tree_edges': result['p211_ball']['non_tree_edges'],
        'exceptional_prime': result['coefficient_check']['only_legal_prime_p_greater_than_10'],
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true',
                        help='Recompute all checks and compare with the existing results.json.')
    args = parser.parse_args()
    if not __debug__:
        raise SystemExit('Run without -O or -OO so that all verification assertions are enabled.')
    if not hasattr(np, 'bitwise_count'):
        raise SystemExit('Install NumPy with bitwise_count support: python -m pip install --upgrade numpy')
    destination = Path(__file__).resolve().with_name('results.json')
    if args.check and not destination.is_file():
        raise SystemExit('The comparison file results.json is missing; run without --check to generate it.')
    saved = json.loads(destination.read_text(encoding='utf-8')) if args.check else None
    result = generate()
    serialized = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)+'\n'
    if args.check and saved != json.loads(serialized):
        raise SystemExit('FAILED: regenerated results differ from results.json; the supplied file was preserved.')
    destination.write_text(serialized, encoding='utf-8')
    if args.check:
        print('All regenerated results agree with results.json.', flush=True)
    print('All checks passed. Output: results.json', flush=True)


if __name__ == '__main__':
    main()
