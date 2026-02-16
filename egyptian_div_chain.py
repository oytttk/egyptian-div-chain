#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Egyptian fraction search with adjacent divisibility preference.

Given a/b, write:
    a/b = k + r, 0 <= r < 1
    r = sum_i 1/d_i   (distinct denominators, output strictly increasing)

Hard constraint:
    max(d_i) <= 2 * b_orig
where b_orig is the denominator as typed in "a/b" (if not typed as a/b, use reduced denominator).

Two optimization modes:
  --mode shortest (default):
      1) minimize number of terms
      2) maximize longest consecutive run where d_{i+1} % d_i == 0
      3) maximize number of divisible adjacent edges
      4) minimize max denominator
      5) minimize sum of denominators

  --mode longest:
      1) maximize longest divisible run
      2) maximize divisible edges
      3) minimize number of terms
      4) minimize max denominator
      5) minimize sum of denominators

Optional strict constraint:
  --force-chain : require d_{i+1} % d_i == 0 for ALL adjacent pairs.
"""

from __future__ import annotations

import argparse
import random
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple


# -----------------------------
# Parsing / formatting
# -----------------------------

def parse_fraction_with_raw_den(s: str) -> Tuple[Fraction, Optional[int]]:
    """Parse as Fraction; if s is 'a/b', also return raw denominator b (abs, as typed)."""
    s = s.strip()
    if "/" in s:
        a_str, b_str = s.split("/", 1)
        a = int(a_str.strip())
        b = int(b_str.strip())
        if b == 0:
            raise ValueError("denominator cannot be 0")
        return Fraction(a, b), abs(b)
    return Fraction(s), None


def format_egyptian(denoms: Sequence[int]) -> str:
    return " + ".join(f"1/{d}" for d in denoms) if denoms else "0"


def verify(frac: Fraction, denoms: Sequence[int]) -> bool:
    return sum(Fraction(1, d) for d in denoms) == frac


# -----------------------------
# Divisibility scoring
# -----------------------------

def longest_divisible_run(denoms: Sequence[int]) -> int:
    """Longest consecutive run where d_{i+1} % d_i == 0."""
    if not denoms:
        return 0
    best = 1
    cur = 1
    for i in range(1, len(denoms)):
        if denoms[i] % denoms[i - 1] == 0:
            cur += 1
        else:
            best = max(best, cur)
            cur = 1
    return max(best, cur)


def count_divisible_edges(denoms: Sequence[int]) -> int:
    """Count adjacent pairs where d_{i+1} % d_i == 0."""
    return sum(1 for i in range(1, len(denoms)) if denoms[i] % denoms[i - 1] == 0)


# -----------------------------
# Bounds helpers
# -----------------------------

def _lb_next_den(rem: Fraction) -> int:
    """Lower bound for next unit fraction denominator: ceil(b/a)."""
    a, b = rem.numerator, rem.denominator
    return (b + a - 1) // a


@dataclass(frozen=True)
class SearchConfig:
    max_terms: int
    max_den: int
    time_limit_s: float
    node_limit: int
    range_cap: int
    force_chain: bool
    mode: str  # "shortest" or "longest"


# -----------------------------
# Core search
# -----------------------------

def bounded_best_egyptian_div_chain(frac: Fraction, cfg: SearchConfig) -> Tuple[Optional[List[int]], Dict[str, object]]:
    """
    Search for an Egyptian fraction decomposition of 0<frac<1 under cfg.max_den.
    Returns best_denoms or None if no exact solution found within constraints.
    """
    t0 = time.time()
    nodes = 0

    best: Optional[List[int]] = None
    best_score = None  # tuple, depends on mode

    ub_cache: Dict[Tuple[int, int], Fraction] = {}

    def ub_distinct_sum(d0: int, remaining: int) -> Fraction:
        """
        Upper bound (exact) using smallest distinct increasing denominators:
        d0, d0+1, ..., d0+remaining-1. If cannot fit <= max_den, return 0.
        """
        key = (d0, remaining)
        if key in ub_cache:
            return ub_cache[key]
        if remaining <= 0:
            ub_cache[key] = Fraction(0, 1)
            return ub_cache[key]
        if d0 + remaining - 1 > cfg.max_den:
            ub_cache[key] = Fraction(0, 1)
            return ub_cache[key]
        s = Fraction(0, 1)
        for d in range(d0, d0 + remaining):
            s += Fraction(1, d)
        ub_cache[key] = s
        return s

    def make_score(denoms: List[int]) -> Tuple[int, ...]:
        run = longest_divisible_run(denoms)
        edges = count_divisible_edges(denoms)
        terms = len(denoms)
        mx = max(denoms)
        sm = sum(denoms)
        if cfg.mode == "shortest":
            # primary: fewer terms
            return (-terms, run, edges, -mx, -sm)
        else:
            # primary: longer run
            return (run, edges, -terms, -mx, -sm)

    def update_best(denoms: List[int]):
        nonlocal best, best_score
        sc = make_score(denoms)
        if best_score is None or sc > best_score:
            best = denoms[:]
            best_score = sc

    def gen_candidates(rem: Fraction, start_d: int, remaining: int, denoms: List[int]) -> List[int]:
        """
        Candidate denominators d in [d0, dmax].
        If force_chain and denoms non-empty: only consider multiples of last denom.
        Otherwise: consider interval with bias toward multiples of last denom.
        """
        d0 = max(start_d, _lb_next_den(rem))
        if d0 > cfg.max_den:
            return []

        a, b = rem.numerator, rem.denominator
        dmax = min(cfg.max_den, (remaining * b) // a)
        if dmax < d0:
            return []

        # force full chain: only multiples of last
        if cfg.force_chain and denoms:
            last = denoms[-1]
            first = ((d0 + last - 1) // last) * last
            if first > dmax:
                return []
            cand = list(range(first, dmax + 1, last))
            # prefer smaller denominators first
            return cand[: cfg.range_cap]

        width = dmax - d0 + 1
        if width <= cfg.range_cap:
            cand = list(range(d0, dmax + 1))
        else:
            s = set()
            s.add(d0)
            # near d0
            near = min(250, cfg.range_cap // 10)
            for i in range(1, near):
                v = d0 + i
                if v <= dmax:
                    s.add(v)
            # multiples of last to encourage divisibility
            if denoms:
                last = denoms[-1]
                m = ((d0 + last - 1) // last) * last
                for k in range(0, 600):
                    v = m + k * last
                    if v > dmax:
                        break
                    s.add(v)
            # random sprinkle
            for _ in range(150):
                s.add(random.randint(d0, dmax))
            cand = sorted(s)[: cfg.range_cap]

        # ordering: multiples of last first
        if denoms:
            last = denoms[-1]
            cand.sort(key=lambda d: (0 if d % last == 0 else 1, d))
        else:
            cand.sort()
        return cand

    def dfs(rem: Fraction, start_d: int, remaining: int, denoms: List[int]):
        nonlocal nodes

        if rem.numerator == 0:
            update_best(denoms)
            return
        if remaining == 0:
            return
        if nodes >= cfg.node_limit:
            return
        if (time.time() - t0) > cfg.time_limit_s:
            return

        nodes += 1

        d0 = max(start_d, _lb_next_den(rem))
        if d0 > cfg.max_den:
            return
        # must be able to place remaining distinct denominators <= max_den
        if d0 + remaining - 1 > cfg.max_den:
            return
        # sum upper bound
        if ub_distinct_sum(d0, remaining) < rem:
            return

        for d in gen_candidates(rem, start_d, remaining, denoms):
            if denoms and d <= denoms[-1]:
                continue
            unit = Fraction(1, d)
            if unit > rem:
                continue
            # if force_chain, the candidate generator already enforced it
            dfs(rem - unit, d + 1, remaining - 1, denoms + [d])

            if nodes >= cfg.node_limit or (time.time() - t0) > cfg.time_limit_s:
                return

    # iterative deepening by term count
    for t in range(1, cfg.max_terms + 1):
        dfs(frac, 2, t, [])
        # For shortest mode: once we have any solution at term count t, do not try larger t.
        if cfg.mode == "shortest" and best is not None:
            break
        if nodes >= cfg.node_limit or (time.time() - t0) > cfg.time_limit_s:
            break

    diag = {
        "nodes": nodes,
        "time_s": time.time() - t0,
        "best_score": best_score,
        "cfg": cfg,
    }
    return best, diag


# -----------------------------
# CLI
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Egyptian fractions with adjacent divisibility preference; max denom <= 2*b_orig."
    )
    ap.add_argument("fraction", type=str, help="Input rational, e.g. 13/32 or 325/799")
    ap.add_argument("--max-terms", type=int, default=None, help="Max number of unit fractions for fractional part")
    ap.add_argument("--time", type=float, default=2.0, help="Search time limit (seconds)")
    ap.add_argument("--nodes", type=int, default=800000, help="Search node limit")
    ap.add_argument("--range-cap", type=int, default=12000, help="Candidate cap per depth")
    ap.add_argument("--seed", type=int, default=0, help="Random seed (0 uses time-based seed)")
    ap.add_argument("--force-chain", action="store_true",
                    help="Require ALL adjacent pairs satisfy divisibility (whole list is a chain).")
    ap.add_argument("--mode", choices=["shortest", "longest"], default="shortest",
                    help="Optimization mode: shortest (fewest terms) or longest (longest divisible run first).")
    args = ap.parse_args()

    if args.seed == 0:
        random.seed(int(time.time() * 1e6) & 0xFFFFFFFF)
    else:
        random.seed(args.seed)

    x, raw_den = parse_fraction_with_raw_den(args.fraction)
    if x == 0:
        print("0 = 0")
        return

    sign = -1 if x < 0 else 1
    x = abs(x)

    b_orig = raw_den if raw_den is not None else x.denominator
    max_den = 2 * b_orig

    k = x.numerator // x.denominator
    r = x - k

    print(f"[input] x = {('-' if sign < 0 else '')}{x.numerator}/{x.denominator}")
    print(f"[bound] max_den = {max_den} (=2*b_orig, b_orig={b_orig})")

    if r == 0:
        val = sign * k
        print(f"[result] {val} = {val}/1")
        return

    # max feasible terms (distinct increasing denominators starting from >=2)
    hard_cap = max(0, max_den - 1)
    if args.max_terms is None:
        # conservative default; raise if needed
        max_terms = min(30, hard_cap)
        if max_den <= 400:
            max_terms = min(50, hard_cap)
    else:
        max_terms = min(args.max_terms, hard_cap)

    cfg = SearchConfig(
        max_terms=max_terms,
        max_den=max_den,
        time_limit_s=args.time,
        node_limit=args.nodes,
        range_cap=args.range_cap,
        force_chain=args.force_chain,
        mode=args.mode,
    )

    print(f"[search] mode={cfg.mode}, force_chain={cfg.force_chain}, "
          f"max_terms={cfg.max_terms}, time={cfg.time_limit_s}s, nodes={cfg.node_limit}")

    best, diag = bounded_best_egyptian_div_chain(r, cfg)

    if best is None:
        print("[result] No exact Egyptian expansion found under max_den constraint within current search limits.")
        print("[diag]", {k: diag[k] for k in ("nodes", "time_s", "best_score")})
        return

    ok = verify(r, best)
    run = longest_divisible_run(best)
    edges = count_divisible_edges(best)

    prefix = "-" if sign < 0 else ""
    if k > 0:
        print(f"[result] {prefix}{k} + ({r.numerator}/{r.denominator}) = {prefix}{k} + {format_egyptian(best)}")
    else:
        print(f"[result] {prefix}{r.numerator}/{r.denominator} = {prefix}({format_egyptian(best)})")

    print("[result] denominators:", best)
    print("[result] max_den_used:", max(best))
    print("[result] longest_divisible_run:", run)
    print("[result] divisible_edges:", edges, f"(max possible is {len(best)-1})")
    print("[result] verified:", ok)
    print("[diag]", {k: diag[k] for k in ("nodes", "time_s", "best_score")})


if __name__ == "__main__":
    main()