"""
lineage_survival_model.py

Companion code for: Choueiri A. (2026). "The Jesus DNA Hypothesis: A Quantitative
Reassessment of Patrilineal Survival and Carrier Numbers for a 1st-Century Judean Founder."

Models survival and size of a single founder's male line as a Galton-Watson
branching process (each man has an independent random number of sons).

  * Survival probabilities are computed exactly by iterating the offspring
    probability generating function.
  * Lineage sizes are estimated by Monte Carlo simulation (2,000,000 founders
    per scenario, fixed seed).
  * Offspring distributions: Poisson (variance = mean) and negative binomial
    with shape k (variance = mean + mean^2 / k).
  * Population growth over the period is applied as a constant per-generation
    mean number of sons m = G ** (1 / n).

Also prints the expected number of Y-chromosome and mtDNA point mutations along
one lineage over 1,996 years (30 AD to 2026).

Requirements: Python 3.9+, NumPy.
Usage:        python3 lineage_survival_model.py > results.json
"""

import json

import numpy as np

SEED = 20261111
N_FOUNDERS = 2_000_000
REPORT_GENERATIONS = [1, 10, 20, 50, 65, 80]

# (label, generations, growth factor over the period, NB shape k or None for Poisson)
SCENARIOS = [
    ("constant_poisson_80gen", 80, 1.0, None),
    ("constant_poisson", 65, 1.0, None),
    ("constant_overdispersed", 65, 1.0, 1.0),
    ("growth10_poisson", 65, 10.0, None),
    ("growth10_overdispersed", 65, 10.0, 1.0),
    ("growth30_poisson", 65, 30.0, None),
]

# Mutation-rate inputs
YEARS = 1996                    # 30 AD to 2026
Y_LENGTH_BP = 62.46e6           # T2T-Y (Rhie et al. 2023)
Y_RATE_PER_SITE_YEAR = 8.71e-10 # Helgason et al. 2015
Y_WELL_MEASURED_BP = 15.2e6     # region measured by Helgason et al. 2015
MT_LENGTH_BP = 16569
MT_RATE_PER_SITE_YEAR = 1.665e-8  # Soares et al. 2009, whole mitogenome


def survival_curve(n, m, k=None):
    """Exact P(line alive) after each of n generations via PGF iteration."""
    s, out = 0.0, []
    for _ in range(n):
        if k is None:
            s = np.exp(m * (s - 1.0))
        else:
            s = (1.0 + m * (1.0 - s) / k) ** (-k)
        out.append(1.0 - s)
    return out


def simulate_sizes(n, m, k, rng, n_founders=N_FOUNDERS):
    """Number of living male-line males after n generations, per founder."""
    z = np.ones(n_founders, dtype=np.int64)
    for _ in range(n):
        alive = z > 0
        if k is None:
            # sum of z iid Poisson(m) = Poisson(m * z)
            z[alive] = rng.poisson(m * z[alive])
        else:
            # sum of z iid NB(k, p) = NB(k * z, p), p = k / (k + m)
            z[alive] = rng.negative_binomial(k * z[alive], k / (k + m))
    return z


def main():
    rng = np.random.default_rng(SEED)
    results = {"seed": SEED, "founders_per_scenario": N_FOUNDERS, "scenarios": {}}
    for label, n, growth, k in SCENARIOS:
        m = growth ** (1.0 / n)
        surv = survival_curve(n, m, k)
        z = simulate_sizes(n, m, k, rng)
        alive = z[z > 0]
        results["scenarios"][label] = {
            "generations": n,
            "growth_factor": growth,
            "mean_sons_per_man": round(m, 5),
            "nb_shape_k": k,
            "survival_by_generation": {g: round(surv[g - 1], 4) for g in REPORT_GENERATIONS if g <= n},
            "p_survive_final": round(surv[-1], 4),
            "mean_size_all_founders": round(float(z.mean()), 3),
            "surviving_lines": {
                "mean": round(float(alive.mean()), 1),
                "median": float(np.median(alive)),
                "p99": float(np.percentile(alive, 99)),
                "max": int(alive.max()),
            },
            "share_founders_ge_1000": float((z >= 1000).mean()),
            "share_founders_ge_200000": float((z >= 200_000).mean()),
        }

    # Full survival curves for Figure 1 (exact; no random draws, so results above are unchanged)
    curves = {}
    for label, n_curve, growth, k, m_gens in [
        ("constant_poisson", 80, 1.0, None, 65),
        ("growth10_poisson", 80, 10.0, None, 65),
        ("constant_overdispersed", 80, 1.0, 1.0, 65),
    ]:
        m = growth ** (1.0 / m_gens)
        curves[label] = [1.0] + [round(p, 5) for p in survival_curve(n_curve, m, k)]
    results["survival_curves_gen0_to_80"] = curves

    y_full = Y_LENGTH_BP * Y_RATE_PER_SITE_YEAR * YEARS
    results["mutations"] = {
        "years": YEARS,
        "y_expected_full_chromosome": round(y_full, 1),
        "y_expected_well_measured_region": round(Y_WELL_MEASURED_BP * Y_RATE_PER_SITE_YEAR * YEARS, 1),
        "y_sequence_identity": 1 - y_full / Y_LENGTH_BP,
        "mt_expected": round(MT_LENGTH_BP * MT_RATE_PER_SITE_YEAR * YEARS, 2),
        "mt_to_y_per_site_rate_ratio": round(MT_RATE_PER_SITE_YEAR / Y_RATE_PER_SITE_YEAR, 1),
        "generations_at_31y": round(YEARS / 31, 1),
        "generations_at_25y": round(YEARS / 25, 1),
    }
    print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
