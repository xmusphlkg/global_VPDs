"""Shared deterministic daily SEIRV transitions for fitted and fixed analyses."""
from __future__ import annotations

import math
import numpy as np
from numba import njit, prange


def _beta(r0: float, infectious_days: float, contacts: np.ndarray) -> tuple[float, float]:
    rho = float(np.max(np.linalg.eigvals(contacts).real))
    if rho <= 0:
        raise ValueError("Contact matrix must have a positive spectral radius")
    return float(r0 * -math.expm1(-1.0 / infectious_days) / rho), rho


@njit(cache=True)
def _simulate_daily_one(population, births, contacts, contact_rho, gaps, r0,
                   latent_days, infectious_days, acquisition_ratio,
                   clinical_fraction, importation, reporting, contact_scale,
                   natural_immunity_years=0.0, source_reduction=0.0):
    """Return annual expected reports and diagnostics; rates are per day.

    Allocate vaccine protection within the noninfected, nonrecovered pool,
    assuming independence of vaccine protection and prior natural immunity.
    """
    p_ei = -math.expm1(-1.0 / latent_days)
    p_ir = -math.expm1(-1.0 / infectious_days)
    p_r_wane = (0.0 if natural_immunity_years <= 0.0 else
                -math.expm1(-1.0 / (365.0 * natural_immunity_years)))
    S = population[0].copy()
    E = np.zeros(16); I = np.zeros(16); R = np.zeros(16); V = np.zeros(16)
    seed = min(max(1.0, 1e-7 * S[1]), S[1])
    I[1] = seed; S[1] -= seed
    annual = np.zeros(len(gaps))
    peak_daily = 0.0; cumulative = 0.0
    infectious_density = np.zeros(16)
    for yi in range(len(gaps)):
        protected = 1.0 - gaps[yi]
        susceptible = np.maximum(1.0 - protected * acquisition_ratio, 1e-12)
        clinical = np.minimum(gaps[yi] / susceptible, 1.0)
        protected_infected = protected * (1.0 - acquisition_ratio) / susceptible
        source = np.maximum(1.0 - source_reduction * protected_infected, 0.0)
        available = S + V
        V = available * (1.0 - susceptible)
        S = available * susceptible
        beta = r0 * p_ir / contact_rho[yi]
        for _ in range(365):
            for age in range(16):
                total = max(S[age] + E[age] + I[age] + R[age] + V[age], 1.0)
                infectious_density[age] = I[age] * source[age] / total
            count = 0.0
            for age in range(16):
                exposure = 0.0
                for donor in range(16):
                    exposure += contacts[yi, age, donor] * infectious_density[donor]
                force = beta * contact_scale[yi] * exposure + importation
                newE = S[age] * min(-math.expm1(-force), 1.0)
                newI = E[age] * p_ei
                newR = I[age] * p_ir
                waned = R[age] * p_r_wane
                count += newI
                annual[yi] += newI * clinical_fraction * clinical[age] * reporting
                S[age] += waned - newE
                E[age] += newE - newI
                I[age] += newI - newR
                R[age] += newR - waned
            cumulative += count; peak_daily = max(peak_daily, count)
            S[0] += births[yi] / 365.0
            for state in (S, E, I, R, V):
                # Descending order keeps all ageing flows based on the same
                # pre-ageing state without allocating a daily flow vector.
                for age in range(15, 0, -1):
                    flow = state[age - 1] / (5.0 * 365.0)
                    state[age - 1] -= flow
                    state[age] += flow
        if yi + 1 < len(population):
            current = S + E + I + R + V
            scale = population[yi + 1] / np.maximum(current, 1e-12)
            S *= scale; E *= scale; I *= scale; R *= scale; V *= scale
    return annual, peak_daily, cumulative, np.vstack((S, E, I, R, V))


def simulate_daily(population, births, contacts, contact_rho, gaps, r0,
                   latent_days, infectious_days, acquisition_ratio,
                   clinical_fraction, importation, reporting, contact_scale,
                   natural_immunity_years=0.0, source_reduction=0.0):
    """Run the shared daily transition kernel for one country."""
    return _simulate_daily_one(population, births, contacts, contact_rho, gaps, r0,
                               latent_days, infectious_days, acquisition_ratio,
                               clinical_fraction, importation, reporting, contact_scale,
                               natural_immunity_years, source_reduction)


@njit(cache=True, parallel=True)
def simulate_daily_panel(population, births, contacts, contact_rho, gaps, r0,
                         latent_days, infectious_days, acquisition_ratio,
                         clinical_fraction, importation, reporting, contact_scale,
                         natural_immunity_years=0.0, source_reduction=0.0):
    """Run identical daily transitions for a stack of countries."""
    n_country = population.shape[0]
    n_year = gaps.shape[1]
    annual = np.zeros((n_country, n_year))
    for ci in prange(n_country):
        result, _, _, _ = _simulate_daily_one(
            population[ci], births[ci], contacts[ci], contact_rho[ci], gaps[ci], r0,
            latent_days, infectious_days, acquisition_ratio, clinical_fraction,
            importation, reporting[ci], contact_scale, natural_immunity_years,
            source_reduction)
        annual[ci] = result
    return annual


def simulate_with_diagnostics(population, births, contacts, targets, r0,
                              latent_days, infectious_days, acquisition_ratio,
                              clinical_fraction=1.0, importation=0.0, reporting=1.0):
    """Fixed-input years through the same daily transition kernel."""
    population = np.asarray(population, dtype=float)
    targets = np.asarray(targets, dtype=float)
    contacts = np.asarray(contacts, dtype=float)
    if population.shape != (16,) or targets.ndim != 2 or targets.shape[1] != 16:
        raise ValueError("Expected sixteen age bands")
    _, rho = _beta(r0, infectious_days, contacts)
    annual, peak, cumulative, states = simulate_daily(
        population[None, :], np.asarray(births, dtype=float),
        np.repeat(contacts[None, :, :], len(targets), axis=0),
        np.full(len(targets), rho), targets, r0, latent_days, infectious_days,
        acquisition_ratio, clinical_fraction, importation, reporting,
        np.ones(len(targets)))
    return annual, {"peak_daily_infections": float(peak),
                    "cumulative_infections": float(cumulative),
                    "final_population_by_age": states.sum(axis=0),
                    "final_states": states}
