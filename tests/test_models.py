"""Rate laws: limits, special cases and the refusals that keep fits honest."""

from __future__ import annotations

import numpy as np
import pytest

from kinetics.models import (
    INHIBITION_MODELS,
    competitive,
    hill,
    michaelis_menten,
    mixed_inhibition,
    noncompetitive,
    uncompetitive,
)

VMAX, KM, KI = 10.0, 25.0, 8.0


def test_zero_substrate_gives_zero_rate():
    assert michaelis_menten(0.0, VMAX, KM) == 0.0


def test_khalf_is_the_substrate_giving_half_vmax():
    assert michaelis_menten(KM, VMAX, KM) == pytest.approx(VMAX / 2)


def test_rate_approaches_vmax_but_never_reaches_it():
    assert michaelis_menten(1e9, VMAX, KM) < VMAX
    assert michaelis_menten(1e9, VMAX, KM) > 0.999 * VMAX


def test_rate_is_monotone_in_substrate():
    S = np.array([0.1, 1, 5, 10, 25, 100, 1000])
    v = michaelis_menten(S, VMAX, KM)
    assert np.all(np.diff(v) > 0)


def test_scalar_input_returns_python_float():
    """numpy scalars would replay as ``np.float64(5.0)`` and break json.dumps."""
    out = michaelis_menten(25.0, VMAX, KM)
    assert type(out) is float
    assert f"{out}" == "5.0"


def test_array_input_returns_array():
    out = michaelis_menten([0.0, 25.0], VMAX, KM)
    assert isinstance(out, np.ndarray)
    assert out.shape == (2,)


@pytest.mark.parametrize("km", [0.0, -1.0])
def test_non_positive_km_is_refused(km):
    with pytest.raises(ValueError, match="km must be positive"):
        michaelis_menten(10.0, VMAX, km)


def test_negative_substrate_is_refused():
    with pytest.raises(ValueError, match="cannot be negative"):
        michaelis_menten(-1.0, VMAX, KM)


def test_hill_with_n_equal_one_is_michaelis_menten():
    S = np.array([0.5, 2.0, 25.0, 250.0])
    assert np.allclose(hill(S, VMAX, KM, 1.0), michaelis_menten(S, VMAX, KM))


def test_hill_n_above_one_is_steeper_than_michaelis_menten():
    """Both reach half Vmax at S = K½; the Hill curve gets there faster."""
    assert hill(5.0, VMAX, KM, 3.0) < michaelis_menten(5.0, VMAX, KM)
    assert hill(125.0, VMAX, KM, 3.0) > michaelis_menten(125.0, VMAX, KM)


@pytest.mark.parametrize("func", [competitive, noncompetitive, uncompetitive])
def test_no_inhibitor_reproduces_michaelis_menten(func):
    S = np.linspace(0.1, 200, 25)
    assert np.allclose(func(S, VMAX, KM, 0.0, KI), michaelis_menten(S, VMAX, KM))


def test_competitive_inhibitor_raises_apparent_km_and_spares_vmax():
    """At I = Ki the apparent Km doubles and the saturating rate is unchanged."""
    # Half Vmax is reached at twice the substrate concentration: the apparent Km
    # has doubled from 25 to 50 µM.
    assert competitive(2 * KM, VMAX, KM, KI, KI) == pytest.approx(VMAX / 2)
    assert competitive(KM, VMAX, KM, KI, KI) < VMAX / 2
    saturated_no_inhibitor = competitive(1e9, VMAX, KM, 0.0, KI)
    saturated_inhibited = competitive(1e9, VMAX, KM, KI, KI)
    assert saturated_inhibited == pytest.approx(saturated_no_inhibitor, rel=1e-6)


def test_noncompetitive_inhibitor_halves_vmax_at_i_equals_ki():
    assert noncompetitive(1e9, VMAX, KM, KI, KI) == pytest.approx(VMAX / 2)


def test_uncompetitive_inhibitor_lowers_the_apparent_km():
    """v at S = Km falls below Vmax/2: the curve shifts left, not just down."""
    assert uncompetitive(KM, VMAX, KM, KI, KI) < VMAX / 2


def test_mixed_with_alpha_one_is_noncompetitive():
    S = np.linspace(0.1, 200, 25)
    for I in (0.0, 3.0, 12.0):
        assert np.allclose(
            mixed_inhibition(S, VMAX, KM, I, KI, alpha=1.0), noncompetitive(S, VMAX, KM, I, KI)
        )


def test_mixed_with_large_alpha_approaches_competitive():
    S = np.linspace(0.1, 200, 25)
    assert np.allclose(
        mixed_inhibition(S, VMAX, KM, 8.0, KI, alpha=1e9),
        competitive(S, VMAX, KM, 8.0, KI),
        rtol=1e-6,
    )


@pytest.mark.parametrize("func", [competitive, noncompetitive, uncompetitive])
def test_negative_inhibitor_is_refused(func):
    with pytest.raises(ValueError, match="inhibitor concentration cannot be negative"):
        func(10.0, VMAX, KM, -1.0, KI)


@pytest.mark.parametrize("func", [competitive, noncompetitive, uncompetitive])
def test_non_positive_ki_is_refused(func):
    with pytest.raises(ValueError, match="ki must be positive"):
        func(10.0, VMAX, KM, 1.0, 0.0)


def test_inhibition_models_registry_is_complete_and_callable():
    assert set(INHIBITION_MODELS) == {"competitive", "noncompetitive", "uncompetitive"}
    for name, spec in INHIBITION_MODELS.items():
        assert spec["params"] == ("vmax", "km", "ki"), name
        assert spec["func"](10.0, VMAX, KM, 0.0, KI) > 0


def test_models_do_not_convert_units():
    """Feeding mM instead of µM scales Km by 1000 and leaves Vmax alone."""
    assert michaelis_menten(0.025, VMAX, KM / 1000) == pytest.approx(
        michaelis_menten(25.0, VMAX, KM)
    )
