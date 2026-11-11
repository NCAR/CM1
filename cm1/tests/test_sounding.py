# tests/test_sounding.py
import pytest
from cm1.input.sounding import from_txt

SIMPLE_SOUNDING = """\
1000.00 300.00 14.00
0.0 300.00 14.00 5.0 0.0
500.0 302.00 12.00 6.0 1.0
1000.0 304.00 10.00 7.0 2.0
"""


def test_round_trip():
    """from_txt -> to_txt -> from_txt should reproduce the sounding."""
    original = from_txt(SIMPLE_SOUNDING)
    txt = original.to_txt()
    recovered = from_txt(txt)

    # Check surface header values survive the round-trip
    assert abs(original.SP.item().m_as("hPa") - recovered.SP.item().m_as("hPa")) < 0.01

    # Check profile values at every level
    import numpy as np

    np.testing.assert_allclose(original.Z.values, recovered.Z.values, rtol=1e-4)
    np.testing.assert_allclose(original["U"].values, recovered["U"].values, rtol=1e-4)


def test_from_txt_missing_file():
    with pytest.raises(FileNotFoundError):
        from_txt("/nonexistent/path")


def test_sounding_contains_expected_variables():
    s = from_txt(SIMPLE_SOUNDING)
    for var in ["T", "P", "Z", "U", "V", "Q"]:
        assert var in s, f"Missing variable: {var}"


def test_pressure_decreases_with_height():
    """Pressure should be monotonically decreasing from surface to top."""
    import numpy as np

    s = from_txt(SIMPLE_SOUNDING)
    p = s.P.values
    assert np.all(np.diff(p) < 0), "Pressure is not monotonically decreasing"
