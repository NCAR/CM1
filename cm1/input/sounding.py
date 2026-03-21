import io
import logging
import os
from pathlib import Path
import textwrap
from typing import Optional, TextIO, Tuple, Callable

from matplotlib.figure import Figure
import matplotlib.transforms as transforms
import metpy.calc as mpcalc
import metpy.constants
from metpy.interpolate import interpolate_1d
from metpy.plots import Hodograph, SkewT
import numpy as np
import pandas as pd
import xarray as xr
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pint import Quantity

import cm1.input.era5
from cm1.utils import get_ofile, parse_args

# Assuming this script is located in a subdirectory of the repository
repo_base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
soundings_path = os.path.join(repo_base_path, "soundings")


class Sounding:
    """
    A composition-based representation of an atmospheric sounding,
    tailored for CM1 model input/output.

    Wraps an xarray.Dataset rather than subclassing it, avoiding the fragile
    xarray subclass contract (e.g. __init__ bypass by .sel/.map, missing
    'shape' on scalar Quantities, etc.).

    Attribute access is delegated to the underlying dataset, so
    ``sounding.T``, ``sounding.sel(...)``, etc. all work as expected.
    """

    def __init__(self, ds: xr.Dataset, **attrs):
        """
        Parameters
        ----------
        ds : xr.Dataset
            Dataset to wrap.  Will be quantified (metpy) and sorted by
            descending pressure (surface → top) if not already.
        **attrs
            Extra metadata to store in ``ds.attrs`` (e.g. ``case``,
            ``source_file``).
        """
        if not isinstance(ds, xr.Dataset):
            raise TypeError(f"Expected xr.Dataset, got {type(ds)}")

        ds = ds.metpy.quantify()
        ds = ds.sortby("P", ascending=False)
        ds = ds.load()  # always a single profile. Just load eagerly.
        ds.attrs.update(attrs)
        self._ds = ds

    # ------------------------------------------------------------------
    # Attribute delegation
    # ------------------------------------------------------------------

    def __getattr__(self, name: str):
        # Avoid infinite recursion for private / dunder names
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self._ds, name)

    def __getitem__(self, key):
        return self._ds[key]

    def __setitem__(self, key, value):
        self._ds[key] = value

    def __contains__(self, key):
        return key in self._ds

    def __repr__(self) -> str:
        return self._repr_str()

    def _repr_html_(self) -> str:
        """Rich HTML summary for Jupyter notebooks."""
        lines = self._summary_lines()
        rows = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in lines)
        return (
            "<table><thead><tr><th>Sounding</th><th></th></tr></thead>"
            f"<tbody>{rows}</tbody></table>"
        )

    def _repr_str(self) -> str:
        lines = self._summary_lines()
        width = max(len(k) for k, _ in lines)
        return "Sounding\n" + "\n".join(f"  {k:<{width}} : {v}" for k, v in lines)

    def _summary_lines(self) -> list[tuple[str, str]]:
        """Key/value pairs shown in __repr__ and _repr_html_."""
        ds = self._ds
        lines = []

        if "time" in ds:
            lines.append(
                ("time", pd.Timestamp(ds.time.values).strftime("%Y-%m-%d %H:%M"))
            )
        if "latitude" in ds:
            lines.append(
                ("lat/lon", f"{ds.latitude.item():.3f} / {ds.longitude.item():.3f}")
            )
        if "SP" in ds:
            lines.append(("sfc pressure", f"{ds.SP.item().m_as('hPa'):.1f} hPa"))
        if "P" in ds:
            p = ds.P.data
            lines.append(
                (
                    "levels",
                    f"{p.size}  ({float(p.max().m_as('hPa')):.0f} – {float(p.min().m_as('hPa')):.0f} hPa)",
                )
            )
        if "case" in ds.attrs:
            lines.append(("case", ds.attrs["case"]))
        if "source_file" in ds.attrs:
            lines.append(("source", ds.attrs["source_file"]))

        # CAPE/CIN
        cape, cin = _sharp_mucape_cin(ds.P.data, ds.T.data, ds.Q.data, ds.Z.data)
        lines.append(("mucape", f"{cape:.0f} J/kg"))
        lines.append(("mucin", f"{cin:.0f} J/kg"))

        return lines

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def dataset(self) -> xr.Dataset:
        """Return the underlying xr.Dataset."""
        return self._ds

    def plot(self, fig=None, subplot=None, **kwargs) -> "SkewT":
        """Plot the sounding as a Skew-T diagram."""
        return skewt(self._ds, fig=fig, subplot=subplot, **kwargs)

    def to_txt(self) -> str:
        """
        Convert the sounding to a CM1 input_sounding text string.
        """
        ds = self._ds
        sfc = ds.level.sel(level=ds.P.idxmax())
        sfc_pres = ds.SP if "SP" in ds else ds.P.sel(level=sfc)

        theta = mpcalc.potential_temperature(ds.P, ds.T).metpy.convert_units("K")

        sfc_theta_K = (
            ds["surface_potential_temperature"].metpy.convert_units("K")
            if "surface_potential_temperature" in ds
            else theta.sel(level=sfc)
        )

        if "qv" in ds:
            qv = ds["qv"].metpy.convert_units("g/kg")
        elif "Q" in ds:
            qv = mpcalc.mixing_ratio_from_specific_humidity(
                ds["Q"]
            ).metpy.convert_units("g/kg")
        else:
            raise KeyError("Sounding must contain 'qv' or 'Q'.")

        sfc_qv_gkg = (
            ds["surface_mixing_ratio"].metpy.convert_units("g/kg")
            if "surface_mixing_ratio" in ds
            else qv.sel(level=sfc)
        )

        header = (
            f"{sfc_pres.item().m_as('hPa'):.2f} "
            f"{sfc_theta_K.item().m_as('K'):.2f} "
            f"{sfc_qv_gkg.item().m_as('g/kg'):.2f}\n"
        )

        # Build a clean local dataset for CSV serialisation
        out_ds = xr.Dataset(
            {
                "Z": ds["Z"],
                "theta": theta,
                "qv": qv,
                "U": ds["U"],
                "V": ds["V"],
            }
        )
        columns = ["Z", "theta", "qv", "U", "V"]
        body = (
            out_ds.sortby("Z")
            .to_dataframe()
            .to_csv(
                sep=" ", columns=columns, header=False, index=False, float_format="%.2f"
            )
        )
        return header + body


def _sharp_mucape_cin(p, t, q, z):
    """Fast CAPE/CIN via SHARPlib. Inputs are plain float32 numpy arrays in MKS."""
    from nwsspc.sharp.calc import layer, parcel, thermo

    p_pa = np.ascontiguousarray(p.m_as("Pa"), dtype=np.float32)
    t_k = np.ascontiguousarray(t.m_as("K"), dtype=np.float32)
    q = np.ascontiguousarray(q.m_as("kg/kg"), dtype=np.float32)
    hght = np.ascontiguousarray(z.m_as("m"), dtype=np.float32)

    mixr = (q / (1.0 - q)).astype(np.float32)
    td_k = thermo.temperature_at_mixratio(mixr, p_pa)
    td_k = np.minimum(td_k, t_k)  # clamp supersaturation
    mixr = thermo.mixratio(p_pa, td_k).astype(np.float32)
    vtmp = thermo.virtual_temperature(t_k, mixr)
    lft = parcel.lifter_cm1()
    mu_lyr = layer.PressureLayer(p_pa[0], p_pa[0] - 30000.0)
    mu_pcl = parcel.Parcel.most_unstable_parcel(
        mu_lyr, lft, p_pa, hght, t_k, vtmp, td_k
    )
    pcl_vtmp = mu_pcl.lift_parcel(lft, p_pa)
    buoy = thermo.buoyancy(pcl_vtmp, vtmp)
    mu_pcl.find_lfc_el(p_pa, hght, buoy)
    return mu_pcl.cape_cinh(p_pa, hght, buoy)


# ---------------------------------------------------------------------------
# Pressure integration (vectorised)
# ---------------------------------------------------------------------------


def _integrate_pressure(
    Z: Quantity,
    theta: Quantity,
    qv: Quantity,
    sfc_pres: Quantity,
    sfc_theta: Quantity,
    sfc_mix: Quantity,
) -> Quantity:
    """
    Integrate the hydrostatic equation upward to get pressure at every level.

    Vectorised implementation using ``np.cumsum`` instead of a Python loop,
    matching the physics in George Bryan's getcape.f90.

    Parameters
    ----------
    Z, theta, qv : Pint Quantity arrays, shape (N,), sorted surface→top.
    sfc_pres, sfc_theta, sfc_mix : scalar Pint Quantities for the surface.

    Returns
    -------
    P : Pint Quantity array, shape (N,)
    """
    g = metpy.constants.g
    Cp = metpy.constants.Cp_d
    kappa = metpy.constants.kappa
    P0 = 1000.0 * units.hPa

    theta_v = mpcalc.virtual_temperature(theta, qv)
    theta_v_sfc = mpcalc.virtual_temperature(sfc_theta, sfc_mix)

    # Layer boundaries: [sfc, Z[0], Z[1], ..., Z[N-1]]
    z_bounds = np.concatenate([[0.0 * units.m], Z])
    tv_bounds = np.concatenate([[theta_v_sfc], theta_v])

    # Layer-mean virtual potential temperature and thickness
    tv_mean = 0.5 * (tv_bounds[:-1] + tv_bounds[1:])  # shape (N,)
    dz = np.diff(z_bounds)  # shape (N,)

    # Exner increments per layer, then cumulative sum from surface
    dpi = -(g * dz) / (Cp * tv_mean)
    pi_sfc = (sfc_pres / P0) ** kappa
    pi = pi_sfc + np.cumsum(dpi)

    P = P0 * (pi ** (1.0 / kappa))
    return P.to(sfc_pres.units)


# ---------------------------------------------------------------------------
# CM1 text file parser
# ---------------------------------------------------------------------------


def _parse_cm1_txt_stream(
    stream: TextIO,
    case_name_hint: str,
    source_file: str = None,
) -> Sounding:
    """Parse a CM1 input_sounding text stream into a :class:`Sounding`."""
    header = stream.readline().strip()
    surface_pressure, surface_theta, surface_mixing_ratio = map(float, header.split())

    column_names = ["Z", "theta", "qv", "U", "V"]
    df = pd.read_csv(stream, sep=r"\s+", names=column_names, engine="python")
    df = df.rename_axis("level")
    ds = df.to_xarray()
    ds["level"] = ds.level.values

    ds["qv"] = xr.Variable("level", ds["qv"].data * units("g/kg"))
    ds["theta"] = xr.Variable("level", ds["theta"].data * units.K)
    ds["Z"] = xr.Variable("level", ds["Z"].data * units.m)
    ds["U"] = xr.Variable("level", ds["U"].data * units("m/s"))
    ds["V"] = xr.Variable("level", ds["V"].data * units("m/s"))
    ds["Q"] = mpcalc.specific_humidity_from_mixing_ratio(ds["qv"])

    ds["SP"] = np.array(surface_pressure) * units.hPa

    P = _integrate_pressure(
        ds.Z.data,
        ds.theta.data,
        ds.qv.data,
        ds.SP.data,
        surface_theta * units.K,
        surface_mixing_ratio * (units.g / units.kg),
    )

    ds["P"] = ("level", P)
    ds["T"] = mpcalc.temperature_from_potential_temperature(ds["P"], ds["theta"])
    ds["Tv"] = mpcalc.virtual_temperature(ds.T, ds.qv)

    ds["surface_potential_temperature"] = np.array(surface_theta) * units.K
    ds["surface_mixing_ratio"] = np.array(surface_mixing_ratio) * (units.g / units.kg)
    ds["surface_geopotential_height"] = np.array(0.0) * units.m

    case_name = case_name_hint.replace("input_sounding_", "")
    attrs = {"case": case_name}
    if source_file:
        attrs["source_file"] = source_file

    return Sounding(ds, **attrs)


def from_txt(source: str | os.PathLike) -> Sounding:
    """
    Create a :class:`Sounding` by parsing a CM1 input_sounding file or string.

    Parameters
    ----------
    source : str or path-like
        Either a file path or a multi-line string containing the sounding data.
    """
    if isinstance(source, str) and "\n" in source:
        return _parse_cm1_txt_stream(io.StringIO(source), case_name_hint="from_string")

    path = Path(source)
    if not path.is_file():
        raise FileNotFoundError(f"Sounding file not found: {path}")

    with path.open("r") as f:
        return _parse_cm1_txt_stream(
            f, case_name_hint=path.stem, source_file=str(path.absolute())
        )


def get_case(case: str) -> Sounding:
    """Retrieve a predefined sounding case from the repository ``soundings/`` directory."""
    file_path = os.path.join(soundings_path, f"input_sounding_{case}")
    return from_txt(file_path)


# ---------------------------------------------------------------------------
# ERA5 retrieval helpers
# ---------------------------------------------------------------------------


def ensure_quantity(val: float | Quantity, default_unit) -> Quantity:
    """
    Ensure *val* is a Pint Quantity.

    If a plain float is passed, wraps it in *default_unit* and logs at DEBUG
    level (missing units is normal caller behaviour, not a warning).
    """
    if not hasattr(val, "units"):
        logging.debug(f"Input {val!r} has no units; assuming {default_unit}.")
        return val * default_unit
    return val


def _fetch_era5_sounding(
    fetch_func: Callable,
    time: pd.Timestamp | np.datetime64 | str,
    lat: Quantity,
    lon: Quantity,
    **kwargs,
) -> Sounding:
    """
    Shared ERA5 retrieval logic: unit validation, longitude normalization,
    spatial selection, and Sounding construction.
    """
    lat = ensure_quantity(lat, units.degreeN)
    lon = ensure_quantity(lon, units.degreeE)

    valid_time = pd.Timestamp(time)
    ds = fetch_func(valid_time, **kwargs)

    lon = lon % (360 * units.degreeE)
    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    ds = ds.sortby("P", ascending=False)

    return Sounding(ds)


def era5_aws(time, lat, lon, **kwargs) -> Sounding:
    """Retrieve an ERA5 sounding for *time* / *lat* / *lon* from AWS S3."""
    return _fetch_era5_sounding(cm1.input.era5.aws, time, lat, lon, **kwargs)


def era5_model_level(time, lat, lon, **kwargs) -> Sounding:
    """Retrieve an ERA5 model-level sounding from GLADE campaign storage."""
    return _fetch_era5_sounding(cm1.input.era5.model_level, time, lat, lon, **kwargs)


def era5_pressure_level(time, lat, lon, **kwargs) -> Sounding:
    """Retrieve an ERA5 pressure-level sounding from GLADE campaign storage."""
    return _fetch_era5_sounding(cm1.input.era5.pressure_level, time, lat, lon, **kwargs)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    args = parse_args()
    valid_time = pd.to_datetime(args.time)
    ofile = get_ofile(args)

    if os.path.exists(ofile):
        logging.warning(f"Reading from cache: {ofile}")
        sounding = Sounding(xr.open_dataset(ofile))
    else:
        if os.path.exists("/glade/campaign"):
            sounding = era5_model_level(valid_time, args.lat, args.lon)
        else:
            logging.warning(
                "No campaign storage — fetching pressure-level data from AWS."
            )
            sounding = era5_aws(valid_time, args.lat, args.lon)

        logging.warning(f"Caching data to: {ofile}")
        sounding.dataset.metpy.dequantify().to_netcdf(ofile)

    print(sounding.to_txt())


# ---------------------------------------------------------------------------
# Skew-T plotting
# ---------------------------------------------------------------------------


def skewt(
    ds: xr.Dataset,
    fig: Optional[Figure] = None,
    subplot: Optional[Tuple[int, int, int]] = None,
    rotation: int = 35,
    ptop: Quantity = 100 * units.hPa,
    xlim: Optional[Tuple[float, float]] = (-40, 55),
) -> SkewT:
    """
    Generate a Skew-T / hodograph diagram from a sounding dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Sounding data.  Required variables: ``T``, ``Q``, ``U``, ``V``,
        ``Z``, ``P``, ``SP``, ``surface_geopotential_height``.
    fig : Figure, optional
        Existing Matplotlib figure.  Created if *None*.
    subplot : (nrows, ncols, index), optional
        Subplot spec; passed directly to :class:`~metpy.plots.SkewT`.
    rotation : int
        Skew-T rotation angle (degrees).
    ptop : Quantity
        Upper pressure limit for the plot.
    xlim : (float, float), optional
        Temperature axis limits (°C).

    Returns
    -------
    skew : SkewT
    """
    assert "surface_geopotential_height" in ds, (
        "skewt needs surface_geopotential_height"
    )
    assert "SP" in ds, "skewt needs surface pressure SP"

    if fig is not None and (not subplot or subplot == (1, 1, 1)):
        fig.clear()

    logging.info("sorting by descending pressure")
    ds = ds.sortby("P", ascending=False)
    old_nlevel = ds.level.size

    mask = (ds.P >= 10 * units.hPa) & (ds.P < ds.SP)
    ds = ds.sel(level=mask)
    logging.info(f"After masking: {ds.level.size}/{old_nlevel} levels remain")

    p = ds.P.data
    T = ds.T.data
    Td = mpcalc.dewpoint_from_specific_humidity(p, ds.Q.data)

    # Clip supersaturated levels using the unit-safe xarray .where()
    exceedance = (Td - T) > 0 * units.delta_degC
    num_affected = int(exceedance.sum())
    if num_affected > 0:
        max_excess = (Td - T).max()
        logging.warning(
            f"Supersaturation in {num_affected} levels; max Td–T = {max_excess:.2f}. Clipping."
        )
        Td = np.minimum(Td, T)
    else:
        logging.info("All levels Td ≤ T.")

    u = ds.U.data
    v = ds.V.data
    height = ds.Z.data
    sfc_hgt = ds.surface_geopotential_height.data
    agl = (height - sfc_hgt).to("m")
    sp = ds.SP.data

    barb_increments = {"flag": 25, "full": 5, "half": 2.5}
    plot_barbs_units = "m/s"

    # --- Skew-T background ---
    skew = SkewT(fig, subplot=subplot, rotation=rotation)
    skew.plot_dry_adiabats(lw=0.75, alpha=0.5)
    skew.plot_moist_adiabats(lw=0.75, alpha=0.25)

    mixing_ratio = np.array([0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032])
    lines = skew.plot_mixing_lines(mixing_ratio=mixing_ratio, alpha=0.5)
    for i, seg in enumerate(lines.get_segments()):
        top_idx = np.argmin(seg[:, 1])
        x_pos, y_pos = seg[top_idx]
        skew.ax.text(
            x_pos,
            y_pos,
            f"{mixing_ratio[i] * 1000:.0f}",
            fontsize="xx-small",
            color="green",
            ha="center",
            va="bottom",
            alpha=0.5,
            clip_on=True,
        )

    skew.ax.axvline(0, color="c", linestyle="--", linewidth=1.5, alpha=0.5)

    # --- Temperature traces ---
    if "Tv" in ds:
        Tv = ds["Tv"].data
    else:
        logging.info("Deriving Tv from p, T, Td")
        Tv = mpcalc.virtual_temperature_from_dewpoint(p, T, Td)

    skew.plot(p, T, "r")
    skew.plot(p, Tv, "r", lw=0.5, linestyle="dashed")
    skew.plot(p, Td, "g")

    # --- Parcel path ---
    p_mu, t_mu, td_mu, idx_mu = mpcalc.most_unstable_parcel(p, T, Td)
    logging.info(
        f"MU parcel: {p_mu:~.1f}  {t_mu.to('degC'):~.2f}  {td_mu:~.2f}  i={idx_mu}"
    )

    lcl_p, lcl_temperature = mpcalc.lcl(p_mu, t_mu, td_mu)
    logging.info(f"LCL: {lcl_p:~.1f}  {lcl_temperature:~.2f}")

    trans = transforms.blended_transform_factory(skew.ax.transAxes, skew.ax.transData)
    hpos = 0.82
    skew.ax.plot([hpos, hpos + 0.03], 2 * [lcl_p.m_as("hPa")], transform=trans)
    skew.ax.text(
        hpos + 0.035,
        lcl_p,
        f" LCL {lcl_p.to('hPa').item():~.0f}",
        transform=trans,
        horizontalalignment="left",
        verticalalignment="center",
        fontsize="x-small",
    )

    p_path = p[idx_mu:]
    Tv_env_path = Tv[idx_mu:]
    T_parcel_path = mpcalc.parcel_profile(p_path, t_mu, td_mu)

    e_mu = mpcalc.saturation_vapor_pressure(td_mu)
    w_parcel_start = mpcalc.mixing_ratio(e_mu, p_mu)
    w_parcel_path = mpcalc.saturation_mixing_ratio(p_path, T_parcel_path)
    w_parcel_path[p_path >= lcl_p] = w_parcel_start
    Tv_parcel_path = mpcalc.virtual_temperature(T_parcel_path, w_parcel_path)

    skew.plot(p_path, T_parcel_path, "k", ls="--")
    skew.plot(p_path, Tv_parcel_path, "k", lw=0.5, ls="--")
    skew.shade_cin(p_path, Tv_env_path, Tv_parcel_path)
    skew.shade_cape(p_path, Tv_env_path, Tv_parcel_path)

    skew.ax.set_xlim(xlim)
    skew.ax.set_ylim(None, ptop)

    # --- AGL height labels on y-axis ---
    label_hgts = np.array([0, 1, 3, 6, 9, 12, 15]) * units.km
    for label_hgt in label_hgts:
        if label_hgt >= agl.min() and label_hgt <= agl.max():
            (agl2p,) = interpolate_1d(label_hgt, agl, p)
        s = f"{label_hgt:~.0f}"
        if label_hgt == 0 * units.km:
            agl2p = sp.item()
            s = f"SFC {ds.surface_geopotential_height.item():~.0f}"
        skew.ax.plot([0, 0.01], 2 * [agl2p.m_as("hPa")], transform=trans, color="k")
        skew.ax.text(
            0.015,
            agl2p,
            s,
            transform=trans,
            horizontalalignment="left",
            verticalalignment="center",
            fontsize="x-small",
            clip_on=True,
        )

    # --- Build title once; set it once at the end ---
    cape, cin = _sharp_mucape_cin(p, T, ds.Q.data, height)
    title_parts = []
    if "case" in ds.attrs:
        title_parts.append(ds.attrs["case"])
    if "time" in ds:
        title_parts.append(ds.time.dt.strftime("%Y-%m-%d %H:%M:%S").item())
    else:
        logging.info("No 'time' variable in sounding")
    if "longitude" in ds:
        title_parts.append(f"{ds.longitude.item():.3f} E {ds.latitude.item():.3f} N")
    title_parts.append(f"mucape: {cape:.0f} J/kg  mucin: {cin:.0f} J/kg")

    skip_winds = not u.any() and not v.any()
    if not skip_winds:
        logging.info("Computing wind kinematics")
        right_mover, left_mover, wind_mean = mpcalc.bunkers_storm_motion(
            p, u, v, height
        )
        storm_u, storm_v = wind_mean
        srh03_pos, srh03_neg, srh03_tot = mpcalc.storm_relative_helicity(
            height, u, v, 3 * units.km, storm_u=storm_u, storm_v=storm_v
        )
        title_parts.append(
            f"0-3km srh+: {srh03_pos:~.0f}   srh-: {srh03_neg:~.0f}   srh(tot): {srh03_tot:~.0f}"
        )

    full_title_str = " | ".join(title_parts)
    wrapped_title = textwrap.fill(full_title_str, width=75)
    skew.ax.set_title(wrapped_title, fontsize="x-small", loc="center")

    if skip_winds:
        return skew

    # --- Wind barbs ---
    agl_colors = ["red", "red", "lime", "green", "blueviolet", "cyan"]
    indices = np.searchsorted(label_hgts.m_as("m"), agl.m_as("m"), side="right") - 1
    indices = np.clip(indices, 0, len(agl_colors) - 1)
    barbcolor = [agl_colors[i] for i in indices]

    skew.plot_barbs(
        p,
        u,
        v,
        barbcolor=barbcolor,
        length=6,
        plot_units=plot_barbs_units,
        linewidth=0.6,
        xloc=1.05,
        barb_increments=barb_increments,
    )
    skew.ax.text(
        1.07,
        0,
        plot_barbs_units,
        va="top",
        ha="right",
        transform=skew.ax.transAxes,
        fontsize="xx-small",
    )

    # --- Hodograph inset ---
    logging.info("Drawing hodograph")
    ax_hod = inset_axes(skew.ax, "40%", "40%", loc=1)
    h = Hodograph(ax_hod)
    h.add_grid(increment=10, linewidth=0.75)

    u_center = (u.min() + u.max()) / 2
    v_center = (v.min() + v.max()) / 2
    max_range = max(u.max() - u.min(), v.max() - v.min()).to(plot_barbs_units).m
    plot_radius = (max_range / 2) * 1.1

    ax_hod.set_xlim(u_center.m - plot_radius, u_center.m + plot_radius)
    ax_hod.set_ylim(v_center.m - plot_radius, v_center.m + plot_radius)
    ax_hod.set_ylabel("")
    ax_hod.xaxis.label.set_fontsize("xx-small")
    ax_hod.tick_params(axis="both", which="major", labelsize="xx-small")

    for label_hgt in label_hgts:
        ax_hod.text(
            np.interp(label_hgt, agl, u),
            np.interp(label_hgt, agl, v),
            label_hgt.to("km").m,
            fontsize=7,
        )

    h.plot_colormapped(u, v, agl, intervals=label_hgts, colors=agl_colors)
    ax_hod.plot(storm_u, storm_v, "x")

    return skew


if __name__ == "__main__":
    main()
