"""The 40 Sentinel-2 regions shared with the healpix-resample paper.

The region list, their coordinates and the exact Sentinel-2 product each patch
came from are read from healpix-resample's frozen `data_manifest.csv`, so this
study cannot drift away from the published dataset: same sites, same products,
same radiometric convention.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

import paper_config as cfg

REGION_SUBDIR = "multi_patch_sentinel2"
_KV = re.compile(r"(\w+)=([^;]+)")


def _parse_source_identifier(text: str) -> dict:
    return {k: v for k, v in _KV.findall(str(text))}


def load_sites(manifest: Path | None = None) -> pd.DataFrame:
    """Return one row per region: scene id, class, coordinates, pinned product.

    Columns: scene, scene_class, region_id, lat, lon, product_id, catalog,
             archive_doi, resample_relative_path.
    """
    root = cfg.resample_root()
    manifest = Path(manifest) if manifest else root / "data_manifest.csv"
    df = pd.read_csv(manifest)

    sel = df[
        (df["category"] == "primary_input")
        & df["relative_path"].astype(str).str.contains(REGION_SUBDIR)
    ].copy()
    if sel.empty:
        raise RuntimeError(
            f"No region patches found in {manifest}. Expected rows whose "
            f"relative_path contains {REGION_SUBDIR!r}."
        )

    rows = []
    for _, r in sel.iterrows():
        stem = Path(str(r["relative_path"])).name
        scene = stem[:-len("_data.zarr")] if stem.endswith("_data.zarr") else stem
        parts = scene.split("__")
        if len(parts) != 3:
            raise RuntimeError(f"Unexpected region scene id {scene!r} in the manifest.")
        kv = _parse_source_identifier(r["source_identifier"])
        if "lat" not in kv or "lon" not in kv:
            raise RuntimeError(
                f"{scene}: manifest source_identifier carries no lat/lon "
                f"({r['source_identifier']!r})."
            )
        source = str(r["source"])
        catalog = "eopf" if "EOPF" in source.upper() else "earth-search"
        rows.append(dict(
            scene=scene,
            scene_class=parts[1],
            region_id=parts[2],
            lat=float(kv["lat"]),
            lon=float(kv["lon"]),
            lattice=kv.get("lattice", ""),
            product_id=kv.get("product", ""),
            catalog=catalog,
            archive_doi=str(r.get("archive_doi", "") or ""),
            resample_relative_path=str(r["relative_path"]),
        ))

    out = pd.DataFrame(rows).sort_values(["scene_class", "region_id"]).reset_index(drop=True)
    out = _merge_resample_status(out, root)
    out["abs_lat"] = out["lat"].abs()
    out["hpx_regime"] = out["abs_lat"].apply(
        lambda a: "equatorial" if a < cfg.EQUATORIAL_LAT_DEG else "cap"
    )
    return out


def _merge_resample_status(df: pd.DataFrame, root: Path) -> pd.DataFrame:
    """Annotate each region with how healpix-resample's own acquisition went.

    Four regions (black_forest, finland_boreal, lake_constance, lake_geneva)
    are recorded there as `failed`: earth-search served their `red` asset from
    the requester-pays JP2 archive, so the companion study fell back to the
    EOPF Zarr Sample Service for them. Knowing this up front is what tells the
    difference between "our acquisition is broken" and "this site is not
    served as a public COG".
    """
    df["resample_status"] = ""
    products = root / "notebooks" / "tables" / "real_groundtruth_multiregion_products.csv"
    if not products.exists():
        return df
    try:
        prod = pd.read_csv(products)[["scene", "status", "product_id"]]
    except Exception:                                     # noqa: BLE001
        return df
    prod = prod.rename(columns={"status": "resample_status",
                                "product_id": "resample_product_id"})
    df = df.drop(columns=["resample_status"]).merge(prod, on="scene", how="left")
    df["resample_status"] = df["resample_status"].fillna("")
    return df


def select(df: pd.DataFrame, scenes=None, classes=None, limit=None) -> pd.DataFrame:
    if scenes:
        wanted = set(scenes)
        df = df[df.scene.isin(wanted) | df.region_id.isin(wanted)]
    if classes:
        df = df[df.scene_class.isin(set(classes))]
    if limit:
        df = df.head(int(limit))
    return df.reset_index(drop=True)


def summary(df: pd.DataFrame) -> str:
    per_class = df.scene_class.value_counts().to_dict()
    per_regime = df.hpx_regime.value_counts().to_dict()
    txt = (f"{len(df)} regions | classes {per_class} | HEALPix regime {per_regime} | "
           f"lat {df.lat.min():.1f}..{df.lat.max():.1f} deg | "
           f"pinned products {int((df.product_id != '').sum())}/{len(df)}")
    if "resample_status" in df:
        failed = df[df.resample_status == "failed"]
        if len(failed):
            txt += (f"\n  note: healpix-resample's own earth-search acquisition "
                    f"failed for {len(failed)} of these and fell back to EOPF: "
                    + ", ".join(failed.region_id))
    return txt


if __name__ == "__main__":
    sites = load_sites()
    print(summary(sites))
    print(sites[["scene", "lat", "lon", "hpx_regime", "catalog", "product_id"]]
          .to_string(index=False))
