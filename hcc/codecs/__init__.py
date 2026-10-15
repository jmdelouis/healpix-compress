from hcc.codecs.packing import pack_array, unpack_array  # noqa: F401


def _register_zarr_codecs() -> None:
    """Explicitly register hcc's Zarr v3 codecs, instead of relying solely
    on the ``[project.entry-points."zarr.codecs"]`` table in
    ``pyproject.toml`` being picked up from installed package metadata.

    That entry-points mechanism only works if the package has been
    (re)installed -- ``pip install -e .`` / ``pixi install`` -- in the
    *exact* environment/kernel currently running, so that Python's
    importlib metadata actually lists the entry points. If it hasn't been
    (a fresh venv, a Jupyter kernel that doesn't match the environment
    ``pixi install`` ran in, an editable install predating a
    pyproject.toml change, ``hcc`` imported from a source tree), Zarr raises
    ``UnknownCodecError("Unknown codec: 'hcc-decompose'")`` (or
    'hcc-grad-dwt', 'hcc-sentinel2', etc.) even though ``import hcc``
    itself succeeds. Registering explicitly at import time removes that
    dependency on install-time metadata entirely. If the entry points *are*
    already correctly registered, re-registering the same class under the
    same name here is a harmless no-op.
    """
    try:
        from zarr.registry import register_codec
    except ImportError:
        # Zarr API without a programmatic registry (older/newer than
        # expected); entry points remain the only registration path.
        return

    from hcc.codecs.dct import HccDctCodec
    from hcc.codecs.decompose import HccDecomposeCodec
    from hcc.codecs.grad_dct import HccGradDctCodec
    from hcc.codecs.grad_dwt import HccGradDwtCodec
    from hcc.codecs.sentinel2 import HccSentinel2Codec
    from hcc.codecs.sentinel2_pyramid import HccSentinel2PyramidCodec

    for _name, _cls in (
        ("hcc-dct", HccDctCodec),
        ("hcc-grad-dct", HccGradDctCodec),
        ("hcc-grad-dwt", HccGradDwtCodec),
        ("hcc-decompose", HccDecomposeCodec),
        (HccSentinel2Codec.name, HccSentinel2Codec),
        (HccSentinel2PyramidCodec.name, HccSentinel2PyramidCodec),
    ):
        try:
            register_codec(_name, _cls)
        except Exception as _e:  # pragma: no cover -- defensive, log and continue
            print(f"[hcc.codecs] could not explicitly register {_name!r}: {_e!r}")


_register_zarr_codecs()
