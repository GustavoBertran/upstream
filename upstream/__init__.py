try:
    from importlib.metadata import PackageNotFoundError, version as _version

    try:
        __version__ = _version("upstream")
    except PackageNotFoundError:
        __version__ = "0.1.0"
except Exception:  # pragma: no cover
    __version__ = "0.1.0"
