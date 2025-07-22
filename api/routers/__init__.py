from importlib import import_module

# Expose all router subpackages for convenient access,
# e.g. `routers.prompt.router`.
for _name in ("prompt", "geometry", "render", "rcsb"):
    import_module(f"{__name__}.{_name}")

__all__ = ["prompt", "geometry", "render", "rcsb"]
