



def wrap_lon(lon, to="-180"):
    """Convert longitudes between 0..360 and -180..180 conventions."""
    if to == "-180":
        return ((lon + 180) % 360) - 180
    return lon % 360


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (works on numpy arrays)."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))

def _in_box(lat, lon, box):
    lon_min, lon_max, lat_min, lat_max = box
    return (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)