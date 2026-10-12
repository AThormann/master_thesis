def select_track(ds, track_id):
    """Return the points of a single track."""
    return ds.where(ds.track_id == track_id, drop=True)