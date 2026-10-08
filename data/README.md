
##FILES
site_summary.csv     Elevation and slope statistics computed over the Site07 ROI polygon.
site_coordinates.csv Center point and bounding box of the ROI, in lat/long.
site07_profiles.csv - elevation profiles from the Site07 DEM

##DATA SOURCE
NASA LOLA / LRO site-specific 5 m per pixel DEM and slope products (PGDA).
Site07_final_adj_5mpp_surf.tif (elevation), Site07_final_adj_5mpp_slp.tif (slope).

##UNITS AND REFERENCE
Elevation: meters, relative to the lunar reference sphere (radius 1,737.4 km).
Slope: degrees, measured at 5 m resolution. Slope values depend on the measurement
baseline, so keep that in mind when comparing to lander footprint and tolerance.
Coordinates: decimal degrees, GCS_Moon_2000, longitude as returned by ArcGIS (positive = east).
Near the pole, the bounding-box longitude range is wide even for a ~16 x 16 km site.

##HOW THE STATISTICS WERE COMPUTED
Zonal Statistics as Table (ArcGIS Pro) over the ROI polygon, ignoring NoData.
pct_area_slope_lt5deg / lt10deg: slope raster reclassified to 1 where slope is below
the threshold and 0 elsewhere, then the zonal mean multiplied by 100.
Center point: centroid of the ROI polygon (Feature to Point).
Pixel count 10,227,204 x 25 m^2 = about 255.68 km^2.


#site07_profiles.csv - elevation profiles from the Site07 DEM
  LINE_ID 1 = Transect A (screen-horizontal: left edge midpoint to right edge midpoint of the ROI)
  LINE_ID 3 = Transect B (screen-vertical: top edge midpoint to bottom edge midpoint of the ROI)
  (LINE_ID 2 does not appear: it was a draft line that was deleted and redrawn as ID 3.)
  DISTANCE = meters along the line from its start point; sample spacing about 5 m
  Z = elevation in meters relative to the 1,737.4 km lunar reference sphere
  Each transect is about 16 km long and passes through the site center.
  "Horizontal" and "vertical" refer to the polar stereographic map orientation, not true
  compass directions.
