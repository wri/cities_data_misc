import os, sys, requests, json, pathlib, datetime, boto3
import pandas as pd
import geopandas as gpd
import osmnx as ox
import osmium
import geocube
from botocore.exceptions import ClientError
import r5py
os.chdir('C:/Users/tgwon/wri/cif/cities-cif')
sys.path.append('C:/Users/tgwon/wri/cif/cities-cif')
os.environ['GOOGLE_APPLICATION_USER'] = 'script-service-account@wri-gee.iam.gserviceaccount.com'
os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = 'C:/Users/tgwon/.google/credkey.json'
sys.path = ['C:\\Users\\tgwon\\wri\\indicators', 'C:\\Users\\tgwon\\wri\\cif\\cities-cif', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\python310.zip', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\DLLs', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib\\site-packages', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib\\site-packages\\win32', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib\\site-packages\\win32\\lib', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib\\site-packages\\Pythonwin',]

from city_metrix.layers.layer_geometry import GeoExtent
from city_metrix.layers import OpenStreetMap, OpenStreetMapClass, UrbanExtents, WorldPop

SESSION = boto3.Session(profile_name='CitiesUserPermissionSet-540362055257')
BUCKET = 'wri-cities-indicators'
S3_ISOCHRONE_PREFIX = 'devdata/inputdata'
S3_AMENITYPOINTS_PREFIX = 'devdata/inputdata/amenitypoints'
LOCAL_AMENITYPOINTS_PREFIX = 'C:/Users/tgwon/wri/indicators/amenitypoints'
LOCAL_ISOCHRONE_PREFIX = 'C:/Users/tgwon/wri/indicators/isochrones'
LOCAL_PBF_PREFIX = 'C:/Users/tgwon/wri/indicators/pbf'

AMENITIES = [
    # (name, osm_class, travel_times, travel_distances)
    ('goods_services', [OpenStreetMapClass.COMMERCE], [15, 30], []),
    ('potential_employment', [OpenStreetMapClass.COMMERCE, OpenStreetMapClass.HEALTHCARE_SOCIAL, OpenStreetMapClass.AGRICULTURE, OpenStreetMapClass.GOVERNMENT, OpenStreetMapClass.INDUSTRY, OpenStreetMapClass.TRANSPORTATION_LOGISTICS, OpenStreetMapClass.EDUCATION], [15, 30], []),
    ('healthcare', [OpenStreetMapClass.MEDICAL], [15, 30], [500, 1000]),
    ('education_primary_secondary', [OpenStreetMapClass.PRIMARY_SECONDARY_EDUCATION], [15, 30], [500, 1000, 1250]),
    ('public_transportation', [OpenStreetMapClass.TRANSIT_STOP], [15, 30], [500]),
    ('public_open_space', [OpenStreetMapClass.OPEN_SPACE], [15, 30], [300, 500])
]

TRAVEL_MODES = [r5py.TransportMode.WALK, r5py.TransportMode.BICYCLE]
TRAVEL_SPEEDS = {'walk': 3.6, 'bicycle': 12}

def upload_s3(session, file_name, bucket, object_name):
    s3_client = session.client('s3')
    try:
        response = s3_client.upload_file(file_name, bucket, object_name, ExtraArgs={'ACL': 'public-read'})
    except ClientError as e:
        logging.error(e)
        return False
    return True

def class_name_from_osmclass(osm_class):
    return osm_class.__str__().lower().split('.')[1].replace('_', '-')

def get_data_from_polygon(polygon, osm_class):
    # Set the OSMnx configuration to disable caching
    ox.settings.use_cache = False
    try:
        #osm_feature = ox.features_from_bbox(bbox=(min_lon, min_lat, max_lon, max_lat), tags=self.osm_class.value)
        osm_feature = ox.features_from_polygon(polygon, tags=osm_class.value)
    # When no feature in bbox, return an empty gdf
    except ox._errors.InsufficientResponseError as e:
        osm_feature = gpd.GeoDataFrame(pd.DataFrame(columns=['id', 'geometry']+list(osm_class.value.keys())), geometry='geometry')
        osm_feature.crs = "EPSG:4326"
    osm_feature[f'amenityclass_{class_name_from_osmclass(osm_class)}'] = True

    # Filter by geo_type
    if osm_class == OpenStreetMapClass.ROAD:
        # Filter out Point
        osm_feature = osm_feature[osm_feature.geom_type != 'Point']
    elif osm_class == OpenStreetMapClass.TRANSIT_STOP:
        # Keep Point
        osm_feature = osm_feature[osm_feature.geom_type == 'Point']
    else:
        # Filter out Point and LineString
        osm_feature = osm_feature[osm_feature.geom_type.isin(['Polygon', 'MultiPolygon'])]

    # keep only columns desired to reduce file size
    keep_col = ['id', 'geometry']
    for key in osm_class.value:
        if key in osm_feature.columns:
            keep_col.append(key)
    # keep 'lanes' for 'highway'
    if 'highway' in keep_col and 'lanes' in osm_feature.columns:
        keep_col.append('lanes')
    osm_feature = osm_feature.reset_index()[keep_col]
    osm_feature['id'] = list(range(len(osm_feature)))
    return osm_feature

def get_perimeter_points(city_gdf, osm_class):
    openspace_polys = get_data_from_polygon(city_gdf.dissolve().geometry[0], osm_class)
    osm_class = OpenStreetMapClass.ROAD
    road_lines = get_data_from_polygon(city_gdf.dissolve().geometry[0], osm_class)
    road_lines = road_lines.loc[
        (road_lines.highway != 'motorway') & 
        (road_lines.highway != 'motorway_link') &
        (road_lines.highway != 'trunk') &
        (road_lines.highway != 'trunk_link') &
        (road_lines.highway != 'primary') &
        (road_lines.highway != 'primary_link') &
        (road_lines.highway != 'secondary_link') &
        (road_lines.highway != 'tertiary_link') &
        (road_lines.highway != 'passing_place') &
        (road_lines.highway != 'busway')
    ]
    bbox = GeoExtent(city_gdf.total_bounds)
    utm_epsg = bbox.as_utm_bbox().epsg_code
    buffered_openspace_uu = openspace_polys.to_crs(f'EPSG:{utm_epsg}').buffer(10).to_crs('EPSG:4326').union_all()
    roads_filtered = road_lines.intersects(buffered_openspace_uu)
    road_uu = road_lines.loc[roads_filtered].union_all()
    perimeter_pts = []
    print(len(openspace_polys), end=': ', flush=True)
    for rownum in list(range(len(openspace_polys))):
        if rownum % 1000 == 0:
            print(rownum, end=' ', flush=True)
        multi = openspace_polys.iloc[[rownum]].to_crs(f'EPSG:{utm_epsg}').buffer(5).to_crs('EPSG:4326').boundary.intersection(road_uu)[rownum]
        if multi.geom_type=='MultiPoint':
            perimeter_pts += multi.geoms
        else:
            perimeter_pts.append(multi)
    res = gpd.GeoDataFrame({'id': range(len(perimeter_pts)), 'geometry': perimeter_pts}).set_crs('EPSG:4326').set_geometry('geometry')
    adminboundary_gdf = city_gdf.dissolve()
    a_points_filtered = res.loc[res.within(adminboundary_gdf.geometry[0])]
    a_points_filtered_pointsonly = a_points_filtered.loc[a_points_filtered.geometry.geom_type=='Point']
    a_points_filtered_pointsonly[f'amenityclass_{class_name_from_osmclass(osm_class)}'] = True
    return a_points_filtered_pointsonly.reset_index()

def get_amenities_points(city_gdf, osm_classes):
    bbox = GeoExtent(city_gdf.total_bounds)
    utm_epsg = bbox.as_utm_bbox().epsg_code
    polygon = city_gdf.dissolve().geometry[0]
    osm_results = []
    columns_to_keep = ['geometry']
    ox.settings.use_cache = False  # Set the OSMnx configuration to disable caching
    for osm_class in osm_classes:
        class_name = class_name_from_osmclass(osm_class)
        columns_to_keep.append(f'amenityclass_{class_name}')
        try:
            res = ox.features_from_polygon(polygon, tags=osm_class.value)
        except ox._errors.InsufficientResponseError as e:
            res = gpd.GeoDataFrame(pd.DataFrame(columns=['osmid', 'geometry']+list(merged_osm_dicts.keys())), geometry='geometry')
            res.crs = "EPSG:4326"
        for osm_class_all in osm_classes:
            class_name_all = class_name_from_osmclass(osm_class_all)
            res[f'amenityclass_{class_name}'] = class_name == class_name_all

        osm_results.append(res)
    osm_feature = pd.concat(osm_results, axis=0)
    osm_feature = osm_feature[osm_feature.geom_type.isin(['Point', 'Polygon', 'MultiPolygon'])]#.rename(mapper={'osmid': 'id'}, axis=1)
    #osm_feature['id'] = list(range(len(osm_feature)))
    return osm_feature[columns_to_keep].reset_index().set_crs('EPSG:4326')

def buffered_bbox_as_geog(bbox_wsen, buffer_distance_meters):
    bbox = GeoExtent(bbox_wsen)
    utm_crs = bbox.as_utm_bbox().crs
    bbox_utm = bbox.as_utm_bbox()
    buffered_utm = [
        bbox_utm.bbox[0] - buffer_distance_meters,
        bbox_utm.bbox[1] - buffer_distance_meters,
        bbox_utm.bbox[2] + buffer_distance_meters,
        bbox_utm.bbox[3] + buffer_distance_meters
    ]
    return GeoExtent(buffered_utm, crs=utm_crs).as_geographic_bbox().bbox



CITYDATA_URL = 'https://cities-data-api.wri.org/cities'
citydata = requests.get(CITYDATA_URL).json()
cities = citydata['cities']


for amenity_info in AMENITIES:
    amenity_name, amenity_list, travel_times, travel_distances = amenity_info
    print(f'  {amenity_name}')


    for city in cities:
        city_id = city['id']
        print(city_id)

        # Get boundary
        try:
            url = city['layers_url']['geojson']
            city_admin = gpd.GeoDataFrame.from_file(url)
        except:
            try:
                url = f"https://wri-cities-data-api.s3.us-east-1.amazonaws.com/data/prd/boundaries/geojson/{city['id']}__{city['admin_levels'][0]}.geojson"
                city_admin = gpd.GeoDataFrame.from_file(url)
            except:
                try:
                    url = f"https://wri-cities-data-api.s3.us-east-1.amazonaws.com/data/prd/boundaries/geojson/{city['id']}__{city['admin_levels'][1]}.geojson"
                    city_admin = gpd.GeoDataFrame.from_file(url)
                except:
                    raise Exception(f"No boundary in API for {city['id']}")
        the_city_admin = gpd.GeoDataFrame.from_file(url)
        the_city_urbext = UrbanExtents().get_data(GeoExtent(city_admin.total_bounds)).to_crs('EPSG:4326')
        the_city_urbext = the_city_urbext.dissolve()
        utm_crs = GeoExtent(the_city_admin.total_bounds).as_utm_bbox().crs

    
        LEVELS_TO_DO = ['urbextbound', 'adminbound']
        for level_name in LEVELS_TO_DO:

            print(f'    {level_name}')
            levels = {'urbextbound': the_city_urbext, 'adminbound': the_city_admin}
            id_cols = {'urbextbound': 'city_names', 'adminbound': 'geo_id'}
            level_zones = levels[level_name]


            # Create geodataframe of population-pixel points by vectorizing WorldPop raster. Include only those within the boundary of interest.
            bbox = GeoExtent(level_zones.total_bounds)
            utm_crs = bbox.as_utm_bbox().crs
            worldpop_data = WorldPop(agesex_classes=[]).get_data(bbox)
            wp_df = worldpop_data.drop_vars(['time']).to_dataframe().reset_index()
            pop_points = gpd.GeoDataFrame(wp_df.population, geometry=gpd.points_from_xy(wp_df.x,wp_df.y))
            pop_points_geogr = pop_points.set_crs(utm_crs).to_crs('EPSG:4326')
            # Clip pop_pixels to boundary
            pop_points_clipped = pop_points_geogr.loc[pop_points_geogr.intersects(level_zones.dissolve().geometry[0])]

            # Clip and save road-network PBF to buffered bbox of boundary of interest
            # Requires country-level PBFs to be downloaded and named like XXX.osm.pbf where XXX is three-letter country code
            # Also requires copy of osmconvert in the same directory as the source PBF files
            country_code = [city['country_code_iso3'], 'CONGO'][int(city['country_code_iso3'] in ['COD', 'COG'])]
            filtered_clipped_pbf_filename = f'{city_id}__{level_name}.osm.pbf'
            if not filtered_clipped_pbf_filename in os.listdir(LOCAL_PBF_PREFIX):
                clipped_pbf_filename = f'{city_id}__{level_name}__unfiltered.osm.pbf'
                bbox = buffered_bbox_as_geog(level_zones.total_bounds, 1000)
                print(bbox)
                os.chdir(LOCAL_PBF_PREFIX)
                if not clipped_pbf_filename in os.listdir(LOCAL_PBF_PREFIX):
                    cmd = f"osmconvert64-0.8.8p {country_code}.osm.pbf -b={','.join([str(i) for i in bbox])} --out-pbf -o={filtered_clipped_pbf_filename}"
                    os.system(cmd)

                # Filter so only highway network is included
                # Removed because it seems to cause a file-read error

                #w/highway w/public_transport=platform w/railway=platform w/park_ride r/type=restriction

                # pfile = osmium.FileProcessor(f'{LOCAL_PBF_PREFIX}/{clipped_pbf_filename}')\
                #   .with_filter(osmium.filter.EntityFilter(osmium.osm.WAY))\
                #   .with_filter(osmium.filter.KeyFilter('highway'))
                # if not filtered_clipped_pbf_filename in os.listdir(LOCAL_PBF_PREFIX):
                #     with osmium.SimpleWriter(f'{LOCAL_PBF_PREFIX}/{filtered_clipped_pbf_filename}') as  writer:
                #         for o in pfile:
                #             writer.add(o)
            
            # GET AND SAVE AMENITY POINTS FROM OSM
            oname = f'points__{amenity_name}__{city_id}__{level_name}.geojson'
            if not oname in os.listdir(LOCAL_AMENITYPOINTS_PREFIX):
                # Clip to buffered bbox using command-line osmconvert
                if amenity_name == 'public_open_space':
                    polygon = level_zones.dissolve().reset_index().geometry[0]
                    amenity_points = get_perimeter_points(level_zones, amenity_list[0])  # Change here if we ever have more than one osm_class for this
                else:
                    amenity_points = get_amenities_points(level_zones, amenity_list)
                with open(f'{LOCAL_AMENITYPOINTS_PREFIX}/{oname}', 'w') as ofile:
                    ofile.write(amenity_points.to_json())
                upload_s3(SESSION, f'{LOCAL_AMENITYPOINTS_PREFIX}/{oname}', BUCKET, f'{S3_AMENITYPOINTS_PREFIX}/{oname}')
            else:
                amenity_points = gpd.GeoDataFrame.from_file(f'{LOCAL_AMENITYPOINTS_PREFIX}/{oname}')

            # GENERATE AND SAVE/UPLOAD ISOCHRONE FILE
            transport_network = r5py.TransportNetwork(pathlib.Path(f'{LOCAL_PBF_PREFIX}/{filtered_clipped_pbf_filename}'))
            for idx, travel_mode in enumerate(TRAVEL_MODES):
                travel_mode_name = ['walk', 'bicycle'][idx]
                print(f'      {travel_mode_name}')
                travel_speed = TRAVEL_SPEEDS[travel_mode_name]
                the_travel_times = travel_times[:]
                for travel_dist in travel_distances:
                    the_travel_times.append(travel_dist / travel_speed)
                for time_idx, travel_time in enumerate(the_travel_times):
                    print(travel_time)
                    print(f'{len(amenity_points)}: ', end=' ', flush=True)
                    CHUNKSIZE = 100
                    # res is running total of number-of-accessible-points for each WorldPop pixel point in the AOI
                    res = pd.DataFrame({'to_id': pop_points_clipped.index, 'num_accessible': [0] * len(pop_points_clipped)})

                    # Iterate over amenity points, CHUNKSIZE at a time
                    # Generate TravelTimeMatrix for with those points as origin, and WorldPop pixel points as destinations
                    # Reshape resulting dataframe and update res
                    for chunk in range((len(amenity_points) // CHUNKSIZE) + 1):
                        print(chunk * CHUNKSIZE, end=' ', flush=True)
                        i, j = chunk * CHUNKSIZE, min((chunk + 1) * CHUNKSIZE, len(amenity_points))
                        d = r5py.TravelTimeMatrix(transport_network=transport_network, origins=gpd.GeoDataFrame({'id': range(i, j), 'geometry': amenity_points.iloc[i:j].centroid}), destinations=gpd.GeoDataFrame({'id': range(len(pop_points_clipped)), 'geometry': pop_points_clipped.geometry}), transport_modes=[travel_mode], max_time=datetime.timedelta(minutes=travel_time))
                        d['num_accessible'] = (d.travel_time / d.travel_time).fillna(0)  # Convert non-NAN values to 1, NANs to zero
                        num_accessible = d.pivot(index='to_id', columns='from_id', values='num_accessible').sum(axis=1, skipna=True)
                        res.num_accessible = res.num_accessible + num_accessible

                    # Add num-accessible-point data to the WorldPop pixel points, including those initially clipped off
                    res = res.set_index('to_id')
                    pop_points['accessible_points'] = 0
                    pop_points.loc[pop_points_clipped.index, 'accessible_points'] = res['num_accessible']
                    pop_points.accessible_points = pop_points.accessible_points.fillna(0).astype(int)

                    # Convert to raster and store/upload
                    geo_grid = make_geocube(
                        vector_data=pop_points.set_crs(utm_crs),
                        measurements=['accessible_points'],
                        like = worldpop_data,
                        rasterize_function=rasterize_points_griddata,
                    )

                    if time_idx + 1 > len(travel_times):
                        threshold = travel_time * travel_speed
                        unit = 'meters'
                    else:
                        threshold = travel_times
                        unit = 'minutes'
                    oname = f'{amenity_name}__{city_id}__{travel_mode_name}__{threshold}__{unit}.tif'

                    geo_grid.rio.to_raster(f'{LOCAL_ISOCHRONE_PREFIX}/{oname}', engine='GeoTIFF')
                    upload_s3(SESSION, f'{LOCAL_ISOCHRONE_PREFIX}/{oname}', BUCKET, f'{S3_ISOCHRONE_PREFIX}/{oname}')
                    print()
            print()
