import os, sys, json, requests
os.environ['GOOGLE_APPLICATION_USER'] = 'script-service-account@wri-gee.iam.gserviceaccount.com'
os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = 'C:/Users/tgwon/.google/credkey.json'
import geopandas as gpd
import pandas as pd
import osmnx as ox
#os.chdir('C:/Users/tgwon/wri/cif/cities-cif')
sys.path = ['C:\\Users\\tgwon\\wri\\indicators', 'C:\\Users\\tgwon\\wri\\cif\\cities-cif', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\python310.zip', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\DLLs', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib\\site-packages', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib\\site-packages\\win32', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib\\site-packages\\win32\\lib', 'C:\\Users\\tgwon\\anaconda3\\envs\\cities-cif\\lib\\site-packages\\Pythonwin',]
from city_metrix.metrix_model import GeoExtent
from city_metrix.layers import OpenStreetMap, OpenStreetMapClass, UrbanExtents
FILEPATH = 'C:/Users/tgwon/wri/indicators/demodata'

TARGET_AMENITYTYPES = [OpenStreetMapClass.COMMERCE]


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

    return osm_feature

def get_park_perimeter_points(city_gdf):
    osm_class = OpenStreetMapClass.OPEN_SPACE
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
    print(len(openspace_polys), end=': ')
    sys.stdout.flush()
    for rownum in list(range(len(openspace_polys))):
        if rownum % 1000 == 0:
            print(rownum, end=' ')
            sys.stdout.flush()
        multi = openspace_polys.iloc[[rownum]].to_crs(f'EPSG:{utm_epsg}').buffer(5).to_crs('EPSG:4326').boundary.intersection(road_uu)[rownum]
        if multi.geom_type=='MultiPoint':
            perimeter_pts += multi.geoms
        else:
            perimeter_pts.append(multi)
    res = gpd.GeoDataFrame({'id': range(len(perimeter_pts)), 'geometry': perimeter_pts}).set_crs('EPSG:4326').set_geometry('geometry')
    adminboundary_gdf = city_gdf.dissolve()
    a_points_filtered = res.loc[res.within(adminboundary_gdf.geometry[0])].reset_index()
    a_points_filtered_pointsonly = a_points_filtered.loc[a_points_filtered.geometry.geom_type=='Point']
    return a_points_filtered_pointsonly

def merge_osm_classes(osm_classes):
    result = {}
    for d in osm_classes:
        d = d.value
        for k in d:
            if not (k in list(result.keys())):
                result[k] = []
            if result[k] == True or d[k] == True:
                result[k] = True
            else:
                result[k] += d[k]
    return result

def get_amenities_pointsonly(city_gdf, osm_classes):
    bbox = GeoExtent(city_gdf.total_bounds)
    utm_epsg = bbox.as_utm_bbox().epsg_code
    polygon = city_gdf.dissolve().geometry[0]
    merged_osm_dicts = merge_osm_classes(osm_classes)
    # Set the OSMnx configuration to disable caching
    ox.settings.use_cache = False
    try:
        #osm_feature = ox.features_from_bbox(bbox=(min_lon, min_lat, max_lon, max_lat), tags=self.osm_class.value)
        osm_feature = ox.features_from_polygon(polygon, tags=merged_osm_dicts)
    # When no feature in bbox, return an empty gdf
    except ox._errors.InsufficientResponseError as e:
        osm_feature = gpd.GeoDataFrame(pd.DataFrame(columns=['osmid', 'geometry']+list(merged_osm_dicts.keys())), geometry='geometry')
        osm_feature.crs = "EPSG:4326"

    osm_feature = osm_feature[osm_feature.geom_type.isin(['Point', 'Polygon', 'MultiPolygon'])].reset_index().rename(mapper={'id': 'osmid'}, axis=1)
    # keep only columns desired to reduce file size
    keep_col = ['osmid', 'geometry'] + list(merged_osm_dicts.keys())
    for col in keep_col:
        if not col in osm_feature.columns:
            osm_feature[col] = [pd.NA] * len(osm_feature)

    osm_feature = osm_feature[keep_col]
    osm_feature.geometry = osm_feature.to_crs(utm_epsg).centroid.to_crs('epsg:4326')
    
    result = {}
    for osm_class in osm_classes:
        class_name = osm_class.__str__().lower().split('.')[1].replace('_', '-')
        osm_dict = osm_class.value
        to_keep = ['osmid', 'geometry'] + list(osm_dict.keys())
        result[class_name] = gpd.GeoDataFrame(columns=to_keep, geometry='geometry').set_crs('EPSG:4326')
        for k in osm_dict:
            if osm_dict[k] == True:
                to_append = osm_feature.loc[osm_feature[k].notnull()][to_keep]
            else:
                to_append = osm_feature.loc[osm_feature[k].isin(osm_dict[k])][to_keep]
            result[class_name] = pd.concat([result[class_name], to_append])

    return result


citydata_url = 'https://cities-data-api.wri.org/cities'
citydata = requests.get(citydata_url).json()

dd_cities = [c for c in citydata['cities'] if 'deepdive' in c['projects']]
easy_cities = [c for c in dd_cities if (not c['country_code_iso3'] in ['CHN', 'IDN', 'IND'])]

amenityname = 'commerce'

print(amenityname)
#for city_idx in range(len(easy_cities)-1, -1, -1):
#    city = easy_cities[city_idx]
for city in citydata['cities']:
    print(city['name'])
    try:
        url = gpd.GeoDataFrame.from_file(city['layers_url']['geojson'])
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
    city_urbext = UrbanExtents().get_data(GeoExtent(city_admin.total_bounds)).to_crs('epsg:4326')
    for boundary in [('admin', city_admin), ('urbext', city_urbext)]:
        boundaryname = boundary[0]
        boundary_gdf = boundary[1]
        ofilename = f"{amenityname}__{boundaryname}bound__{city['id']}.geojson"
        if not ofilename in os.listdir(FILEPATH + "/amenitypoints"):
            with open(f"{FILEPATH}/boundaries/{boundaryname}bound__{city['id']}.geojson", 'w') as ofile:
                ofile.write(boundary_gdf.to_json())
                all_results = []
                for amenity in TARGET_AMENITYTYPES:
                    amenity_result_oneamenity = list(get_amenities_pointsonly(boundary_gdf, [amenity]).values())[0]
                    amenity_result_oneamenity['amenity_class'] = amenityname
                    amenity_result_oneamenity = amenity_result_oneamenity[['osmid', 'geometry']].reset_index()
                all_results.append(amenity_result_oneamenity)
                amenity_results = pd.concat(all_results)
                filtered_results = amenity_results.loc[amenity_results.within(boundary_gdf.dissolve().geometry[0])].reset_index()
                with open(f"{FILEPATH}/amenitypoints/{ofilename}", 'w') as ofile:
                    ofile.write(filtered_results.to_json())
        print()
