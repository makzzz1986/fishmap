# for the debug purpose
from time import time
########

import os
import numpy
import matplotlib.pyplot as plt
import overpy
from geopandas import GeoDataFrame, read_file
from pandas import concat
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString, Polygon
from scipy import special as sc
from matplotlib.colors import ColorConverter, LinearSegmentedColormap
from typing import List
from datetime import datetime
# import descartes
from helpers import *


# Thanks @Ivan.Baklanov
def get_sequence(start, end, precision):
    if start > end:
        precision = -precision
    return list(numpy.arange(start, end, precision))


class WaveMap():
    bbox = None
    # bbox_real = None
    bbox_broadened = None

    wave_spec = {
        'angle': 0, 
        'height': 0,
        'period': 0,
        'dang': 100
    }
    wind_spec = {
        'angle': 0, 
        'height': 0,
        'period': 0,
        'dang': 100
    }
    precision = 0
    geo_all = []
    coastline_union = None
    coastline_geo = None
    waves_geo = None
    ocean_geo = None
    cmap = None
    tile_dens = 0.


    def __init__(self, file_path, bbox):
        self.bbox = Bbox(bbox, 'source_bbox')
        coastline_geo_temp = read_file(file_path, bbox=bbox, ignore_fields=['FID'])
        coastline_polygon = coastline_geo_temp.unary_union.intersection(self.bbox.polygon)
        self.coastline_geo = GeoDataFrame({'geometry': coastline_polygon})
        self.geo_all.append(self.coastline_geo)
        self.coastline_union = self.coastline_geo.unary_union
        self.cmap = LinearSegmentedColormap.from_list("", ["green","yellow","red"])
        # self.bbox_real = Bbox(self.coastline_union.bounds, 'bbox_real')
        print(self.bbox)


    def check_incapsulation(self, bbox_1, bbox_2) -> Bbox:
        # bbox_1 encapsulates bbox_2
        if  (bbox_1.xmin < bbox_2.xmin) and \
            (bbox_1.ymin < bbox_2.ymin) and \
            (bbox_1.xmax > bbox_2.xmax) and \
            (bbox_1.ymax > bbox_2.ymax):
            # print(bbox_1, '>', bbox_2)
            return bbox_1
        # bbox_2 encapsulates bbox_1
        elif (bbox_2.xmin < bbox_1.xmin) and \
             (bbox_2.ymin < bbox_1.ymin) and \
             (bbox_2.xmax > bbox_1.xmax) and \
             (bbox_2.ymax > bbox_1.ymax):
            # print(bbox_2, '>', bbox_1)
            return bbox_2
        # else:
            # print(bbox_1, '~', bbox_2)


    def wave_line(self, xstart, ystart, xend, yend, angle, quart) -> List:
        # print(xstart, ystart, xend, yend, wave_spec['angle'], sc.tandg(wave_spec['angle']), sc.cotdg(wave_spec['angle']))           
        # the I and II quarters
        if (quart == 1) or (quart == 2):
            end_point_x = xend
            end_point_y = ((xend - xstart) * sc.tandg(angle)) + ystart
            # if Y coord out of frame - draw from cotn
            if (end_point_y > yend):
                end_point_y = yend
                end_point_x = ((yend - ystart) * sc.cotdg(angle)) + xstart
        # the III quarter 
        elif quart == 3:
            end_point_x = xend
            end_point_y = ((xstart - xend) * sc.tandg(angle)) + ystart
            # if Y coord out of frame - draw from cotn
            if (end_point_y > yend):
                end_point_y = yend
                end_point_x = ((yend - ystart) * sc.cotdg(angle)) + xstart
            # if X coord out of frame - draw from tan from... smth
            if (end_point_x < xend):
                end_point_x = xend
                end_point_y = ((xend - xstart) * sc.tandg(angle)) + ystart
        # the IX quarter
        elif quart == 4:
            end_point_y = yend
            end_point_x = ((yend - ystart) * sc.cotdg(angle)) + xstart
            if (end_point_x > xend):
                end_point_x = xend
                end_point_y = ((xend - xstart) * sc.tandg(angle)) + ystart
        return [(xstart, ystart), (end_point_x, end_point_y)]


    # draw wave lines
    def wave_draw(self, bound, wave_spec, precision) -> GeoDataFrame:
        waves_geo = GeoDataFrame({'geometry': [], 'wave_dang': []}, crs="EPSG:4326")
        if (0 <= wave_spec['angle'] < 90):
            xstart = bound[0]
            ystart = bound[1]
            xend = bound[2]
            yend = bound[3]
            quart = 1
        elif (90 <= wave_spec['angle'] < 180):
            xstart = bound[2]
            ystart = bound[1]
            xend = bound[0]
            yend = bound[3]
            quart = 2
        elif (180 <= wave_spec['angle'] < 270):
            xstart = bound[2]
            ystart = bound[3]
            xend = bound[0]
            yend = bound[1]
            quart = 3
        elif (270 <= wave_spec['angle'] <= 360):
            xstart = bound[0]
            ystart = bound[3]
            xend = bound[2]
            yend = bound[1]
            quart =4

        for x in get_sequence(xstart, xend, precision):
            waves_geo.loc[len(waves_geo)] = {'geometry': LineString(self.wave_line(x, ystart, xend, yend, wave_spec['angle'], quart)),
                                             'wave_dang': wave_spec['dang']}
            # print('X', x, xstart, xend)

        for y in get_sequence(ystart, yend, precision):
            waves_geo.loc[len(waves_geo)] = {'geometry': LineString(self.wave_line(xstart, y, xend, yend, wave_spec['angle'], quart)),
                                             'wave_dang': wave_spec['dang']}
            # print('Y', y, ystart, yend)

        # print(waves_geo)
        return waves_geo


    def coords_list(self, obj) -> List:
        if obj.type == 'LineString':
            return [coord for coord in obj.coords]
        elif obj.type == 'MultiLineString':
            temp_list = []
            for line in obj:
                temp_list.extend([coord for coord in line.coords])
            return temp_list
        else:
            print('WTF?', obj.type, obj)
            return None


    def intersection(self, waves, coastline) -> GeoDataFrame:
        waves_parted = []
        wave_dang = []
        for wave in waves.itertuples():
            diff = wave.geometry.difference(coastline)
            if not diff.is_empty:
                if diff.type == 'LineString':
                    waves_parted.append(diff)
                    wave_dang.append(wave.wave_dang)
                elif diff.type == 'MultiLineString':
                    waves_parted.extend(diff)
                    wave_dang.extend([0 if x>0 else wave.wave_dang for x in range(len(diff))])
        return GeoDataFrame({'wave_dang': wave_dang, 'type': ['wave' for i in range(len(waves_parted))], 'geometry': waves_parted})


    def combination(self, geos) -> GeoDataFrame:
        # print(geos)
        return GeoDataFrame(concat(geos, ignore_index=True))


    def set_towns(self, bbox, place_regexp='city|town|village|hamlet') -> GeoDataFrame:
        api = overpy.Overpass()
        print(f'''
(
  node
  ["place"~"{place_regexp}"]
    ({bbox.osm_coords});
)->._;
(._;>;);
out;''')

        result = api.query(f'''
(
  node
  ["place"~"{place_regexp}"]
    ({bbox.osm_coords});
)->._;
(._;>;);
out;''')
        towns_points_coord = []
        towns_points_names = []
        print(f'Grab from Overpass {str(len(result.nodes))} objects')
        for node in result.nodes:
            # print(node.tags['name'], node.lat, node.lon)
            # print(node.tags)
            towns_points_names.append(node.tags['name'])
            towns_points_coord.append(Point(node.lon,node.lat))
        return GeoDataFrame({'name': towns_points_names, 
                            'type': ['town' for i in range(len(towns_points_names))], 
                            'geometry': towns_points_coord})


    def tiling(self, bounds, tile_dens, filter=0.01) -> GeoDataFrame:
        # creating the matrix of tiles through all the map
        parts_matrix = GeoDataFrame({'geometry': [], 'type': [], 'name': []})
        # finding the largest axis and divide it to density
        if abs(bounds[0]-bounds[2])/abs(bounds[1]-bounds[3]) > 1:
            side_length = abs(bounds[0]-bounds[2])/tile_dens
        else:
            side_length = abs(bounds[1]-bounds[3])/tile_dens
        # print('Side length', side_length)
        ## get length segments from horizontal and vertical
        x_notch = get_sequence(bounds[0], bounds[2], side_length)
        x_notch.append(bounds[2])
        y_notch = get_sequence(bounds[1], bounds[3], side_length)
        y_notch.append(bounds[3])
        ## fill all map by tiles
        matrix_counter = 0
        for y in range(len(y_notch) - 1):
            for x in range(len(x_notch) - 1):
                # print(x_notch[x], y_notch[y], x_notch[x+1], y_notch[y+1])
                # filter out very narrow tiles
                if (abs(x - (x+1)) > filter) or (abs(y - (y+1)) > filter):
                    # print(f'Tile: ({x_notch[x]}, {y_notch[y]}), ({x_notch[x+1]}, {y_notch[y]}), ({x_notch[x+1]}, {y_notch[y+1]}), ({x_notch[x]}, {y_notch[y+1]})')
                    parts_matrix.loc[matrix_counter] = {'geometry': Polygon(((x_notch[x]  , y_notch[y]), 
                                                                             (x_notch[x+1], y_notch[y]), 
                                                                             (x_notch[x+1], y_notch[y+1]), 
                                                                             (x_notch[x]  , y_notch[y+1]))), 
                                                        'name': f'{str(x)}x{str(y)}',
                                                        'type': 'tile'}
                    matrix_counter += 1
        return parts_matrix


    # in case we have a figure like L, |, Г, - . This figures have perpendicular edges and were appeared after splitting source
    def checking_tetris_shape(self, polygon) -> bool:
        points = list(zip(polygon.exterior.xy[0], polygon.exterior.xy[1]))
        for i, p in enumerate(points[:-1]):
            # if it has only perpendicular sides - the points go one after another with change only _one_ axis
            if p[0] != points[i+1][0] and p[1] != points[i+1][1]:
                return False
        return True


    def tiles_coast_diff(self, matrix, coast_union) -> GeoDataFrame:
        tiles_diff = GeoDataFrame({'geometry': [], 'type': [], 'name': []}, crs='EPSG:4326')
        # start_time = time()
        for row in matrix.itertuples():
            pile_cutted = coast_union.intersection(row.geometry) 
            
            # if the Polygon is a parallelogram (has 5 points (4 vertices + 1 dublicates the start point))
            # we can filter it, because it isn't connected to the ocean, but is on the edge of map
            ## for Polygons filtering is easy:
            if pile_cutted.type == 'Polygon':
                points_quantity = len(pile_cutted.exterior.coords.xy[0])
                # it is just a parallelogram
                if points_quantity == 5:
                    pass
                # it is a shape like L | Г -  - source artifacts
                elif (points_quantity < 10) and (self.checking_tetris_shape(pile_cutted)):
                    pass
                else:
                    tiles_diff.loc[len(tiles_diff)] = {'geometry': pile_cutted, 
                                                    'type': 'coast_cut',
                                                    'name': row.name}
            ## in case we have pile in the middle of an age which consists of a few polygons,
            ## we should check that all polygons in that tile have the len of 5 vertices
            elif pile_cutted.type == 'MultiPolygon':
                list_pole_length = [len(poly.exterior.coords.xy[0]) for poly in pile_cutted]
                if list_pole_length != [5 for l in range(len(list_pole_length))]:
                    tiles_diff.loc[len(tiles_diff)] = {'geometry': pile_cutted, 
                                                       'type': 'coast_cut',
                                                       'name': row.name}
        # print(f'--- {str(time() - start_time)} seconds ---')
        # for i in tiles_diff.itertuples():
        #     if i.geometry.type == 'Polygon':
        #         print(f'This is {str(i.geometry.type)}, name: {str(i.name)}, lens {str(len(i.geometry.exterior.coords.xy[0]))}')
        #     elif i.geometry.type == 'MultiPolygon':
        #         print(f'This is {str(i.geometry.type)}, name: {str(i.name)}, lens {str([len(poly.exterior.coords.xy[0]) for poly in i.geometry])}')
        return tiles_diff


    def splitting_map(self, geo, bounds, tile_dens=1) -> GeoDataFrame:
        tiles = self.tiling(bounds, tile_dens, self.precision)
        ## no need in soil polygons which don't touch the ocean. The source of soil consist of polygons, splitted
        ## approximetely by 1 degree Longtitude and Latitude. If we have an area of ~1 - this polygon is surrounded by soil,
        ## no need to check waves on it. So, we can exclude them from the map of waves
        ## ignore UserWarning about area and CRS. We don't care about geographical area!
        without_big_soil_filter = geo['geometry'].area<1
        without_big_soil = geo[without_big_soil_filter]
        without_big_soil_union = without_big_soil.unary_union

        if len(tiles) > 1: # it means we have the tile != map
            overlaping = tiles.overlaps(without_big_soil_union)
            tiles_overlaped = tiles[overlaping]
        else: # if we have tile == map
            tiles_overlaped = tiles
        print('Tiles:', str(len(tiles_overlaped)))
        return self.tiles_coast_diff(tiles_overlaped, without_big_soil_union)


    def bbox_broading(self, frame, adding_lenght=0.1, name='bbox enlarged') -> Bbox:
        bbox_enlarge = Bbox((frame[0]-adding_lenght,
                                frame[1]-adding_lenght,
                                frame[2]+adding_lenght,
                                frame[3]+adding_lenght),
                                name=name)
        print('Enlarging from', frame, 'to', bbox_enlarge)
        return bbox_enlarge


    def ocean_calculating(self, precision=0.0001, tile_dens=1, debug=False) -> None:
        # if tile_dens < than 1 drop it to 1
        if tile_dens < 1:
            tile_dens = 1

        self.precision = precision
        self.tile_dens = tile_dens

        tiles = self.splitting_map(self.coastline_geo, self.coastline_union.bounds, self.tile_dens)
        for tile in tiles.geometry:
            # getting coordinates of tile's center
            centroid_coordinates = tile.centroid.coords
            lon = centroid_coordinates.xy[0][0]
            lat = centroid_coordinates.xy[1][0]
            # getting center time zone (API or fixed Portugal TZ)
            tz = Time(lon, lat, {'days': 1, 'hours': 8})
            tz_time = tz.get_static()
            # getting waves' specification - angle, height, period and dangerousness
            wave = Wave(lon, lat, tz_time['timestamp'])
            if debug is True:
                wave_spec = wave.get_random(angle=45)
            else:
                wave_spec = wave.get_stormglass(os.environ['WEATHER_API_TOKEN'], debug=True)

            waves_tile_geo = self.wave_draw(tile.bounds, wave_spec, self.precision)
            waves_parted = self.intersection(waves_tile_geo, tile)
            self.geo_all.append(waves_parted)

        self.ocean_geo = self.combination(self.geo_all)


    def plot(self, show_towns=False, show_bboxes=False) -> None:
        if show_bboxes is True:
            self.ocean_geo[len(self.ocean_geo)] = self.bbox.geo

        if show_towns is True:
            # towns = self.set_towns(self.bbox, place_regexp='city')
            towns = self.set_towns(self.bbox, place_regexp='city|town')
            self.ocean_geo = self.combination([self.ocean_geo, towns])

        # self.ocean_geo.plot(legend=True, column='wave_dang', cmap=self.cmap, vmin=0, vmax=100, missing_kwds = {'color': 'tan', "edgecolor": 'darkgoldenrod'})
        self.ocean_geo.plot(legend=True, column='wave_dang', cmap=self.cmap, vmin=0, vmax=100, missing_kwds = {'color': 'tan', "edgecolor": 'black'})
        plt.annotate(
            text=f'Precision: {str(self.precision)}\nTiling: {str(self.tile_dens)}',
            xy=(self.bbox.xmin, self.bbox.ymax),
            verticalalignment='top'
        )

        # city names
        if show_towns is True:
            for x, y, name in zip(towns.geometry.x, towns.geometry.y, towns.name):
                plt.annotate(name, xy=(x, y), xytext=(3, 3), textcoords='offset points', color='darkblue')

        plt.gca().set_aspect('equal', adjustable='box')
        # plt.axis('scaled')
        plt.title('Waves and the coastline intersection')
        plt.show()


    def save_to_file(self, file_prefix):
        now = datetime.now().isoformat()
        filepath = f'{file_prefix}-{now}.geojson'
        self.ocean_geo.to_file(filepath, driver='GeoJSON')
        print('File saved to', filepath)


# bbox = (-9.48859,38.70044,-9.4717541,38.7284016)
bbox = (-9.8,38.1,-9.1,39.2) # around Lisbon
# bbox = (-9.3, 38.5, -9.1, 38.7) # Costa Caparica
# bbox = (-8.0,36.0,-10.0,42.0)  # VERY BIG!
shape_file = '/home/maksimpisarenko/tmp/osmcoast/land-polygons-split-4326/land_polygons.shp'
portugal = WaveMap(shape_file, bbox)

portugal.ocean_calculating(precision=0.0001, tile_dens=4, debug=True)
portugal.plot(show_towns=True, show_bboxes=False)
portugal.save_to_file('ready_shapes/portugal')
