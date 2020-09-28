from geopandas import GeoDataFrame, read_file
import overpy
# import descartes
import pandas
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString
from scipy.special import tandg, cotdg
import numpy
from matplotlib.colors import ColorConverter, LinearSegmentedColormap


# Thanks @Ivan.Baklanov
def get_sequence(start, end, precision):
    if start > end:
        precision = -precision
    return list(numpy.arange(start, end, precision))


class bbox_box():
    xmin = 0
    ymin = 0
    xmax = 0
    ymax = 0
    tpl = () # xmin, ymin, xmax, ymax
    dct = {}
    geo = None
    osm_coords = ''
    name = ''

    def __init__(self, bbox, name=''):
        self.name = str(name)
        self.tpl = bbox
        self.xmin = bbox[0]
        self.ymin = bbox[1]
        self.xmax = bbox[2]
        self.ymax = bbox[3]
        self.dct = self.bbox2dict(self.tpl)
        self.geo = self.frame_draw(self.dct, name)
        # for OSM OVERPASS API we need change the order of lat, lan
        move_coords = [
            str(self.ymin),
            str(self.xmin),
            str(self.ymax),
            str(self.xmax) 
        ]
        self.osm_coords = ','.join(move_coords)
        # print("BBOX", self.tpl, self.osm_coords)

    def __str__(self):
        return "{}: {}, {}, {}, {}".format(self.name if self.name != '' else 'Unnamed', self.xmin, self.ymin, self.xmax, self.ymax)

    def __repr__(self):
        return self.name if self.name != '' else 'Unnamed'

    def bbox2dict(self, bbox):
        return {
            'xmin': bbox[0],
            'ymin': bbox[1],
            'xmax': bbox[2],
            'ymax': bbox[3],
        }

    def frame_draw(self, bbox_dict, name, type='bbox'):
        temp_geodataframe = GeoDataFrame([], columns=['geometry', 'name', 'type'] , crs='EPSG:4326')
        temp_geodataframe.loc[0] = {'name': name, 'type': type, 'geometry': MultiLineString([
            ((bbox_dict['xmin'], bbox_dict['ymin']), (bbox_dict['xmax'], bbox_dict['ymin'])),
            ((bbox_dict['xmax'], bbox_dict['ymin']), (bbox_dict['xmax'], bbox_dict['ymax'])),
            ((bbox_dict['xmax'], bbox_dict['ymax']), (bbox_dict['xmin'], bbox_dict['ymax'])),
            ((bbox_dict['xmin'], bbox_dict['ymax']), (bbox_dict['xmin'], bbox_dict['ymin']))
        ])}
        return temp_geodataframe
        

class coast_part():
    bbox = None
    bbox_real = None
    bbox_broadened = None
    frame_fids = {}
    frame_clusters = []

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
    coastline_geo = None
    waves_geo = None
    ocean_geo = None
    cmap = None


    def __init__(self, file_path, bbox):
        self.bbox = bbox_box(bbox, 'source_bbox')
        self.coastline_geo = read_file(file_path, bbox=bbox)
        self.geo_all.append(self.coastline_geo)
        self.cmap = LinearSegmentedColormap.from_list("", ["green","yellow","red"])
        print(self.bbox)
        # calculate the frames

        # TBD: maybe we need to operate between total_bounds and bounds without any frames?
        # the real frame which is much bigger than which we requested (I don't know the reason)
        xmin = self.coastline_geo.total_bounds.tolist()[0]
        ymin = self.coastline_geo.total_bounds.tolist()[1]
        xmax = self.coastline_geo.total_bounds.tolist()[2]
        ymax = self.coastline_geo.total_bounds.tolist()[3]
        for i, fid in enumerate(self.coastline_geo.FID):
            # for every FID from the shape file we create their own frame
            xmin_frame = self.coastline_geo.geometry[i].bounds[0]
            ymin_frame = self.coastline_geo.geometry[i].bounds[1]
            xmax_frame = self.coastline_geo.geometry[i].bounds[2]
            ymax_frame = self.coastline_geo.geometry[i].bounds[3]
            # print('FID of the object from shapefile:', fid)
              # adding to the dict FID's frame
            self.frame_fids[fid] = bbox_box((xmin_frame, ymin_frame, xmax_frame, ymax_frame), fid)
            # aggregate FIDs by their frames in case we have islands near the coast line. Save the biggest frame which incapsulates the others
            incapsulating_check_success_flag = False
            for i, cluster in enumerate(self.frame_clusters):
                checking = self.check_incapsulation(cluster, self.frame_fids[fid])
                if checking:
                    self.frame_clusters[i] = checking
                    incapsulating_check_success_flag = True
            if incapsulating_check_success_flag is False:
                self.frame_clusters.append(self.frame_fids[fid])

            # update real bbox maximums and minimums
            if xmin_frame < xmin:
                xmin = xmin_frame
            if ymin_frame < ymin:
                ymin = ymin_frame
            if xmax_frame > xmax:
                xmax = xmax_frame
            if ymax_frame > ymax:
                ymax = ymax_frame
        # remove duplicates
        self.frame_clusters = set(self.frame_clusters)
        # for cl in self.frame_clusters:
            # print(cl)
        self.bbox_real = bbox_box((xmin, ymin, xmax, ymax), 'bbox_real')


    def check_incapsulation(self, bbox_1, bbox_2):
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


    # will be updated, works for the first quarter only 8()
    def wave_line(self, xstart, ystart, xend, yend, wave_spec):
        # print(xstart, ystart, xend, yend, wave_spec['angle'])
        if (0 <= wave_spec['angle'] <= 180):
            end_point_x = xend
            end_point_y = ((xend - xstart) * tandg(wave_spec['angle'])) + ystart
            if (end_point_y > yend):
                end_point_y = yend
                end_point_x = ((yend - ystart) * cotdg(wave_spec['angle'])) + xstart
        elif (180 < wave_spec['angle'] <= 360):
            end_point_y = yend
            end_point_x = ((yend - ystart) * tandg(wave_spec['angle'])) + xstart
            if (end_point_x > xend):
                end_point_x = xend
                end_point_y = ((xend - xstart) * cotdg(wave_spec['angle'])) + ystart
        return [(xstart, ystart), (end_point_x, end_point_y)]


    def set_waves(self, angle=45, height=0, period=0):
        # Will be repaced by API call
        self.wave_spec = {
            'angle': angle, 
            'height': height,
            'period': period
        }
        # dangerousness should be calculated
        self.wave_spec['dang'] = 80


    def set_wind(self, angle=45, height=0, period=0):
        # Will be repaced by API call
        self.wind_spec = {
            'angle': angle, 
            'height': height,
            'period': period
        }
        # dangerousness should be calculated
        self.wind_spec['dang'] = 80


    # draw wave lines
    def wave_draw(self, bbox, wave_spec, precision):
        waves_geo = GeoDataFrame([], columns=['geometry'], crs="EPSG:4326")

        if (0 <= wave_spec['angle'] < 90):
            xstart = bbox.xmin
            ystart = bbox.ymin
            xend = bbox.xmax
            yend = bbox.ymax
        elif (90 <= wave_spec['angle'] < 180):
            xstart = bbox.xmax
            ystart = bbox.ymin
            xend = bbox.xmin
            yend = bbox.ymax
        elif (180 <= wave_spec['angle'] < 270):
            xstart = bbox.xmax
            ystart = bbox.ymax
            xend = bbox.xmin
            yend = bbox.ymin
        elif (270 <= wave_spec['angle'] <= 360):
            xstart = bbox.xmin
            ystart = bbox.ymax
            xend = bbox.xmax
            yend = bbox.ymin

        for x in get_sequence(xstart, xend, precision):
            waves_geo.loc[len(waves_geo), 'geometry'] = LineString(self.wave_line(x, ystart, xend, yend, wave_spec))
            # print('X', x, xstart, xend)

        for y in get_sequence(ystart, yend, precision):
            waves_geo.loc[len(waves_geo), 'geometry'] = LineString(self.wave_line(xstart, y, xend, yend, wave_spec))
            # print('Y', y, ystart, yend)

        # print(waves_geo)
        return waves_geo


    def coords_list(self, obj):
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


    def intersection(self, waves, coastline):
        waves_parted = []
        wave_dang = []
        for wave in waves.geometry:
            diff = wave.difference(coastline.unary_union)
            if not diff.is_empty:
                if diff.type == 'LineString':
                    waves_parted.append(diff)
                    wave_dang.append(self.wave_spec['dang'])
                elif diff.type == 'MultiLineString':
                    waves_parted.extend(diff)
                    wave_dang.extend([0 if x>0 else self.wave_spec['dang'] for x in range(len(diff))])
        return GeoDataFrame({'wave_dang': wave_dang, 'type': ['wave' for i in range(len(waves_parted))], 'geometry': waves_parted})


    def combination(self, geos):
        # print(geos)
        return GeoDataFrame(pandas.concat(geos, ignore_index=True))


    def set_towns(self, bbox, place_regexp='city|town|village|hamlet'):
        api = overpy.Overpass()
#         print(f'''
# (
#   node
#   ["place"~"{place_regexp}"]
#     ({bbox.osm_coords});
# )->._;
# (._;>;);
# out;''')

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


    def ocean_plot(self, precision=0.0001, show_towns=False, show_bboxes=False, show_frames=False):
        self.precision = precision

        # enlarging the full frame or we won't have waves at the protrusive points
        enlarging = 0.01
        self.bbox_broadened = bbox_box(
            (self.bbox_real.xmin - enlarging,
            self.bbox_real.ymin - enlarging, 
            self.bbox_real.xmax + enlarging, 
            self.bbox_real.ymax + enlarging),
            'bbox_broadened'
        )

        for cluster in self.frame_clusters:
            waves_cluster_geo = self.wave_draw(cluster, self.wave_spec, self.precision)
            waves_parted = self.intersection(waves_cluster_geo, self.coastline_geo)
            self.geo_all.append(waves_parted)
        # self.waves_geo = self.wave_draw(self.bbox_broadened, self.wave_spec, self.precision)
        # waves_parted = self.intersection(self.waves_geo, self.coastline_geo)
        # self.geo_all.append(waves_parted)

        if show_bboxes is True:
            self.geo_all.extend([self.bbox_real.geo, self.bbox.geo])

        if show_towns is True:
            towns = self.set_towns(self.bbox_real, place_regexp='city')
            # print(towns)
            self.geo_all.append(towns)

        if show_frames is True:
            # print([self.frame_draw(self.frame_fids[frame], frame) for frame in self.frame_fids])
            self.geo_all.extend([self.frame_fids[frame].geo for frame in self.frame_fids])

        self.ocean_geo = self.combination(self.geo_all)
        # print(self.ocean_geo)
        print(self.bbox_real)

        # self.ocean_geo.plot(legend=True, column='wave_dang', cmap=self.cmap, vmin=0, vmax=100, missing_kwds = {'color': 'tan', "edgecolor": 'darkgoldenrod'})
        self.ocean_geo.plot(legend=True, column='wave_dang', cmap=self.cmap, vmin=0, vmax=100, missing_kwds = {'color': 'tan', "edgecolor": 'black'})
        plt.annotate(
            text='Wave angle: %s\nPrecision: %s' % (self.wave_spec['angle'], self.precision),
            xy=(self.bbox_real.xmin, self.bbox_real.ymax),
            verticalalignment='top'
        )
        
        # city names
        if show_towns is True:
            for x, y, name in zip(towns.geometry.x, towns.geometry.y, towns.name):
                plt.annotate(name, xy=(x, y), xytext=(3, 3), textcoords='offset points', color='darkblue')

        plt.title('Waves and the coastline intersection')
        plt.show()



# bbox = (-9.48859, 38.71225, -9.48369, 38.70596)
bbox = (-9.48859,38.70044,-9.4717541,38.7284016)
# shape_file = '/home/maksimpisarenko/tmp/osmcoast/coastlines-split-4326/lines.shp'
shape_file = '/home/maksimpisarenko/tmp/osmcoast/land-polygons-split-4326/land_polygons.shp'
cascais = coast_part(shape_file, bbox)
cascais.set_waves(angle=220)
cascais.set_wind()
cascais.ocean_plot(precision=0.1, show_towns=True, show_bboxes=False, show_frames=True)
