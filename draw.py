import geopandas
import overpy
import pandas
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString
from scipy.special import tandg, cotdg
import numpy
from matplotlib.colors import ColorConverter, LinearSegmentedColormap


class coast_part():
    bbox = [] # xmin, ymin, xmax, ymax
    bbox_dict = {}
    bbox_geo = None
    bbox_real = [] # xmin, ymin, xmax, ymax
    bbox_real_dict = {}
    bbox_real_geo = None
    frame_fids = {}

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
        self.bbox = bbox
        self.bbox_dict = self.bbox2dict(bbox)
        self.bbox_geo = self.frame_draw(self.bbox_dict)
        self.coastline_geo = geopandas.read_file(file_path, bbox=bbox)
        self.geo_all.append(self.coastline_geo)
        self.cmap = LinearSegmentedColormap.from_list("", ["green","yellow","red"])

        # Getting the real bbox! It is much bigger than bbox
        xmin = bbox[0]
        xmax = bbox[2]
        ymin = bbox[1]
        ymax = bbox[3]

        for _, fid, r in self.coastline_geo.itertuples():
            xmin_frame = xmin
            xmax_frame = xmax
            ymin_frame = ymin
            ymax_frame = ymax
            print('FID of the object from shapefile:', fid)
            for pair in list(r.coords):
                if pair[0] > xmax_frame:
                    xmax_frame = pair[0]
                if pair[0] < xmin_frame:
                    xmin_frame = pair[0]
                if pair[1] > ymax_frame:
                    ymax_frame = pair[1]
                if pair[1] < ymin_frame:
                    ymin_frame = pair[1]
            self.frame_fids[fid] = self.bbox2dict([xmin_frame, ymin_frame, xmax_frame, ymax_frame])
            if xmax_frame > xmax:
                xmax = xmax_frame
            if xmin_frame < xmin_frame:
                xmin_frame = xmin_frame
            if ymax_frame > ymax:
                ymax = ymax_frame
            if ymin_frame < ymin:
                ymin = ymin_frame

        self.bbox_real = (xmin, ymin, xmax, ymax)
        self.bbox_real_dict = self.bbox2dict(self.bbox_real)
        self.bbox_real_geo = self.frame_draw(self.bbox_real_dict)
        

    def frame_draw(self, bbox_dict):
        temp_geodataframe = geopandas.GeoDataFrame([], columns=['geometry', 'name'] , crs='EPSG:4326')
        temp_geodataframe.loc[0] = {'name':'bbox_real_frame', 'geometry': MultiLineString([\
            ((bbox_dict['xmin'], bbox_dict['ymin']), (bbox_dict['xmax'], bbox_dict['ymin'])),\
            ((bbox_dict['xmax'], bbox_dict['ymin']), (bbox_dict['xmax'], bbox_dict['ymax'])),\
            ((bbox_dict['xmax'], bbox_dict['ymax']), (bbox_dict['xmin'], bbox_dict['ymax'])),\
            ((bbox_dict['xmin'], bbox_dict['ymax']), (bbox_dict['xmin'], bbox_dict['ymin']))\
        ])}
        return temp_geodataframe


    def bbox2dict(self, bbox):
        return {
            'xmin': bbox[0],
            'ymin': bbox[1],
            'xmax': bbox[2],
            'ymax': bbox[3],
        }

    # will be updated, works for the first quarter only 8()
    def wave_line(self, start_point, wave, bbox):
        xmax = 0
        ymax = 0
        if (0 < wave['angle'] < 90):
            xmax = bbox['xmax'] # remake them to bbox_dict
            ymax = bbox['ymax']
        elif (90 < wave['angle'] <= 180):
            xmax = bbox['xmin']
            ymax = bbox['ymax']
        elif (180 < wave['angle'] < 270):
            xmax = bbox['xmin']
            ymax = bbox['ymin']
        elif (270 < wave['angle'] <= 360):
            xmax = bbox['xmax']
            ymax = bbox['ymin']
        elif (wave['angle'] == 90) or (wave['angle'] == 270):
            # print([start_point, (start_point[0], bbox['ymax'])])
            return [start_point, (start_point[0], bbox['ymax'])]
        # print(start_point, xmax, ymax)
        end_point_x = xmax
        end_point_y = ((xmax - start_point[0]) * tandg(wave['angle'])) + start_point[1]
        if end_point_y > ymax:
            # print('Too big!')
            end_point_y = ymax
            end_point_x = ((ymax - start_point[1]) * cotdg(wave['angle'])) + start_point[0]
        # print(start_point, (end_point_x, end_point_y), wave, bbox)
        return [start_point, (end_point_x, end_point_y)]


    def wave_parted(self, wave, intersect):
        # drawing parted line
        if intersect.type == 'Point':
            return {'waves': [LineString([wave.coords[0], intersect])], 'wave_dang': [self.wave_spec['dang']]}
        if intersect.type == 'MultiPoint':
            points_list = []
            # multilinestrings = []
            line_list = [] # list with lines
            points_list.append(wave.coords[0])
            points_list.extend([[point.x, point.y] for point in intersect])

            if len(intersect) == 0:
                print('No intersect?')
                return {'waves': [], 'wave_dang': []}
            # if it is odd, than it ends on the ground
            if len(intersect) % 2 == 0:
                print('Odd?')
                for i in intersect:
                    print(i)
                print(points_list)
                points_list.append(wave.coords[-1])
                print(points_list)
            # pair dots to lines!
            for pair in range(0, len(points_list), 2):
                line_list.append(LineString([points_list[pair], points_list[pair+1]]))
            # return list of wave LineStrings and their dangerously, wave_dang for the first one and 0 for others, 
            # because they are after ground.
            return {'waves': line_list, 'wave_dang': [0 if x>0 else self.wave_spec['dang'] for x in range(len(line_list))]} 


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
        waves_geo = geopandas.GeoDataFrame([], columns=['geometry'], crs="EPSG:4326")

        xstep = bbox['xmin']
        ystep = bbox['ymin']
        while True:
            waves_geo.loc[len(waves_geo), 'geometry'] = LineString(self.wave_line((xstep, bbox['ymin']), wave_spec, bbox))
            xstep += precision
            if xstep >= bbox['xmax']:
                break
        while True:
            waves_geo.loc[len(waves_geo), 'geometry'] = LineString(self.wave_line((bbox['xmin'], ystep), wave_spec, bbox))
            ystep += precision
            if ystep >= bbox['ymax']:
                break
        return waves_geo


    def intersection(self, waves, coastline):
        intersected = {'waves': [], 'wave_dang': []}
        for _, wave_line in waves.itertuples(): # checking each line for intersection to every coastline part (it can be islands!)
            full_intersect_points = MultiPoint(()) # all Points and MultiPoints
            for _, _, coast_line in coastline.itertuples():
                intersect = wave_line.intersection(coast_line)
                full_intersect_points = full_intersect_points.union(intersect) # refresh variable with all Points, it is like list.extend 
            if not full_intersect_points.is_empty:
                wave_parts = self.wave_parted(wave_line, full_intersect_points)
                intersected['waves'].extend(wave_parts['waves'])
                intersected['wave_dang'].extend(wave_parts['wave_dang'])
            
        # intersection_points = geopandas.GeoDataFrame(geometry=intersection_list)  # points of intercestion wave and coastline
        intersection = geopandas.GeoDataFrame(intersected['wave_dang'], geometry=intersected['waves'], columns=['wave_dang'])

        return intersection


    def combination(self, geos):
        # print(geos)
        return geopandas.GeoDataFrame(pandas.concat(geos, ignore_index=True))


    # Did I miss something? I have to convert coordinates to lon and lat
    def convert_bbox(self, bbox):
        changed_order = [str(bbox[1]), str(bbox[0]), str(bbox[3]), str(bbox[2])]
        return ','.join(changed_order)


    def set_towns(self, bbox, place_regexp='city|town|village|hamlet'):
        api = overpy.Overpass()
        converted_bbox = self.convert_bbox(bbox)
#         print(f'''
# (
#   node
#   ["place"~"{place_regexp}"]
#     ({converted_bbox});
# )->._;
# (._;>;);
# out;''')
        result = api.query(f'''
(
  node
  ["place"~"{place_regexp}"]
    ({converted_bbox});
)->._;
(._;>;);
out;''')
        towns_points_coord = []
        towns_points_names = []
        for node in result.nodes:
            # print(node.tags['name'], node.lat, node.lon)
            towns_points_names.append(node.tags['name'])
            towns_points_coord.append(Point(node.lon,node.lat))
        return geopandas.GeoDataFrame(towns_points_names, geometry=towns_points_coord, columns=['name'])


    def ocean_plot(self, precision=0.0001, show_towns=False, show_bboxes=False, show_frames=False):
        self.precision = precision
        self.waves_geo = self.wave_draw(self.bbox_real_dict, self.wave_spec, self.precision)
        
        # waves_parted = self.intersection(self.coastline_geo, self.waves_geo)
        waves_parted = self.intersection(self.waves_geo, self.coastline_geo)
        self.geo_all.append(waves_parted)

        if show_bboxes is True:
            self.geo_all.extend([self.bbox_real_geo, self.bbox_geo])

        if show_towns is True:
            towns = self.set_towns(self.bbox_real)
            self.geo_all.append(towns)

        if show_frames is True:
            print([self.frame_draw(self.frame_fids[frame]) for frame in self.frame_fids])
            self.geo_all.extend([self.frame_draw(self.frame_fids[frame]) for frame in self.frame_fids])

        self.ocean_geo = self.combination(self.geo_all)
        # print(self.coastline_geo)
        self.ocean_geo.plot(legend=True, column='wave_dang', cmap=self.cmap, vmin=0, vmax=100, missing_kwds = {'color': 'black', 'label': 'Coast line'})
        # print(self.bbox_real)
        plt.annotate(\
            text='Wave angle: %s\nPrecision: %s' % (self.wave_spec['angle'], self.precision), \
            xy=(self.bbox_real_dict['xmin'], self.bbox_real_dict['ymax']),\
            verticalalignment='top'\
        )

        # city names
        if show_towns is True:
            for x, y, name in zip(towns.geometry.x, towns.geometry.y, towns.name):
                plt.annotate(name, xy=(x, y), xytext=(3, 3), textcoords="offset points")

        plt.title('Waves and the coastline intersection')
        plt.show()



# bbox = (-9.48859, 38.71225, -9.48369, 38.70596)
bbox = (-9.48859,38.70044,-9.4717541,38.7284016)
cascais = coast_part('/home/maksimpisarenko/tmp/osmcoast/coastlines-split-4326/lines.shp', bbox)
cascais.set_waves(angle=40)
cascais.set_wind()
cascais.ocean_plot(precision=0.001, show_bboxes=False, show_frames=True)