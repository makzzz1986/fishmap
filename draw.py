import geopandas
import pandas
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString
import scipy.special as sc
import numpy
from matplotlib.colors import ColorConverter, LinearSegmentedColormap


class coast_part():
    bbox = [] # xmin, ymin, xmax, ymax
    bbox_dict = {}
    bbox_geo = None
    bbox_real = [] # xmin, ymin, xmax, ymax
    bbox_real_dict = {}
    bbox_real_geo = None

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
    coastline_geo = None
    waves_geo = None
    ocean_geo = None
    cmap = None


    def __init__(self, file_path, bbox):
        self.bbox = bbox
        self.bbox_dict = self.bbox2dict(bbox)
        self.coastline_geo = geopandas.read_file(file_path, bbox=bbox)
        self.cmap = LinearSegmentedColormap.from_list("", ["green","yellow","red"])

        # Getting the real bbox! It is much bigger than bbox
        xmin = bbox[0]
        xmax = bbox[2]
        ymin = bbox[1]
        ymax = bbox[3]

        for _, n, r in self.coastline_geo.itertuples(): 
            for pair in list(r.coords):
                if pair[0] > xmax:
                    xmax = pair[0]
                if pair[0] < xmin:
                    xmin = pair[0]
                if pair[1] > ymax:
                    ymax = pair[1]
                if pair[1] < ymin:
                    ymin = pair[1]

        self.bbox_real = (xmin, ymin, xmax, ymax)
        self.bbox_real_dict = self.bbox2dict(self.bbox_real)

        # bboxes frames
        self.bbox_geo = geopandas.GeoDataFrame([], columns=['geometry', 'name'], crs="EPSG:4326")
        self.bbox_geo.loc[0] = {'name':'bbox_frame', 'geometry': MultiLineString([\
            ((self.bbox_dict['xmin'], self.bbox_dict['ymin']), (self.bbox_dict['xmax'], self.bbox_dict['ymin'])),\
            ((self.bbox_dict['xmax'], self.bbox_dict['ymin']), (self.bbox_dict['xmax'], self.bbox_dict['ymax'])),\
            ((self.bbox_dict['xmax'], self.bbox_dict['ymax']), (self.bbox_dict['xmin'], self.bbox_dict['ymax'])),\
            ((self.bbox_dict['xmin'], self.bbox_dict['ymax']), (self.bbox_dict['xmin'], self.bbox_dict['ymin']))\
        ])}

        self.bbox_real_geo = geopandas.GeoDataFrame([], columns=['geometry', 'name'] , crs='EPSG:4326')
        self.bbox_real_geo.loc[0] = {'name':'bbox_real_frame', 'geometry': MultiLineString([\
            ((self.bbox_real_dict['xmin'], self.bbox_real_dict['ymin']), (self.bbox_real_dict['xmax'], self.bbox_real_dict['ymin'])),\
            ((self.bbox_real_dict['xmax'], self.bbox_real_dict['ymin']), (self.bbox_real_dict['xmax'], self.bbox_real_dict['ymax'])),\
            ((self.bbox_real_dict['xmax'], self.bbox_real_dict['ymax']), (self.bbox_real_dict['xmin'], self.bbox_real_dict['ymax'])),\
            ((self.bbox_real_dict['xmin'], self.bbox_real_dict['ymax']), (self.bbox_real_dict['xmin'], self.bbox_real_dict['ymin']))\
        ])}
        

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
        end_point_y = ((xmax - start_point[0]) * sc.tandg(wave['angle'])) + start_point[1]
        if end_point_y > ymax:
            # print('Too big!')
            end_point_y = ymax
            end_point_x = ((ymax - start_point[1]) * sc.cotdg(wave['angle'])) + start_point[0]
        # print(start_point, (end_point_x, end_point_y), wave, bbox)
        return [start_point, (end_point_x, end_point_y)]

    def wave_parted(self, wave, intersect):
        # drawing parted line
        if intersect.type == 'Point':
            return {'waves': [LineString([wave.coords[0], intersect])], 'wave_dang': [self.wave_spec['dang']]}
        if intersect.type == 'MultiPoint':
            line_list = []
            # multilinestrings = []
            result_list = [] # colored lines
            line_list.append(wave.coords[0])
            line_list.extend([[point.x, point.y] for point in intersect])
            # if it is odd, than it ends on the ground
            if len(intersect) == 0:
                print('No intersect?')
                return {'waves': [], 'wave_dang': []}
            elif len(intersect) % 2 == 0:
                line_list.append(wave.coords[-1])
            # pair dots to lines!
            for pair in range(0, len(line_list), 2):
                result_list.append(LineString([line_list[pair], line_list[pair+1]]))
            # return list of wave LineStrings and their dangerously, wave_dang for the first one and 0 for others, 
            # because they are after ground.
            return {'waves': result_list, 'wave_dang': [0 if x>0 else self.wave_spec['dang'] for x in range(len(result_list))]} 

    def waves_set(self, angle=45, height=0, period=0):
        # Will be repaced by API call
        self.wave_spec = {
            'angle': angle, 
            'height': height,
            'period': period
        }
        # dangerousness should be calculated
        self.wave_spec['dang'] = 80

    def wind_set(self, angle=45, height=0, period=0):
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

    def intersection(self, geo1, geo2):
        # print('waves')
        # print(waves)
        # points of intersection
        # intersection_list = []
        intersected = {'waves': [], 'wave_dang': []}
        for _, _, geo1_line in geo1.itertuples():
            for _, geo2_line in geo2.itertuples():
                intersect = geo1_line.intersection(geo2_line)
                # removing not intersected:
                if not intersect.is_empty:
                    # intersection_list.append(intersect)
                    # drawing parted line
                    wave_parts = self.wave_parted(geo2_line, intersect)
                    intersected['waves'].extend(wave_parts['waves'])
                    intersected['wave_dang'].extend(wave_parts['wave_dang'])
            
        # intersection_points = geopandas.GeoDataFrame(geometry=intersection_list)  # points of intercestion wave and coastline
        intersection = geopandas.GeoDataFrame(intersected['wave_dang'], geometry=intersected['waves'], columns=['wave_dang'])

        return intersection

    def combination(self, geos):
        return geopandas.GeoDataFrame(pandas.concat(geos, ignore_index=True))


    def ocean_plot(self, precision=0.0001):
        self.precision = precision
        self.waves_geo = self.wave_draw(self.bbox_real_dict, self.wave_spec, self.precision)
        intersection = self.intersection(self.coastline_geo, self.waves_geo)
        self.ocean_geo = self.combination([self.coastline_geo, intersection, self.bbox_real_geo, self.bbox_geo])
        # print(self.ocean_geo)
        self.ocean_geo.plot(legend=True, column='wave_dang', cmap=self.cmap, vmin=0, vmax=100, missing_kwds = {'color': 'black', 'label': 'Coast line'})
        plt.annotate(\
            text='Wave angle: %s\nPrecision: %s' % (self.wave_spec['angle'], self.precision), \
            xy=(self.bbox_real_dict['xmin'], self.bbox_real_dict['ymax']),\
            verticalalignment='top'\
        )
        plt.title('Waves and the coastline intersection')
        plt.show()



bbox = (-9.48859, 38.71225, -9.48369, 38.70596)
cascais = coast_part('/home/maksimpisarenko/tmp/osmcoast/coastlines-split-4326/lines.shp', bbox)
cascais.waves_set(angle=50)
cascais.wind_set()
cascais.ocean_plot()