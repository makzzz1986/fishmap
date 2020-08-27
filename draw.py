import geopandas
import overpy
import pandas
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString
from scipy.special import tandg, cotdg
import numpy
from matplotlib.colors import ColorConverter, LinearSegmentedColormap


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
        move_coords = [\
            str(self.ymin), \
            str(self.xmin), \
            str(self.ymax), \
            str(self.xmax)  \
        ]
        self.osm_coords = ','.join(move_coords)

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

    def frame_draw(self, bbox_dict, name):
        temp_geodataframe = geopandas.GeoDataFrame([], columns=['geometry', 'name'] , crs='EPSG:4326')
        temp_geodataframe.loc[0] = {'name': name, 'geometry': MultiLineString([\
            ((bbox_dict['xmin'], bbox_dict['ymin']), (bbox_dict['xmax'], bbox_dict['ymin'])),\
            ((bbox_dict['xmax'], bbox_dict['ymin']), (bbox_dict['xmax'], bbox_dict['ymax'])),\
            ((bbox_dict['xmax'], bbox_dict['ymax']), (bbox_dict['xmin'], bbox_dict['ymax'])),\
            ((bbox_dict['xmin'], bbox_dict['ymax']), (bbox_dict['xmin'], bbox_dict['ymin']))\
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
        self.coastline_geo = geopandas.read_file(file_path, bbox=bbox)
        self.geo_all.append(self.coastline_geo)
        self.cmap = LinearSegmentedColormap.from_list("", ["green","yellow","red"])

        # calculate the frames

        # the real frame which is much bigger than which we requested (I don't know the reason)
        xmin = self.bbox.xmax
        xmax = self.bbox.xmin
        ymin = self.bbox.ymax
        ymax = self.bbox.ymin
        for _, fid, r in self.coastline_geo.itertuples():
            # for every FID from the shape file we create their own frame
            xmin_frame = self.bbox.xmax
            xmax_frame = self.bbox.xmin
            ymin_frame = self.bbox.ymax
            ymax_frame = self.bbox.ymin
            print('FID of the object from shapefile:', fid)
            # find left bottom and right upper points
            for pair in list(r.coords):
                if pair[0] > xmax_frame:
                    xmax_frame = pair[0]
                if pair[0] < xmin_frame:
                    xmin_frame = pair[0]
                if pair[1] > ymax_frame:
                    ymax_frame = pair[1]
                if pair[1] < ymin_frame:
                    ymin_frame = pair[1]
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
            if xmax_frame > xmax:
                xmax = xmax_frame
            if xmin_frame < xmin:
                xmin = xmin_frame
            if ymax_frame > ymax:
                ymax = ymax_frame
            if ymin_frame < ymin:
                ymin = ymin_frame

        # remove duplicates
        self.frame_clusters = set(self.frame_clusters)
        print(self.frame_clusters)
        self.bbox_real = bbox_box((xmin, ymin, xmax, ymax), 'bbox_real')


    def check_incapsulation(self, bbox_1, bbox_2):
        # bbox_1 encapsulates bbox_2
        if  (bbox_1.xmin < bbox_2.xmin) and \
            (bbox_1.ymin < bbox_2.ymin) and \
            (bbox_1.xmax > bbox_2.xmax) and \
            (bbox_1.ymax > bbox_2.ymax):
            print(bbox_1, '>', bbox_2)
            return bbox_1
        # bbox_2 encapsulates bbox_1
        elif (bbox_2.xmin < bbox_1.xmin) and \
             (bbox_2.ymin < bbox_1.ymin) and \
             (bbox_2.xmax > bbox_1.xmax) and \
             (bbox_2.ymax > bbox_1.ymax):
            print(bbox_2, '>', bbox_1)
            return bbox_2
        else:
            print(bbox_1, '~', bbox_2)


    # will be updated, works for the first quarter only 8()
    def wave_line(self, start_point, wave, bbox):
        xmax = 0
        ymax = 0
        if (0 < wave['angle'] < 90):
            xmax = bbox.xmax # remake them to bbox_dict
            ymax = bbox.ymax
        elif (90 < wave['angle'] <= 180):
            xmax = bbox.xmin
            ymax = bbox.ymax
        elif (180 < wave['angle'] < 270):
            xmax = bbox.xmin
            ymax = bbox.ymin
        elif (270 < wave['angle'] <= 360):
            xmax = bbox.xmax
            ymax = bbox.ymin
        elif (wave['angle'] == 90) or (wave['angle'] == 270):
            # print([start_point, (start_point[0], bbox['ymax'])])
            return [start_point, (start_point[0], bbox.ymax)]
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
                # print('Odd?')
                # for i in intersect:
                #     print(i)
                # print(points_list)
                points_list.append(wave.coords[-1])
                # print(points_list)
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

        xstep = bbox.xmin
        ystep = bbox.ymin
        while True:
            waves_geo.loc[len(waves_geo), 'geometry'] = LineString(self.wave_line((xstep, bbox.ymin), wave_spec, bbox))
            xstep += precision
            if xstep >= bbox.xmax:
                break
        while True:
            waves_geo.loc[len(waves_geo), 'geometry'] = LineString(self.wave_line((bbox.xmin, ystep), wave_spec, bbox))
            ystep += precision
            if ystep >= bbox.ymax:
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
        for node in result.nodes:
            # print(node.tags['name'], node.lat, node.lon)
            towns_points_names.append(node.tags['name'])
            towns_points_coord.append(Point(node.lon,node.lat))
        return geopandas.GeoDataFrame(towns_points_names, geometry=towns_points_coord, columns=['name'])


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


        self.waves_geo = self.wave_draw(self.bbox_broadened, self.wave_spec, self.precision)
        
        # waves_parted = self.intersection(self.coastline_geo, self.waves_geo)
        waves_parted = self.intersection(self.waves_geo, self.coastline_geo)
        self.geo_all.append(waves_parted)

        if show_bboxes is True:
            self.geo_all.extend([self.bbox_real.geo, self.bbox.geo])

        if show_towns is True:
            towns = self.set_towns(self.bbox_real)
            self.geo_all.append(towns)

        if show_frames is True:
            # print([self.frame_draw(self.frame_fids[frame], frame) for frame in self.frame_fids])
            self.geo_all.extend([self.frame_fids[frame].geo for frame in self.frame_fids])

        self.ocean_geo = self.combination(self.geo_all)
        # print(self.coastline_geo)
        self.ocean_geo.plot(legend=True, column='wave_dang', cmap=self.cmap, vmin=0, vmax=100, missing_kwds = {'color': 'black', 'label': 'Coast line'})
        # print(self.bbox_real)
        plt.annotate(\
            text='Wave angle: %s\nPrecision: %s' % (self.wave_spec['angle'], self.precision), \
            xy=(self.bbox_real.xmin, self.bbox_real.ymax),\
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
cascais.ocean_plot(precision=0.01, show_towns=False, show_bboxes=False, show_frames=True)