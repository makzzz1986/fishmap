import geopandas
import pandas
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString
import scipy.special as sc
from matplotlib.colors import ListedColormap


# will be updated, works for the first quarter only 8()
def wave_line(start_point, angle, bbox):
    xmax = 0
    ymax = 0
    if (0 < angle <= 90):
        xmax = bbox[2]
        ymax = bbox[3]
    elif (90 < angle <= 180):
        xmax = bbox[0]
        ymax = bbox[3]
    elif (180 < angle <= 270):
        xmax = bbox[0]
        ymax = bbox[1]
    elif (270 < angle <= 360):
        xmax = bbox[2]
        ymax = bbox[1]
    # print(start_point, xmax, ymax)
    end_point_x = xmax
    end_point_y = ((xmax - start_point[0]) * sc.tandg(angle)) + start_point[1]
    if end_point_y > ymax:
        # print('Too big!')
        end_point_y = ymax
        end_point_x = ((ymax - start_point[1]) * sc.cotdg(angle)) + start_point[0]
    # print(start_point, (end_point_x, end_point_y))
    return [start_point, (end_point_x, end_point_y)]


def wave_parted(wave, intersect, wave_dang=0):
    # drawing parted line
    if intersect.type == 'Point':
        return {'waves': [LineString([wave.coords[0], intersect])], 'wave_dang': [wave_dang]}
    if intersect.type == 'MultiPoint':
        line_list = []
        # multilinestrings = []
        result_list = [] # colored lines
        line_list.append(wave.coords[0])
        line_list.extend([[point.x, point.y] for point in intersect])
        # if it is odd, than it ends on the ground
        if len(intersect) % 2 == 0:
            line_list.append(wave.coords[-1]())
        # pair dots to lines!
        for pair in range(0, len(line_list), 2):
            result_list.append(LineString([line_list[pair], line_list[pair+1]]))
        # return list of wave LineStrings and their dangerously, wave_dang for the first one and 0 for others, 
        # because they are after ground.
        return {'waves': result_list, 'wave_dang': [0 if x>0 else wave_dang for x in range(len(result_list))]} 

bbox = (-9.48859, 38.71225, -9.48369, 38.70596)
coast = geopandas.read_file('/home/maksimpisarenko/tmp/osmcoast/coastlines-split-4326/lines.shp', bbox=bbox)

xmin = bbox[0]
xmax = bbox[2]
ymin = bbox[1]
ymax = bbox[3]

# Getting the real bbox! It is much bigger than bbox
for _, n, r in coast.itertuples(): 
    for pair in list(r.coords):
        if pair[0] > xmax:
            xmax = pair[0]
        if pair[0] < xmin:
            xmin = pair[0]
        if pair[1] > ymax:
            ymax = pair[1]
        if pair[1] < ymin:
            ymin = pair[1]

real_bbox = (xmin, ymin, xmax, ymax)

cmap = ListedColormap(['blue'], name='allblue')  # colormap
wave_angle = 30
wave_dang = 10
precision = 0.001

# draw wave lines
waves = geopandas.GeoDataFrame([], columns=['geometry'], crs="EPSG:4326")

xstep = xmin
ystep = ymin
while True:
    waves.loc[len(waves), 'geometry'] = LineString(wave_line((xstep, ymin), wave_angle, real_bbox))
    xstep += precision
    if xstep >= xmax:
        break
while True:
    waves.loc[len(waves), 'geometry'] = LineString(wave_line((xmin, ystep), wave_angle, real_bbox))
    ystep += precision
    if ystep >= ymax:
        break


# points of intersection
intersection_list = []
waves_intersected = {'waves': [], 'wave_dang': []}
for _, fid, coastline in coast.itertuples():
    for _, wave in waves.itertuples():
        intersect = coastline.intersection(wave)
        # removing not intersected:
        if not intersect.is_empty:
            intersection_list.append(intersect)
            # drawing parted line
            wave_parts = wave_parted(wave, intersect, wave_dang)
            print(wave_parts, len(wave_parts['waves']), len(wave_parts['wave_dang']))
            waves_intersected['waves'].extend(wave_parts['waves'])
            waves_intersected['wave_dang'].extend(wave_parts['wave_dang'])
      
print(waves_intersected)
print(len(waves_intersected['waves']))
print(len(waves_intersected['wave_dang']))

intersection_points = geopandas.GeoDataFrame(geometry=intersection_list)
waves = geopandas.GeoDataFrame(waves_intersected['wave_dang'], geometry=waves_intersected['waves'], columns=['wave_dang'])
print(waves)

combined = geopandas.GeoDataFrame(pandas.concat([coast, waves], ignore_index=True)).plot(cmap=cmap)
# coast.loc[len(coast), 'geometry'] = intersection

# print('bbox', bbox, 'real_bbox: ', (xmin, ymin), (xmax, ymax))

# add bboxes to plot
coast.loc[len(coast), 'geometry'] = LineString([(bbox[0], bbox[1]), (bbox[2], bbox[1])])
coast.loc[len(coast), 'geometry'] = LineString([(bbox[2], bbox[1]), (bbox[2], bbox[3])])
coast.loc[len(coast), 'geometry'] = LineString([(bbox[2], bbox[3]), (bbox[0], bbox[3])])
coast.loc[len(coast), 'geometry'] = LineString([(bbox[0], bbox[3]), (bbox[0], bbox[1])])

coast.loc[len(coast), 'geometry'] = LineString([(real_bbox[0], real_bbox[1]), (real_bbox[2], real_bbox[1])])
coast.loc[len(coast), 'geometry'] = LineString([(real_bbox[2], real_bbox[1]), (real_bbox[2], real_bbox[3])])
coast.loc[len(coast), 'geometry'] = LineString([(real_bbox[2], real_bbox[3]), (real_bbox[0], real_bbox[3])])
coast.loc[len(coast), 'geometry'] = LineString([(real_bbox[0], real_bbox[3]), (real_bbox[0], real_bbox[1])])


# print(coast)
# print(waves)

combined.plot()
plt.show()