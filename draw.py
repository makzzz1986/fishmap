import geopandas
import pandas
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiPoint
import scipy.special as sc


def height(angle, bottom, max=9999):
    # tangens!
    height = bottom * sc.tandg(angle)
    print(bottom, height, max)
    if height > max:
        return max
    else:
        return height

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

wave_angle = 30
precision = 30

ysteps = (ymax - ymin) / precision
xsteps = (xmax - xmin) / precision

# draw wave lines
waves = geopandas.GeoDataFrame([], columns=['geometry'], crs="EPSG:4326")
for step in range(precision):
    waves.loc[len(waves), 'geometry'] = LineString(wave_line((xmin, ymin + (ysteps * step)), wave_angle, real_bbox))
    waves.loc[len(waves), 'geometry'] = LineString(wave_line((xmin + (xsteps * step), ymin), wave_angle, real_bbox))


# points of intersection
intersection_list = []
waves_ocean_list = []
for _, fid, coastline in coast.itertuples():
    for _, wave in waves.itertuples():
        intersect = coastline.intersection(wave)
        # removing not intersected:
        if not intersect.is_empty:
            intersection_list.append(intersect)
            # print(intersect)
            # drawing parted line
            if intersect.type == 'Point':
                waves_ocean_list.append(LineString([wave.coords[0], intersect]))
            if intersect.type == 'MultiPoint':
                line_list = []
                print(wave)
                line_list.append(wave.coords[0])
                line_list.extend([[point.x, point.y] for point in intersect])
                print(line_list)
                # if it is odd, than it ends on the ground
                if len(intersect) % 2 == 0:
                    line_list.append(wave.coords[-1]())
                # print(line_list)
                # print(len(line_list))
                
                # pair dots to lines!
                for pair in range(0, len(line_list), 2):
                    waves_ocean_list.append(LineString([line_list[pair], line_list[pair+1]]))
                # waves_ocean_list.append(MultiPoint(line_list))
waves_ocean_list.append(LineString(line_list))
      
intersection_points = geopandas.GeoDataFrame(geometry=intersection_list)
waves_ocean = geopandas.GeoDataFrame(geometry=waves_ocean_list)
# print(waves_ocean)

combined = geopandas.GeoDataFrame(pandas.concat([coast, intersection_points, waves_ocean], ignore_index=True))
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