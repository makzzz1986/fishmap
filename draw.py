import geopandas
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString

bbox = (-9.48859, 38.71225, -9.48369, 38.70596)
coast = geopandas.read_file('/home/maksimpisarenko/tmp/osmcoast/coastlines-split-4326/lines.shp', bbox=bbox)

xmin = bbox[0]
xmax = bbox[2]
ymin = bbox[1]
ymax = bbox[3]

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

ten_buffer = (ymax - ymin) / 10
for ten in range(11):
    coast.loc[len(coast), 'geometry'] = LineString([(xmin, ymin + (ten_buffer * ten)), (xmax, ymin + (ten_buffer * ten))])
    print(ymin + (ten_buffer * ten))


print(xmin, xmax, ymin, ymax)
print(coast)

coast.plot()
plt.show()