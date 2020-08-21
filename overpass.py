import overpy
import geopandas
from shapely.geometry import Point
import matplotlib.pyplot as plt

api = overpy.Overpass()

result = api.query("""
(
  node
  ["place"~"city|town|village|hamlet"]
    (38.6893216,-9.5005921,38.8869669,-9.3720427);
)->._;
(._;>;);
out;""")
# result = api.query("""
# (
#   node
#   ["place"~"city|town|village|hamlet"]
#     (38.725965379021, -9.525146484375, 38.856820134744, -9.411506652832);
# )->._;
# (._;>;);
# out;
#     """)
print(result)

towns_points_coord = []
towns_points_names = []
for node in result.nodes:
    print(node.tags['name'], node.lat, node.lon)
    towns_points_names.append(node.tags['name'])
    towns_points_coord.append(Point(node.lon,node.lat))

points = geopandas.GeoDataFrame(towns_points_names, geometry=towns_points_coord, columns=['name'])
# print(points)
points.plot()
for x, y, name in zip(points.geometry.x, points.geometry.y, points.name):
    plt.annotate(name, xy=(x, y), xytext=(3, 3), textcoords="offset points")
plt.show()