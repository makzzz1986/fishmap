import overpy
import geopandas
import shapely.geometry as geometry
import matplotlib.pyplot as plt

api = overpy.Overpass()

# fetch all ways and nodes
result = api.query("""
(
  node
  ["place"~"city|town|village|hamlet"]
    (38.725965379021, -9.525146484375, 38.856820134744, -9.411506652832);
)->._;
(._;>;);
out;
    """)
print(result)

points_towns = []
for node in result.nodes:
    print(node.tags['name'], node.lat, node.lon)
    points_towns.append(geometry.Point(node.lon,node.lat))


points = geopandas.GeoDataFrame(points_towns, columns=['geometry'], crs="EPSG:4326")
print(points)
points.plot()
plt.show()