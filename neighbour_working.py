################################################################################
from qgis.utils import iface
from PyQt5.QtCore import QVariant # Updated for PyQt5
from qgis.core import QgsField, QgsFeature, QgsSpatialIndex

_NAME_FIELD = 'GEOID'
_SUM_FIELD = 'AWATER'

_NEW_NEIGHBORS_FIELD = 'NEIGHBORS'
_NEW_SUM_FIELD = 'SUM'

layer = iface.activeLayer()

layer.startEditing()
layer.dataProvider().addAttributes(
[QgsField(_NEW_NEIGHBORS_FIELD, QVariant.String),
QgsField(_NEW_SUM_FIELD, QVariant.Int)])
layer.updateFields()

feature_dict = {f.id(): f for f in layer.getFeatures()}

index = QgsSpatialIndex()
for f in feature_dict.values():
    index.insertFeature(f)

for f in feature_dict.values():
    print(f'Working on {f[_NAME_FIELD]}')
    geom = f.geometry()
    intersecting_ids = index.intersects(geom.boundingBox())
    
    neighbors = []
    neighbors_sum = 0
    for intersecting_id in intersecting_ids:
        intersecting_f = feature_dict[intersecting_id]

        if (f != intersecting_f and not intersecting_f.geometry().disjoint(geom)):
            neighbors.append(intersecting_f[_NAME_FIELD])
            neighbors_sum += intersecting_f[_SUM_FIELD]

    print(neighbors)
    print(neighbors_sum)
    f.setAttribute(_NEW_NEIGHBORS_FIELD, ','.join(neighbors))
    f.setAttribute(_NEW_SUM_FIELD, neighbors_sum)
    layer.updateFeature(f)
    
layer.commitChanges()
print('Processing complete.')