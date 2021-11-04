# ArcGIS script to fill glacier outline attributes
# T Black, 29 July 2020

import arcpy
import attributes as atb

# Set parameters. Outlines are chosen by user in ArcGIS GUI.
outlines = arcpy.GetParameterAsText(0)
python_version = 'PYTHON_9.3'
field_names = [f.name for f in arcpy.ListFields(outlines)]

# Create a cursor object to loop through and update each row in the table.
with arcpy.da.UpdateCursor(outlines, field_names) as cursor:
    # For each row, get the Image_ID and use it to update each other field.
    for row in cursor:
        image_id = row[field_names.index('Image_ID')]

        row[field_names.index('Sensor')] = atb.getSensor(image_id)
        row[field_names.index('Image_Tile_Coordinate')] = \
            atb.getTileCoordinate(image_id)
        row[field_names.index('Image_Reference_System')] = \
            atb.getReferenceSystem(image_id)
        row[field_names.index('Source_Date')] = atb.getSourceDate(image_id)
        #row[field_names.index('Circadian_Date')] = \
        #    atb.getCircadianDate(image_id)
        row[field_names.index('Year_')] = atb.getYear(image_id)
        #row[field_names.index('Season')] = atb.getSeason(image_id)
              
        cursor.updateRow(row)

del cursor, row
