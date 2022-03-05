#!/usr/bin/env python3
# For Kenai Fjords National Park - Coastal Glaciers project (2020)
# Import glacier geospatial information, calculate changes, and plot
# Taryn Black, August 2020

# %% Import modules
import geopandas as gpd
import pandas as pd
import manage
import glaciermetrics as gm
import matplotlib.pyplot as plt
import gplots as gp

# %% User-specified parameters
# Path for file geodatabase containing glacier information
fgdb = '../data/CoastalGlaciers.gdb'

# Path to save analysis output
outpath = '../output/'

# Path/filename to save area/length spreadsheet
data_spreadsheet = '../output/CoastalGlaciers_data.xlsx'

# File geodatabase layer names
outlines_layer = 'Glacier_Outlines'
reflines_layer = 'Glacier_Reference_Lines'
points_layer = 'Glacier_Points'
centerlines_layer = 'Glacier_Centerlines'
boxes_layer = 'Glacier_Boxes'

gp.designProperties()
fig_width = 3.38 # JOG specifies one-column figures <= 86mm, or 3.385 inches

# %% Load glacier layer information
outlines = gpd.read_file(fgdb, layer=outlines_layer, driver='FileGDB')

reflines = gpd.read_file(fgdb, layer=reflines_layer, driver='FileGDB')
reflines.sort_values(by='Glacier_ID', inplace=True)
reflines.set_index('Glacier_ID', drop=False, inplace=True)

points = gpd.read_file(fgdb, layer=points_layer, driver='FileGDB')
points.sort_values(by='Glacier_ID', inplace=True)
points.set_index('Glacier_ID', drop=False, inplace=True)

centerlines = gpd.read_file(fgdb, layer=centerlines_layer, driver='FileGDB')
centerlines.sort_values(by='Glacier_ID', inplace=True)
centerlines.set_index('Glacier_ID', drop=False, inplace=True)

refboxes = gpd.read_file(fgdb, layer=boxes_layer, driver='FileGDB')
refboxes.sort_values(by='Glacier_ID', inplace=True)
refboxes.set_index('Glacier_ID', drop=False, inplace=True)

# %% Get scope of glaciers and time
GIDS = points.Glacier_ID.values

outlines['Year'] = pd.to_datetime(outlines['Year_'], format='%Y')
del outlines['Year_']
YEAR_START = outlines.Year.min().year
YEAR_END = outlines.Year.max().year
YEARS = range(YEAR_START, YEAR_END+1)

DATE_START = outlines.Source_Date.min()
DATE_END = outlines.Source_Date.max()

START_DECADE = YEAR_START - YEAR_START%10
if YEAR_END%10 == 0:
    END_DECADE = YEAR_END
else:
    END_DECADE = YEAR_END + (10 - YEAR_END%10)
RANGE_DECADES = range(START_DECADE, END_DECADE, 10)
DECADES = [pd.to_datetime('{}-01-01'.format(y)).date() for y in RANGE_DECADES]

# %% Initialize dictionary of Glacier objects to store glacier information
all_glaciers = {id: manage.Glacier(id) for id in GIDS}
for id in all_glaciers:
    all_glaciers[id].refline = reflines.loc[id].geometry
    all_glaciers[id].refbox = refboxes.loc[id].geometry
    if centerlines.loc[id].ndim > 1:
        all_glaciers[id].centerline = centerlines.loc[id].geometry.iloc[0]
    else:
        all_glaciers[id].centerline = centerlines.loc[id].geometry
    all_glaciers[id].officialname = points.loc[id].Official_Name
    all_glaciers[id].unofficialname = points.loc[id].Unofficial_Name
    all_glaciers[id].fjordname = points.loc[id].Fjord_Name
    all_glaciers[id].armname = points.loc[id].Arm_Name

# %% Construct an observation time series for each glacier
for id in GIDS:
    
    # Get reference line and all observations for glacier ID
    glacier = outlines.query('Glacier_ID == @id')
    refline = all_glaciers[id].refline
    refbox  = all_glaciers[id].refbox
    cenline = all_glaciers[id].centerline
    
    # Loop through all observations and process data
    for n in range(len(glacier)):
        observation = glacier.iloc[n]

        # Create a Terminus Observation for a row in geodataframe
        obs = manage.TerminusObservation(gid=observation.Glacier_ID,
                                         qflag=observation.Quality_Flag,
                                         termination=observation.Termination_Type,
                                         imageid=observation.Image_ID,
                                         sensor=observation.Sensor,
                                         date=observation.Source_Date,
                                         circadiandate=observation.Circadian_Date,
                                         year=observation.Year,
                                         season=observation.Season,
                                         geometry=observation.geometry)
        
        # Calculate glacier area against reference line
        obs.area = gm.glacierArea(obs, refline)

        # Calculate glacier terminus area (excluding sides) against ref box
        obs.termarea = gm.glacierArea(obs, refbox)

        # Find outline intersection with centerline and resulting length
        inx_point, sublength = gm.centerlineIntersection(obs, cenline)
        obs.centerlineintersection = inx_point
        obs.length = sublength
        
        # Add glacier observation to time series
        all_glaciers[id].add_observation(obs)
    
    # Ensure that all observations are sorted by date
    all_glaciers[id].sort_by_date()

    # Extract time series of area, length, and dates from all observations
    all_glaciers[id].areas = all_glaciers[id].extract('area')
    all_glaciers[id].termareas = all_glaciers[id].extract('termarea')
    all_glaciers[id].lengths = all_glaciers[id].extract('length')
    all_glaciers[id].dates = all_glaciers[id].extract('date')

# %% For each glacier, determine termination type (for plotting groups)
lake_terminating = []
tidewater = []
land_terminating = []
mixed_terminating = []

for id in GIDS:
    termination_types = all_glaciers[id].extract('termination').unique()
    if len(termination_types) == 1:
        if termination_types == 'LAKE':
            lake_terminating.append(id)
        elif termination_types == 'TW':
            tidewater.append(id)
        elif termination_types == 'LAND':
            land_terminating.append(id)
    elif len(termination_types) > 1:
        mixed_terminating.append(id)

# %% For each glacier, create plots and calculate other metrics

for id in GIDS:
    glacier = all_glaciers[id]
    print('Analyzing glacier #{}: {}'.format(id, gp.getGlacierName(glacier)))

    # Plot relative area over time
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.totalRelativeMeasure(ax, glacier, 'area')
    ax.set_xticks(pd.to_datetime(['1985-01-01', '1990-01-01', '1995-01-01', '2000-01-01', '2005-01-01', '2010-01-01', '2015-01-01', '2020-01-01']))
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_relativearea.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)
    
    # Plot relative length over time
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.totalRelativeMeasure(ax, glacier, 'length')
    ax.set_xticks(pd.to_datetime(['1985-01-01', '1990-01-01', '1995-01-01', '2000-01-01', '2005-01-01', '2010-01-01', '2015-01-01', '2020-01-01']))
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_relativelength.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)
    
    # Plot relative area (with/without sides) and length over time, together
    fig, ax = plt.subplots(figsize=(fig_width, fig_width))
    gp.totalRelativeMeasureCompare(ax, glacier)
    ax.set_xticks(pd.to_datetime(['1985-01-01', '1990-01-01', '1995-01-01', '2000-01-01', '2005-01-01', '2010-01-01', '2015-01-01', '2020-01-01']))
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_sizecompare.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)

    # Plot spring and fall area separately
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.seasonRelativeMeasure(ax, glacier, 'area', spring=True, autumn=True)
    ax.set_xticks(pd.to_datetime(['1985-01-01', '1990-01-01', '1995-01-01', '2000-01-01', '2005-01-01', '2010-01-01', '2015-01-01', '2020-01-01']))
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_seasonarea.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)
    
    # Plot spring and fall length separately
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.seasonRelativeMeasure(ax, glacier, 'length', spring=True, autumn=True)
    ax.set_xticks(pd.to_datetime(['1985-01-01', '1990-01-01', '1995-01-01', '2000-01-01', '2005-01-01', '2010-01-01', '2015-01-01', '2020-01-01']))
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_seasonlength.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)

    # Plot seasonal area change between each measurement
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.seasonMeasureChange(ax, glacier, 'area', spring=True, autumn=True)
    ax.set_xticks(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_seasonareachange.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)

    # Plot seasonal length change between each measurement
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.seasonMeasureChange(ax, glacier, 'length', spring=True, autumn=True)
    ax.set_xticks(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_seasonlengthchange.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)

    # Plot difference between lengths, discontinuous across seasonal gaps
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(2*fig_width, fig_width))
    ax = ax.flatten()
    gp.lengthDiff(ax, glacier, YEARS)
    ax[0].set_xticks(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    ax[1].set_xticks(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    ax[0].set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    ax[1].set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_discontinuous_seasonlengthchange.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)

    # Plot decadal net area change
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.decadalMeasureChange(ax, glacier, 'area', DECADES)
    plt.savefig('{}{:02}_decadalareachange.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)

    # Plot decadal net length change
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.decadalMeasureChange(ax, glacier, 'length', DECADES)
    plt.savefig('{}{:02}_decadallengthchange.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300)

    # # Plot annual area change for past two decades
    # fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    # gp.individualMeasureChange(ax, glacier, 'area')
    # plt.savefig('{}{:02}_indivareachange.png'.format(outpath, id), \
    #     bbox_inches='tight', dpi=300)

    # # Plot annual length change for past two decades
    # fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    # gp.individualMeasureChange(ax, glacier, 'length')
    # plt.savefig('{}{:02}_indivlengthchange.png'.format(outpath, id), \
    #     bbox_inches='tight', dpi=300)

    # Plot rate of length change and rolling mean
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.meanChangeRate(ax, glacier, 'length', rolling_years=5)
    ax.set_xticks(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_lengthchangerate.png'.format(outpath, id), bbox_inches='tight', dpi=300)
    
    # Plot rate of area change and rolling mean
    fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
    gp.meanChangeRate(ax, glacier, 'area', rolling_years=5)
    ax.set_xticks(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    plt.savefig('{}{:02}_areachangerate.png'.format(outpath, id), bbox_inches='tight', dpi=300)

    # Plot relative length over time - in miles for StoryMap
    fig, ax = plt.subplots(figsize=(10,10))
    dates = glacier.extract('date') 
    cumul_length, _, _ = gm.netMeasureChange(glacier, 'length')
    cumul_length_mi = cumul_length / 1.609
    name = gp.getGlacierName(glacier)
    ax.plot(dates, cumul_length, 'o-', color='mediumblue', linewidth=4, markersize=10)
    ax.set_title(f'{name} Length Change', fontsize=24, fontweight='bold')
    ax.set_xlabel('Date', fontsize=20)
    ax.set_ylabel('Miles', fontsize=20)
    ax.set_xticks(pd.to_datetime(['1985-01-01', '1990-01-01', '1995-01-01', '2000-01-01', '2005-01-01', '2010-01-01', '2015-01-01', '2020-01-01']))
    ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
    ax.tick_params(axis='both', labelsize=16)
    plt.savefig('{}{:02}_storymap_length.png'.format(outpath, id), \
        bbox_inches='tight', dpi=300, facecolor='white', transparent=False)

    plt.close('all')
    
    # Calculate output metrics and save to file         
    with open('{}{:02}_metrics.txt'.format(outpath, id), 'w+') as f:
        name = gp.getGlacierName(glacier)
        f.write('ID #{}: {}\n'.format(id, name))
        f.write('{} total observations\n'.format(len(glacier.extract('area'))))
        f.write('Analyzed on: {}\n\n'.format(pd.Timestamp.utcnow().round('s')))

        # Calculate overall net area change and average change rate
        cumul_areachange, change_dates, num_obs = gm.netMeasureChange(
            glacier, 'area')
        cumul_termareachange, change_dates, num_obs = gm.netMeasureChange(
            glacier, 'termarea')
        areachange_rate_total = gm.rateMeasureChange(glacier, 'area')

        cumul_lenchange, _, _ = gm.netMeasureChange(glacier, 'length')
        lenchange_rate_total = gm.rateMeasureChange(glacier, 'length')

        # Calculate seasonal changes
        all_lengths = gm.fullTimeSeries(glacier, YEARS)
        aut_med = all_lengths[all_lengths.season == 'AUT'].ldiff.median(skipna=True)
        aut_num = all_lengths[all_lengths.season == 'AUT'].ldiff.count()
        spr_med = all_lengths[all_lengths.season == 'SPR'].ldiff.median(skipna=True)
        spr_num = all_lengths[all_lengths.season == 'SPR'].ldiff.count()

        f.write(f'Observed date range: {change_dates[0]} to {change_dates[1]}\n')
        f.write(f'Number of observations: {num_obs}\n')
        f.write(f'Net length change: {cumul_lenchange.iloc[-1]:.3f} km\n')
        f.write(f'Net area change: {cumul_areachange.iloc[-1]:.3f} km2\n')
        f.write(f'Net termarea change: {cumul_termareachange.iloc[-1]:.3f} km2\n')
        f.write(f'Rate of area change: {areachange_rate_total:.3f} km2/yr\n')
        f.write(f'Rate of length change: {lenchange_rate_total:.3f} km/yr\n')
        f.write(f'Median autumn difference (rel to spring): {aut_med*1000:.0f}m (n={aut_num})\n')
        f.write(f'Median spring difference (rel to autumn): {spr_med*1000:.0f}m (n={spr_num})\n\n')

        # Calculate net area change and average change rate for each decade
        for start_year in DECADES:
            end_year = gm.addDecade(start_year)
            cumul_areachange, change_dates, num_obs = gm.netMeasureChange(
                glacier, 'area', start_year, end_year)
            areachange_rate = gm.rateMeasureChange(
                glacier, 'area', start_year, end_year)

            cumul_lenchange, _, _ = gm.netMeasureChange(
                glacier, 'length', start_year, end_year)
            lenchange_rate = gm.rateMeasureChange(
                glacier, 'length', start_year, end_year)
            
            # Write outputs to text file
            f.write('Subset date range: {} to {}\n'.format(
                start_year, end_year))
            f.write('Observed date range: {} to {}\n'.format(
                change_dates[0], change_dates[1]))
            f.write('Number of observations: {}\n'.format(num_obs))
            f.write('Net area change: {:.3f} km2\n'.format(
                cumul_areachange.iloc[-1]))
            f.write('Net length change: {:.3f} km\n'.format(
                cumul_lenchange.iloc[-1]))
            f.write('Rate of area change: {:.3f} km2/yr\n'.format(
                areachange_rate))
            f.write('Rate of length change: {:.3f} km/yr\n\n'.format(
                lenchange_rate))
            
            plt.close('all')

# %% Create spreadsheet of glacier length/area data
with pd.ExcelWriter(data_spreadsheet, date_format='YYYY-MM-DD') as writer:
    for id in GIDS:
        glacier = all_glaciers[id]
        time_diff = pd.to_timedelta(glacier.dates.diff().values).days / 365
        area_changerate = glacier.areas.diff() / time_diff
        termarea_changerate = glacier.termareas.diff() / time_diff
        length_changerate = glacier.lengths.diff() / time_diff
        ss_data = {'area (km2)': glacier.areas,
                   'area change rate (km2/yr)': area_changerate,
                   'terminus area (km2)': glacier.termareas,
                   'terminus area change rate (km2/yr)': termarea_changerate,
                   'length (km)': glacier.lengths,
                   'length change rate (m/yr)': length_changerate * 1000,
                   'date': glacier.dates}
        ss_df = pd.DataFrame(ss_data)
        ss_df.to_excel(writer, \
            sheet_name='{} {}'.format(id, gp.getGlacierName(glacier)))

# %% Create summary plots of all glaciers

# Plot all observations per glacier over time (scatter plot)
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.glacierObservations(ax, all_glaciers)
plt.savefig('{}observation_timeseries.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot area change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'area', group=lake_terminating, \
    subtitle='\nLake-terminating Glaciers')
plt.savefig('{}summary_lake_areachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot area change for land-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'area', group=land_terminating, \
    subtitle='\nLand-terminating Glaciers')
plt.savefig('{}summary_land_areachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot area change for mixed-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'area', group=mixed_terminating, \
    subtitle='\nMixed-terminating Glaciers')
plt.savefig('{}summary_mixed_areachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot area change for tidewater glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'area', group=tidewater, \
    subtitle='\nTidewater Glaciers')
plt.savefig('{}summary_tidewater_areachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot terminus area change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'termarea', group=lake_terminating, \
    subtitle='\nLake-terminating Glaciers')
plt.savefig('{}summary_lake_termareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot terminus area change for land-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'termarea', group=land_terminating, \
    subtitle='\nLand-terminating Glaciers')
plt.savefig('{}summary_land_termareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot terminus area change for mixed-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'termarea', group=mixed_terminating, \
    subtitle='\nMixed-terminating Glaciers')
plt.savefig('{}summary_mixed_termareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot terminus area change for tidewater glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'termarea', group=tidewater, \
    subtitle='\nTidewater Glaciers')
plt.savefig('{}summary_tidewater_termareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot length change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'length', group=lake_terminating, \
    subtitle='\nLake-terminating Glaciers')
plt.savefig('{}summary_lake_lengthchange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot length change for land-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'length', group=land_terminating, \
    subtitle='\nLand-terminating Glaciers')
plt.savefig('{}summary_land_lengthchange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot length change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'length', group=mixed_terminating, \
    subtitle='\nMixed-terminating Glaciers')
plt.savefig('{}summary_mixed_lengthchange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot length change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.measureSummary(ax, all_glaciers, 'length', group=tidewater, \
    subtitle='\nTidewater Glaciers')
plt.savefig('{}summary_tidewater_lengthchange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized area change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'area', group=lake_terminating, \
    subtitle='\nLake-terminating Glaciers')
plt.savefig('{}summary_lake_normareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized area change for land-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'area', group=land_terminating, \
    subtitle='\nLand-terminating Glaciers')
plt.savefig('{}summary_land_normareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized area change for mixed-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'area', group=mixed_terminating, \
    subtitle='\nMixed-terminating Glaciers')
plt.savefig('{}summary_mixed_normareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized area change for tidewater glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'area', group=tidewater, \
    subtitle='\nTidewater Glaciers')
plt.savefig('{}summary_tidewater_normareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized terminus area change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'termarea', group=lake_terminating, \
    subtitle='\nLake-terminating Glaciers')
plt.savefig('{}summary_lake_normtermareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized terminus area change for land-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'termarea', group=land_terminating, \
    subtitle='\nLand-terminating Glaciers')
plt.savefig('{}summary_land_normtermareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized terminus area change for mixed-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'termarea', group=mixed_terminating, \
    subtitle='\nMixed-terminating Glaciers')
plt.savefig('{}summary_mixed_normtermareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized terminus area change for tidewater glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'termarea', group=tidewater, \
    subtitle='\nTidewater Glaciers')
plt.savefig('{}summary_tidewater_normtermareachange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized length change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'length', group=lake_terminating, \
    subtitle='\nLake-terminating Glaciers')
plt.savefig('{}summary_lake_normlengthchange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized length change for land-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'length', group=land_terminating, \
    subtitle='\nLand-terminating Glaciers')
plt.savefig('{}summary_land_normlengthchange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized length change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'length', \
    group=mixed_terminating, subtitle='\nMixed-terminating Glaciers')
plt.savefig('{}summary_mixed_normlengthchange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot normalized length change for lake-terminating glaciers
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
gp.normMeasureSummary(ax, all_glaciers, 'length', group=tidewater, \
    subtitle='\nTidewater Glaciers')
plt.savefig('{}summary_tidewater_normlengthchange.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)


#%% PUBLICATION FIGURES

# Plot all observations per glacier over time (scatter plot)
fig, ax = plt.subplots(figsize=(2*fig_width, 0.85*fig_width))
yticklabs = []
for g in all_glaciers:
    glacier = all_glaciers[len(all_glaciers)+1-g]
    name = gp.getGlacierName(glacier)
    yticklabs.append(name)
    glacierid = glacier.extract('gid')
    obs_dates = glacier.extract('date')
    seasons = glacier.extract('season')
    ax.plot(obs_dates[seasons=='AUT'], len(all_glaciers)+1-glacierid[seasons=='AUT'], '.', color='darkorange', label='Autumn')
    ax.plot(obs_dates[seasons=='SPR'], len(all_glaciers)+1-glacierid[seasons=='SPR'], '.', color='mediumseagreen', label='Spring')
ax.set_xlabel('Date')
ax.legend(handles=ax.get_lines()[0:2], ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.18))
yticklocs = range(1, len(all_glaciers)+1)
plt.yticks(yticklocs, yticklabs)
plt.ylim(bottom=0)

ax.set_xticks(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
ax.set_xticklabels(['1985', '1990', '1995', '2000', '2005', '2010', '2015', '2020'])
ax.set_ylabel(ax.get_ylabel(), visible=False)
ax.set_title(ax.get_title(), visible=False)
plt.savefig('{}fig3_observation_timeseries.png'.format(outpath), \
    bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot lake-terminating area, length, decadal area, decadal length
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(2*fig_width, fig_width))
gp.measureSummary(ax1, all_glaciers, 'area', group=lake_terminating, \
    subtitle='\nLake-terminating Glaciers')
gp.measureSummary(ax2, all_glaciers, 'length', group=lake_terminating, \
    subtitle='\nLake-terminating Glaciers')
ax1.legend().remove()
ax2.legend().remove()
fig.legend(handles=ax1.get_lines(), ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.02))
ax1.set_xlim(left=pd.to_datetime('1980-01-01'), right=pd.to_datetime('2025-01-01'))
ax1.set_xticks(pd.to_datetime(['1980-01-01', '1990-01-01', '2000-01-01', '2010-01-01', '2020-01-01']))
ax1.set_xticklabels(['1980', '1990', '2000', '2010', '2020'])
ax1.set_title(ax1.get_title(), visible=False)
ax2.set_title(ax2.get_title(), visible=False)
ax1.annotate(text='a', xy=(-0.15, 1.05), xycoords='axes fraction')
ax2.annotate(text='b', xy=(-0.15, 1.05), xycoords='axes fraction')
plt.tight_layout()
plt.savefig('{}fig6_summary_subplots_lake.png'.format(outpath), bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot land-terminating area, length, decadal area, decadal length
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(2*fig_width, fig_width))
gp.measureSummary(ax1, all_glaciers, 'area', group=land_terminating, \
    subtitle='\nLand-terminating Glaciers')
gp.measureSummary(ax2, all_glaciers, 'length', group=land_terminating, \
    subtitle='\nLand-terminating Glaciers')
ax1.legend().remove()
ax2.legend().remove()
fig.legend(handles=ax1.get_lines(), ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.02))
ax1.set_xlim(left=pd.to_datetime('1980-01-01'), right=pd.to_datetime('2025-01-01'))
ax1.set_xticks(pd.to_datetime(['1980-01-01', '1990-01-01', '2000-01-01', '2010-01-01', '2020-01-01']))
ax1.set_xticklabels(['1980', '1990', '2000', '2010', '2020'])
ax1.set_title(ax1.get_title(), visible=False)
ax2.set_title(ax2.get_title(), visible=False)
ax1.annotate(text='a', xy=(-0.15, 1.05), xycoords='axes fraction')
ax2.annotate(text='b', xy=(-0.15, 1.05), xycoords='axes fraction')
plt.tight_layout()
plt.savefig('{}fig7_summary_subplots_land.png'.format(outpath), bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot mixed-terminating area, length, decadal area, decadal length
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(2*fig_width, fig_width))
gp.measureSummary(ax1, all_glaciers, 'area', group=mixed_terminating, \
    subtitle='\nMixed-terminating Glaciers')
gp.measureSummary(ax2, all_glaciers, 'length', group=mixed_terminating, \
    subtitle='\nMixed-terminating Glaciers')
ax1.legend().remove()
ax2.legend().remove()
fig.legend(handles=ax1.get_lines(), ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.02))
ax1.set_xlim(left=pd.to_datetime('1980-01-01'), right=pd.to_datetime('2025-01-01'))
ax1.set_xticks(pd.to_datetime(['1980-01-01', '1990-01-01', '2000-01-01', '2010-01-01', '2020-01-01']))
ax1.set_xticklabels(['1980', '1990', '2000', '2010', '2020'])
ax1.set_title(ax1.get_title(), visible=False)
ax2.set_title(ax2.get_title(), visible=False)
ax1.annotate(text='a', xy=(-0.15, 1.05), xycoords='axes fraction')
ax2.annotate(text='b', xy=(-0.15, 1.05), xycoords='axes fraction')
plt.tight_layout()
plt.savefig('{}fig8_summary_subplots_mixed.png'.format(outpath), bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot tidewater area, length, decadal area, decadal length
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(2*fig_width, fig_width))
gp.measureSummary(ax1, all_glaciers, 'area', group=tidewater, \
    subtitle='\nTidewater Glaciers')
gp.measureSummary(ax2, all_glaciers, 'length', group=tidewater, \
    subtitle='\nTidewater Glaciers')
ax1.legend().remove()
ax2.legend().remove()
fig.legend(handles=ax1.get_lines(), ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.02))
ax1.set_xlim(left=pd.to_datetime('1980-01-01'), right=pd.to_datetime('2025-01-01'))
ax1.set_xticks(pd.to_datetime(['1980-01-01', '1990-01-01', '2000-01-01', '2010-01-01', '2020-01-01']))
ax1.set_xticklabels(['1980', '1990', '2000', '2010', '2020'])
ax1.set_title(ax1.get_title(), visible=False)
ax2.set_title(ax2.get_title(), visible=False)
ax1.annotate(text='a', xy=(-0.15, 1.05), xycoords='axes fraction')
ax2.annotate(text='b', xy=(-0.15, 1.05), xycoords='axes fraction')
plt.tight_layout()
plt.savefig('{}fig5_summary_subplots_tidewater.png'.format(outpath), bbox_inches='tight', dpi=300)
plt.close(fig)

# Plot median seasonal differences
fig, ax = plt.subplots(figsize=(fig_width, 0.625*fig_width))
for id in GIDS:
    glacier = all_glaciers[id]
    all_lengths = gm.fullTimeSeries(glacier, YEARS)
    aut_avg = all_lengths[all_lengths.season == 'AUT'].ldiff.median(skipna=True)
    spr_avg = all_lengths[all_lengths.season == 'SPR'].ldiff.median(skipna=True)
    ann_var = spr_avg - aut_avg
    if id in lake_terminating:
        color = 'royalblue'
        marker = 'o'
        label = 'Lake-terminating'
    elif id in tidewater:
        color = 'mediumseagreen'
        marker = 's'
        label = 'Tidewater'
    elif id in land_terminating:
        color = 'darkorange'
        marker = '^'
        label = 'Land-terminating'
    elif id in mixed_terminating:
        color = 'darkorchid'
        marker = 'D'
        label = 'Mixed-terminating'
    ax.scatter(spr_avg*1000, aut_avg*1000, color=color, marker=marker, label=label)
# add shading to indicate image resolution limits
ax.axhspan(ymin=-15, ymax=15, facecolor='silver', alpha=0.5, zorder=0.75)
ax.axvspan(xmin=-15, xmax=15, facecolor='silver', alpha=0.5, zorder=0.75)
ax.set_xlabel('Median winter length change (m)')
ax.set_ylabel('Median summer length change (m)')
ax.axhline(0, color='darkgray', zorder=0.75)
ax.axvline(0, color='darkgray', zorder=0.75)
# custom legend
handles, labels = ax.get_legend_handles_labels()
first_ids = sorted([lake_terminating[0]-1, tidewater[0]-1, land_terminating[0]-1, mixed_terminating[0]-1])
first_ids
unique_handles = [handles[i] for i in first_ids]
unique_labels = [labels[i] for i in first_ids]
ax.legend(unique_handles, unique_labels, loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(f'{outpath}fig9_median_seasonal_variability.png', bbox_inches='tight', dpi=300)


print('Done.')

