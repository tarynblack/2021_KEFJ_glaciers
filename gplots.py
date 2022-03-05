# Functions for common plots for glacier terminus/area analysis

import matplotlib as mpl
import matplotlib.pyplot as plt
import glaciermetrics as gm
import pandas as pd
import numpy as np


season_c = {'spring': 'mediumseagreen',
            'summer': 'darkorchid',
            'autumn': 'darkorange',
            'winter': 'dodgerblue'}

measure_units = {'length': 'km', \
                 'area': '$km^2$',
                 'termarea': '$km^2$'}

# Pick a color scheme for each termination type, and for each color use different marker shapes for each glacier
# lake = dodgerblue
# land = darkorchid
# tidewater = mediumseagreen
# mixed = darkorange
# CHANGING THIS -- PLOT ALL OF ONE TERMINATION TYPE ON SINGLE PLOT, MIX COLORS FOR DISTINCTION
glacier_design = {
    1  : {'c' : 'dodgerblue',\
          's' : 'o-'},              # Bear, lake-terminating
    2  : {'c' : 'dodgerblue',\
          's' : 'o-'},              # Aialik, tidewater
    3  : {'c' : 'darkorchid',\
          's' : '^-'},              # Pedersen, lake-terminating
    4  : {'c' : 'darkorchid',\
          's' : '^-'},              # Holgate, tidewater
    5  : {'c' : 'dodgerblue',\
          's' : 'o-'},              # South Holgate - West, mixed
    6  : {'c' : 'darkorchid',\
          's' : '^-'},              # South Holgate - East, mixed
    7  : {'c' : 'dodgerblue',\
          's' : 'o-'},              # Northeastern, land-terminating
    8  : {'c' : 'mediumseagreen',\
          's' : 's-'},              # Northwestern, tidewater
    9  : {'c' : 'darkorange',\
          's' : 'D-'},              # Ogive, tidewater
    10 : {'c' : 'gray',\
          's' : 'v-'},              # Anchor, tidewater
    11 : {'c' : 'mediumseagreen',\
          's' : 's-'},              # Reconstitution, mixed
    12 : {'c' : 'darkorange',\
          's' : 'D-'},              # Southwestern, mixed
    13 : {'c' : 'darkorchid',\
          's' : '^-'},              # Sunrise, land-terminating 
    14 : {'c' : 'mediumseagreen',\
          's' : 's-'},              # Paguna, land-terminating
    15 : {'c' : 'deeppink',\
          's' : 'p-'},              # McCarty, tidewater
    16 : {'c' : 'darkorange',\
          's' : 'D-'},              # Dinglestadt, land-terminating
    17 : {'c' : 'gray',\
          's' : 'v-'},              # Split, land-terminating
    18 : {'c' : 'mediumseagreen',\
          's' : 's-'},              # Yalik, lake-terminating
    19 : {'c' : 'gray',\
          's' : 'v-'}               # Petrof, mixed
}

def designProperties():
    # set text properties
    mpl.rcParams['axes.titlesize'] = 10
    mpl.rcParams['axes.labelsize'] = 10
    mpl.rcParams['xtick.labelsize'] = 9
    mpl.rcParams['ytick.labelsize'] = 9
    mpl.rcParams['legend.fontsize'] = 9
    mpl.rcParams['legend.title_fontsize'] = 9
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'optima'

    # set line and marker properties
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['lines.markersize'] = 3

    # set axes and grid properties
    mpl.rcParams['axes.axisbelow'] = True
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['grid.color'] = 'lightgray'

    # set mathematical notation to be regular, not italicized
    mpl.rcParams['mathtext.default'] = 'regular'


def zeroLine(ax):
    ax.axhline(linewidth=0.75, color='darkgray', zorder=0.75)


def manageSubplots(fig, subplot_bool, idx):
    if subplot_bool == False:
        ax = fig.add_axes([0., 0., 1., 1.])
    elif subplot_bool == True:
        ax = fig.add_subplot(idx)
    return ax


def getGlacierName(GlacierClass):
    if GlacierClass.officialname:
        name = GlacierClass.officialname
    else:
        name = GlacierClass.unofficialname
    return name


def getSeasonMeasures(dates, measures, seasons, seasonstr):
    dates_season = dates.where(seasons==seasonstr).dropna()
    measures_season = measures.where(seasons==seasonstr).dropna()
    return dates_season, measures_season


def check_measure(measure):
    measure_types = ['length', 'area', 'termarea']
    if measure not in measure_types:
        raise ValueError("Invalid measure type. Expected one of: {}".format(
            measure_types))


def align_yscale(ax1, ax2):
    ax1min, ax1max = ax1.get_ylim()
    ax2min, ax2max = ax2.get_ylim()
    axmin = min(ax1min, ax2min)
    axmax = max(ax1max, ax2max)
    ax1.set_ylim(axmin, axmax)
    ax2.set_ylim(axmin, axmax)


def annotate_bars(ax, anno, x, y):
    for i in range(len(anno)):
        sign = y[i]/abs(y[i])
        ax.annotate(anno[i], 
                    xy=(x[i], 0),#y[i]), 
                    xytext = (0, -sign * 8),
                    textcoords = 'offset points',
                    ha='center', va='center')


def totalRelativeMeasure(ax, glacier, measure, subplots=False, idx=111):
    check_measure(measure)

    # ax = manageSubplots(fig, subplots, idx)
    dates = glacier.extract('date') 
    cumul_measures, _, _ = gm.netMeasureChange(glacier, measure)
    name = getGlacierName(glacier)

    graph, = ax.plot(dates, cumul_measures, 'o-', color='mediumblue')
    ax.set_title('{}: Observed {} Change'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('Cumulative {} Change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    # figureProperties(fig, ax, graph)


def totalRelativeMeasureCompare(ax, glacier, subplots=False, idx=111):
    dates = glacier.extract('date') 
    name = getGlacierName(glacier)
    cumul_area, _, _ = gm.netMeasureChange(glacier, 'area')
    cumul_termarea, _, _ = gm.netMeasureChange(glacier, 'termarea')
    cumul_length, _, _ = gm.netMeasureChange(glacier, 'length')
    
    # Plot areas on one axis
    # ax_area = manageSubplots(fig, subplots, idx)
    gr_area, = ax.plot(dates, cumul_area, 'o-', color='darkblue', \
        label='Terminus and lateral area')
    gr_termarea, = ax.plot(dates, cumul_termarea, 'o-', color='blue', \
        fillstyle='none', label='Terminus area only')
    # figureProperties(fig, ax_area, gr_termarea)
    gr_length, = ax.plot([], [], 's-', color='salmon', \
        label='Centerline length')
    ax.set_title('{}: Observed Size Changes'.format(name))
    ax.set_xlabel('Date')
    ax.set_ylabel('Cumulative Area Change ({})'.format(
        measure_units['area']))
    ax.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    
    # Plot length on other axis
    ax_length = ax.twinx()
    gr_length, = ax_length.plot(dates, cumul_length, 's-', color='salmon', \
        alpha=0.7, label='Centerline length')
    ax_length.set_ylabel('Cumulative Length Change ({})'.format(
        measure_units['length']))
    
    # figureProperties(fig, ax_area, gr_area)
    # figureProperties(fig, ax_area, gr_termarea)
    # figureProperties(fig, ax_length, gr_length)
    ax_length.grid(b=None)
    align_yscale(ax, ax_length)


def seasonRelativeMeasure(ax, glacier, measure, \
    spring=False, summer=False, autumn=False, winter=False, \
    subplots=False, idx=111):
    check_measure(measure)

    # ax = manageSubplots(fig, subplots, idx)
    dates = glacier.extract('date')
    measures = glacier.extract(measure)
    seasons = glacier.extract('season')
    name = getGlacierName(glacier)

    if spring:
        dates_spr, measures_spr = getSeasonMeasures(dates, measures, \
            seasons, 'SPR')
        cumul_measures_spr = measures_spr.diff().cumsum()
        cumul_measures_spr.iloc[0] = 0.0
        graph, = ax.plot(dates_spr, cumul_measures_spr, \
            'o-', color=season_c['spring'], label='Spring')
        # figureProperties(fig, ax, graph)
    if summer:
        dates_sum, measures_sum = getSeasonMeasures(dates, measures, \
            seasons, 'SUM')
        cumul_measures_sum = measures_sum.diff().cumsum()
        cumul_measures_sum.iloc[0] = 0.0
        graph, = ax.plot(dates_sum, cumul_measures_sum, \
            '^-', color=season_c['summer'], label='Summer')
        # figureProperties(fig, ax, graph)
    if autumn:
        dates_aut, measures_aut = getSeasonMeasures(dates, measures, \
            seasons, 'AUT')
        cumul_measures_aut = measures_aut.diff().cumsum()
        cumul_measures_aut.iloc[0] = 0.0
        graph, = ax.plot(dates_aut, cumul_measures_aut, \
            's-', color=season_c['autumn'], label='Autumn')
        # figureProperties(fig, ax, graph)
    if winter:
        dates_win, measures_win = getSeasonMeasures(dates, measures, \
            seasons, 'WIN')
        cumul_measures_win = measures_win.diff().cumsum()
        cumul_measures_win.iloc[0] = 0.0
        graph, = ax.plot(dates_win, cumul_measures_win, \
            'D-', color=season_c['winter'], label='Winter')
        # figureProperties(fig, ax, graph)
    
    ax.set_title('{}: Observed Seasonal {} Change'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('Cumulative {} Change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    ax.legend()
    # figureProperties(fig, ax, graph)


def individualMeasureChange(ax, glacier, measure, \
    date_start=None, date_end=None, subplots=False, idx=111):
    check_measure(measure)
    
    measures, dates = gm.filterByDates(glacier, measure, date_start, date_end)
    measure_change = measures.diff()
    name = getGlacierName(glacier)

    ax.plot(dates, measure_change, color='darkgray')
    ax.set_title('{}: {} Change Between Observations'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('{} change ({})'.format(
        measure.capitalize(), measure_units[measure]))


def decadalMeasureChange(ax, glacier, measure, decade_startyears, \
    subplots=False, idx=111):
    check_measure(measure)
    
    # ax = manageSubplots(fig, subplots, idx)
    name = getGlacierName(glacier)

    decade_startyears = pd.to_datetime(decade_startyears)
    decadal_net_changes = pd.Series(dtype='float64')
    decade_labels = []
    bar_annotations = []
    for startyear in decade_startyears:
        endyear = gm.addDecade(startyear)
        cumul_decadal_measure_change, _, num_obs = gm.netMeasureChange(
            glacier, measure, startyear,endyear)
        net_decadal_measure_change = cumul_decadal_measure_change.iloc[-1]
        midyear = endyear.year - 4
        decadal_net_changes.loc[midyear] = net_decadal_measure_change
        decade_labels.append('{}-{}'.format(startyear.year, startyear.year+9))
        bar_annotations.append('{} obsv'.format(num_obs))
    
    rects = ax.bar(decadal_net_changes.index.values, 
                   decadal_net_changes.values, 
                   color='mediumblue', width=5)
    annotate_bars(ax, bar_annotations, 
                  decadal_net_changes.index.values, decadal_net_changes.values)
    ax.set_title('{}: Observed Decadal {} Change'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Decade')
    ax.set_ylabel('Net {} Change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    ax.set_xticks([y.year+5 for y in decade_startyears])
    ax.set_xticklabels(decade_labels)
    if max(ax.get_ylim()) == 0.0:
        yrange, _ = ax.get_ylim()
        ax.set_ylim(top=abs(0.05*yrange))
    # figureProperties(fig, ax, rects)
    # ax.set_axisbelow(True)
    ax.grid(axis='x')


def seasonMeasureChange(ax, glacier, measure, \
    spring=False, summer=False, autumn=False, winter=False):
    check_measure(measure)
    
    dates = glacier.extract('date')
    measures = glacier.extract(measure)
    measures_change = measures.diff()
    seasons = glacier.extract('season')
    name = getGlacierName(glacier)
    ax.plot(dates, measures_change, color='darkgray')
    
    if spring:
        # dates_spr, measures_spr = getSeasonMeasures(
        #     dates, measures, seasons, 'SPR')
        # dates_spr_change = dates_spr[dates_spr.notna()]
        # measures_spr_change = measures_spr[measures_spr.notna()].diff()
        # ax.bar(dates_spr_change, measures_spr_change, \
        #     width=75, color=season_c['spring'], label='Spring')
        dates_spr_change = dates[seasons == 'SPR']
        measures_spr_change = measures_change[seasons == 'SPR']
        ax.scatter(dates_spr_change, measures_spr_change, color=season_c['spring'], label='Spring')
    if summer:
        # dates_sum, measures_sum = getSeasonMeasures(
        #     dates, measures, seasons, 'SUM')
        # dates_sum_change = dates_sum[dates_sum.notna()]
        # measures_sum_change = measures_sum[measures_sum.notna()].diff()
        # ax.bar(dates_sum_change, measures_sum_change, \
        #     width=75, color=season_c['summer'], label='Summer')
        dates_sum_change = dates[seasons == 'SUM']
        measures_sum_change = measures_change[seasons == 'SUM']
        ax.scatter(dates_sum_change, measures_sum_change, color=season_c
        ['summer'], label='Summer')
    if autumn:
        # dates_aut, measures_aut = getSeasonMeasures(
        #     dates, measures, seasons, 'AUT')
        # dates_aut_change = dates_aut[dates_aut.notna()]
        # measures_aut_change = measures_aut[measures_aut.notna()].diff()
        # ax.bar(dates_aut_change, measures_aut_change, \
        #     width=75, color=season_c['autumn'], label='Autumn')
        dates_aut_change = dates[seasons == 'AUT']
        measures_aut_change = measures_change[seasons == 'AUT']
        ax.scatter(dates_aut_change, measures_aut_change, color=season_c['autumn'], label='Autumn')
    if winter:
        # dates_win, measures_win = getSeasonMeasures(
        #     dates, measures, seasons, 'WIN')
        # dates_win_change = dates_win[dates_win.notna()]
        # measures_win_change = measures_win[measures_win.notna()].diff()
        # ax.bar(dates_win_change, measures_win_change, \
        #     width=75, color=season_c['winter'], label='Winter')
        dates_win_change = dates[seasons == 'WIN']
        measures_win_change = measures_change[seasons == 'WIN']
        ax.scatter(dates_win_change, measures_win_change, color=season_c
        ['winter'], label='Winter')
    
    ax.set_title('{}: Observed Seasonal {} Change'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('{} change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    zeroLine(ax)
    ax.legend()


def lengthDiff(ax, glacier, YEARS):
    name = getGlacierName(glacier)

    all_lengths = gm.fullTimeSeries(glacier, YEARS)

    # plot length and length-difference over time with NaN holes for no observation
    # -- plot length
    ax[0].plot(all_lengths.date, all_lengths.length, color='darkgray')
    ax[0].scatter(all_lengths[all_lengths.season == 'SPR'].date, all_lengths[all_lengths.season == 'SPR'].length, color=season_c['spring'], label='Spring')
    ax[0].scatter(all_lengths[all_lengths.season == 'AUT'].date, all_lengths[all_lengths.season == 'AUT'].length, color=season_c['autumn'], label='Autumn')
    ax[0].set_title(f'{name}: Length')
    ax[0].legend()
    # -- plot diff
    ax[1].plot(all_lengths.date, all_lengths.ldiff, color='darkgray')
    ax[1].scatter(all_lengths[all_lengths.season == 'SPR'].date, all_lengths[all_lengths.season == 'SPR'].ldiff, color=season_c['spring'], label='Spring')
    ax[1].scatter(all_lengths[all_lengths.season == 'AUT'].date, all_lengths[all_lengths.season == 'AUT'].ldiff, color=season_c['autumn'], label='Autumn')
    ax[1].set_title(f'{name}: Diff between lengths')
    zeroLine(ax[1])


def meanChangeRate(ax, glacier, measure, rolling_years):
    check_measure(measure)

    name = getGlacierName(glacier)

    measures = glacier.extract(measure)
    dates = glacier.dates

    df = pd.DataFrame(index=dates, data={'measure': measures.values})
    intermeasure_rate = df.measure.diff() / (dates.diff() / np.timedelta64(1, 'Y')).values
    rolling_rate = df.diff().rolling(str(rolling_years*365)+'D').mean()

    ax.plot(dates, intermeasure_rate, label='Rate of change')
    ax.plot(dates, rolling_rate, label='{}-year rolling mean'.format(rolling_years))
    ax.set_title('{}: Rate of {} Change'.format(name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('{} change ({}/yr)'.format(
        measure.capitalize(), measure_units[measure]))
    ax.legend()


def glacierObservations(ax, all_glaciers):
    yticklabs = []
    for g in all_glaciers:
        glacier = all_glaciers[g]
        name = getGlacierName(glacier)
        yticklabs.append(name)
        glacierid = glacier.extract('gid')
        obs_dates = glacier.extract('date')
        seasons = glacier.extract('season')
        # ax.plot(obs_dates, glacierid, '.', color='mediumblue')
        ax.plot(obs_dates[seasons=='AUT'], glacierid[seasons=='AUT'], '.', color='darkorange', label='Autumn')
        ax.plot(obs_dates[seasons=='SPR'], glacierid[seasons=='SPR'], '.', color='mediumseagreen', label='Spring')
    ax.set_title('Observation Time Series for Each Glacier')
    ax.set_xlabel('Date')
    ax.set_ylabel('Glacier')
    ax.legend(handles=ax.get_lines()[0:2], ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.18))
    yticklocs = range(1, len(all_glaciers)+1)
    plt.yticks(yticklocs, yticklabs)
    plt.ylim(bottom=0)


def measureSummary(ax, all_glaciers, measure, group, subtitle=None):
    check_measure(measure)

    for g in all_glaciers:
        if g in group:
            glacier = all_glaciers[g]
            dates = glacier.extract('date')
            cumul_measures, _, _ = gm.netMeasureChange(glacier, measure)
            name = getGlacierName(glacier)
            graph, = ax.plot(dates, cumul_measures,
                glacier_design[g]['s'], color=glacier_design[g]['c'],
                label=name)
            # figureProperties(fig, ax, graph)
        
    ax.set_title('Summary of observed {} changes{}'.format(
        measure, subtitle.title()))
    ax.set_xlabel('Date')
    ax.set_ylabel('Cumulative {} change ({})'.format(
        measure, measure_units[measure]))
    ax.legend(loc='center left', bbox_to_anchor=(1.04, 0.5))
    zeroLine(ax)


def normMeasureSummary(ax, all_glaciers, measure, group, subtitle=None):
    check_measure(measure)

    for g in all_glaciers:
        if g in group:
            glacier = all_glaciers[g]
            dates = glacier.extract('date')
            scaled_measure = gm.normMeasureChange(glacier, measure)
            name = getGlacierName(glacier)
            graph, = ax.plot(dates, scaled_measure,
                glacier_design[g]['s'], color=glacier_design[g]['c'],
                label=name)
        
    ax.set_title('Summary of Observed {} Changes, Normalized{}'.format(
        measure.capitalize(), subtitle.title()))
    ax.set_xlabel('Date')
    ax.set_ylabel('Scaled Cumulative {} Change'.format(measure.capitalize()))
    plt.ylim(-0.02, 1.02)
    ax.legend(loc='center left', bbox_to_anchor=(1.04, 0.5))
    zeroLine(ax)

