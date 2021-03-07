import numpy as np
import pandas as pd


def lin_trend(y, x):
    """
    function to calculate linear trend (y = mx + c) based on example numpy.linalg.lstsq
    :param y:
    :param x:
    :return: list with value of slope (m) amd intercept (c)
    """
    A = np.vstack([y, np.ones(len(y))]).T
    m, c = np.linalg.lstsq(A, x, rcond=None)[0]
    return m


def sst_trends():
    """
    function to calculate sst trends for different time ranges, scenarios and locations
    :return: dataframe with sst trends
    """
    locations = ['CAR', 'GBR', 'SEA']
    scenarios = ['RCP26', 'RCP45', 'RCP85']
    time_ranges = [['1970-01-01', '2010-12-31'],
                   ['1986-01-01', '2005-12-31'],
                   ['2081-01-01', '2100-12-31']]

    sst_trends_dtf = pd.DataFrame({'Time_range': [], 'Location': [], 'Scenario': [], 'SST_trend': []})

    for loc in locations:
        for sce in scenarios:
            mon = pd.read_csv('./Monthly-SST-scenarios-csv/Months-' + loc + '-' + sce + '-MPI.csv', names=['Month'],
                              header=None)
            sst = pd.read_csv('./Monthly-SST-scenarios-csv/SST-' + loc + '-' + sce + '-MPI.csv', names=['SST'],
                              header=None)

            timeseries_dtf = pd.DataFrame({'Month': mon.values.flatten(),
                                           'SST': sst.values.flatten(),
                                           'Date': pd.date_range("1955-01-01", "2100-12-31", freq="BM")})

            timeseries_dtf = timeseries_dtf.set_index('Date')

            for tr in time_ranges:
                trend = lin_trend(timeseries_dtf[tr[0]:tr[1]].SST, timeseries_dtf[tr[0]:tr[1]].Month)
                time_ranges_dtf = pd.DataFrame({'Time_range': [str(tr[0] + '_' + tr[1])],
                                                'Location': [loc],
                                                'Scenario': [sce],
                                                'SST_trend': [trend]})
                sst_trends_dtf = pd.concat([sst_trends_dtf, time_ranges_dtf])

    return sst_trends_dtf


def state_var_trends():
    """
    function to calculate coral, symbiont and trait trends for different time ranges, scenarios and locations
    :return: dataframe with trends of state variables
    """
    locations = ['CAR275', 'GBR920', 'SEA330']
    scenarios = ['RCP26', 'RCP45', 'RCP85']
    time_ranges = [['1970-01-01', '2010-12-31'],
                   ['1986-01-01', '2005-12-31'],
                   ['2081-01-01', '2100-12-31']]

    stvar_trends_dtf = pd.DataFrame({'Time_range': [], 'Location': [], 'Scenario': [],
                                     'Coral_trend': [], 'Symb_trend': [], 'Trait_trend': []})
    time_years = pd.read_csv('./Results-csv/Time_in_years.csv', names=['Time'], header=None)

    for loc in locations:
        for sce in scenarios:
            coral = pd.read_csv('./Results-csv/CORAL-' + sce + '-' + loc + '.csv', names=['Coral'], header=None)
            symbiont = pd.read_csv('./Results-csv/SYMB-' + sce + '-' + loc + '.csv', names=['Symbiont'], header=None)
            trait = pd.read_csv('./Results-csv/TRAIT-' + sce + '-' + loc + '.csv', names=['Trait'], header=None)

            timeseries_dtf = pd.DataFrame({'Time': time_years[2001 * 12:].values.flatten(),
                                           'Coral': coral[2001 * 12:].values.flatten(),
                                           'Symbiont': symbiont[2001 * 12:].values.flatten(),
                                           'Trait': trait[2001 * 12:].values.flatten(),
                                           'Date': pd.date_range("1955-01-01", "2100-12-31", freq="BM")})

            timeseries_dtf = timeseries_dtf.set_index('Date')

            for tr in time_ranges:
                trend_coral = lin_trend(timeseries_dtf[tr[0]:tr[1]].Coral, timeseries_dtf[tr[0]:tr[1]].Time)
                trend_symbiont = lin_trend(timeseries_dtf[tr[0]:tr[1]].Symbiont, timeseries_dtf[tr[0]:tr[1]].Time)
                trend_trait = lin_trend(timeseries_dtf[tr[0]:tr[1]].Trait, timeseries_dtf[tr[0]:tr[1]].Time)
                time_ranges_dtf = pd.DataFrame({'Time_range': [str(tr[0] + '_' + tr[1])],
                                                'Location': [loc[:-3]],
                                                'Scenario': [sce],
                                                'Coral_trend': [trend_coral],
                                                'Symb_trend': [trend_symbiont],
                                                'Trait_trend': [trend_trait]})
                stvar_trends_dtf = pd.concat([stvar_trends_dtf, time_ranges_dtf])

    return stvar_trends_dtf


sst_dtf = sst_trends()
stvar_dtf = state_var_trends()
dtf_all = sst_dtf.merge(stvar_dtf)
dtf_all = pd.DataFrame(dtf_all, columns=['Time_range', 'Location', 'Scenario',
                                         'SST_trend', 'Coral_trend', 'Symb_trend', 'Trait_trend'])
dtf_all.to_latex()
dtf_all.to_csv('trends.csv')
